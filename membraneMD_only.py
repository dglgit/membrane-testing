from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology as FFTopology
from openmmforcefields.generators import GAFFTemplateGenerator
import sys
import random
import numpy as np
pathRoot="empty1"
# Load the membrane and alcohol system
print("Loading initial structure...")
#pdb = PDBFile('empty.pdb')  # Input file with initial membrane structure
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml','amber14/lipid17.xml')

# Add alcohol molecules (e.g., ethanol)
ethanol = Molecule.from_smiles('CCO')  # Ethanol SMILES
ethanol.generate_conformers(n_conformers=1)
ethanol.assign_partial_charges(partial_charge_method='gasteiger')
ethanol_topology = ethanol.to_topology().to_openmm()
ethanol_positions = ethanol.conformers[0].to_openmm()


# Add GAFF parameters for ethanol
gaff_template_generator = GAFFTemplateGenerator(molecules=[ethanol])
forcefield.registerTemplateGenerator(gaff_template_generator.generator)

# Build the system with a POPC membrane
print("Building the system...")
modeller = Modeller(ethanol_topology, ethanol_positions)
modeller.topology.setUnitCellDimensions([3,3,15]*nanometers)
print("System dimensions:", modeller.topology.getUnitCellDimensions())
modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1*nanometer, membraneCenterZ=0.0*nanometer)
modeller.deleteWater()
positions=modeller.positions
def get2DRad(v):
    return float(v[0]**2+v[1]**2)**0.5
for residue in modeller.topology.residues():
    if residue.name =='UNK':
        print('deleting first ethanol')
        modeller.delete([residue])
    if residue.name in ['POPC', 'POP']:  # Membrane lipid names
        for atom in residue.atoms():
            #print(positions)
            #print(atom.index)
            #print(positions[atom.index])
            pos=positions[atom.index]
            if abs(pos[0])>1*nanometer or abs(pos[1])>1*nanometer:
                modeller.delete([residue])
                break
# Create the system
print("Creating the system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
)

# Add a barostat for pressure control (NPT ensemble)
system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))

# Set up integrator and simulation
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)

# Set initial positions and velocities
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

# Minimize energy
print("Minimizing energy...")
simulation.minimizeEnergy(tolerance=3)

# Equilibrate the system
with open("init_structure_"+pathRoot+'.pdb','w') as pdbout:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=pdbout, keepIds=True)
print("Equilibrating...")

simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(DCDReporter(f'trajectory_{pathRoot}.dcd', 1000, enforcePeriodicBox=True))
simulation.step(50000)  # Equilibration for 100 ps

# Production simulation
print("Running production simulation...")
simulation.step(1000000)  # Run for ~1 ns

# Save final structure
with open(f"final_structure_{pathRoot}.pdb", "w") as output:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), output)

print("Simulation complete!")

