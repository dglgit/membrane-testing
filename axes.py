from pymol.cgo import CYLINDER, cyl_text
from pymol.vfont import plain

obj = [CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0., CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,]

cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])

cmd.load_cgo(obj,'axes')
