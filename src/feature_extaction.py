from pylab import *
from prody import *
import Bio.PDB as pdb
import numpy as np;

backbone = parsePDB('4i1a.pdb', subset='bb');
coordinates = backbone.getCoords()

names = backbone.getNames();
numAngles = int(names.size/4);
caDistance =[];
phi =[];
psi =[];
for i in range(0, numAngles-2):
#	nIndex  = 4*i;
	precaIndex = 4*i+1;
	cpIndex  = 4*i+2;
#	oIndex  = 4*i+3
	nIndex = 4*i+4;
	caIndex =4*i+5;
	cp2Index =4*i+6;
#       o2Index =4*i+7;
	n2Index =4*i+8;
	
	vcp    = pdb.Vector(coordinates[cpIndex, 0:3]);
	vn   = pdb.Vector(coordinates[nIndex, 0:3]);
	vca   = pdb.Vector(coordinates[caIndex, 0:3]);
	vcp2   = pdb.Vector(coordinates[cp2Index, 0:3]);
	vn2  = pdb.Vector(coordinates[n2Index, 0:3]);
	
	alpha1 = coordinates[precaIndex, 0:3];
	alpha2 = coordinates[caIndex, 0:3];
	
	currentDistance = np.linalg.norm(alpha1-alpha2);

	currentPhi = pdb.calc_dihedral(vcp, vn, vca, vcp2);
	currentPsi = pdb.calc_dihedral(vn, vca, vcp2, vn2);
	phi.append(currentPhi);
	psi.append(currentPsi);
	caDistance.append(currentDistance);
print(len(phi));
print(caDistance);
