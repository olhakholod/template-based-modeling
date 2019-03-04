from pylab import *
from prody import *
import Bio.PDB as pdb
import numpy as np;
import matplotlib.pyplot as plt

PDB_FILE='4i1a.pdb';
def histogram(pdb_file):
	backbone = parsePDB(pdb_file, subset='bb');
	coordinates = backbone.getCoords()
	names = backbone.getNames();
	numCA = int(names.size/4);
	caDistances = np.empty([numCA, numCA]);
	for i in range(0, numCA):
		for j in range(0, numCA):
			indexI= 4*i+1
			indexJ= 4*j+1;
			CAi = coordinates[indexI, 0:3];
			CAj = coordinates[indexJ, 0:3];
			currentDistance = np.linalg.norm(CAi-CAj);
			caDistances[i][j] = currentDistance;

	print(caDistances[1:10,1:10]);

	return caDistances;