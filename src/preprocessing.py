# Libraries/imports
import os
import numpy as np
import math

# Global variables

ALIGNMENT_PATH = "../data/alignment/aligned_sequences.ali"
TEMP_PDB_PATH  = "../data/template_4i1a.pdb"
TARG_PDB_PATH  = "../test/target_5z82.pdb"

# Functions
def load_sequences_from_file(alignment_file):
	'''
	This function loads a fasta sequence from a .fa file.
	'''
	names     = []
	sequences = []
	count     = 0
	skip      = 0
	fp = open(alignment_file,'r')
	for line in fp:
		if line[0] == '>':
			names.append(line.split(';')[1].rstrip())
			if count > 0:
				sequences.append(s)
			count += 1
			s = ""
		else:
			if line=='\n' or line.rstrip()[0:9]=="structure":
				continue
			else:
				s += line.rstrip()[:-1]
	sequences.append(s)
	fp.close()
	return names,sequences

def load_coordinates_from_file(file,sequence):
	a,x,y,z,t = [],[],[],[],[]
	fp = open(file,'r')
	for line in fp:
		if line[0:4] == 'ATOM':
			if line[12:16].strip()=='CA':
				a.append(line[17:20])
				x.append(float(line[30:38]))
				y.append(float(line[38:46]))
				z.append(float(line[46:54]))
				t.append(float(line[60:66]))
	fp.close()
	# print(distance(max(x),min(x),max(y),min(y),max(z),min(z)))
	new_a,new_x,new_y,new_z,new_t = [],[],[],[],[]
	j = 0
	for i in sequence:
		if i!='-':
			new_a.append(a[j])
			new_x.append(x[j])
			new_y.append(y[j])
			new_z.append(z[j])
			new_t.append(t[j])
			j += 1
		else:
			new_a.append('UNK')
			new_x.append(math.inf)
			new_y.append(math.inf)
			new_z.append(math.inf)
			new_t.append(0)
	return new_a,new_x,new_y,new_z,new_t

def distance(x1,x2,y1,y2,z1,z2):
	'''
	Calculate euclidian distance explicitly.
	'''
	return math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))

def np_distance(a,b):
	'''
	Calculate euclidean distance using numpy vector norm.
	'''
	return np.linalg.norm(a-b)

def distance_matrix(x,y,z): #this is a silly way to do it, but it works...
	'''
	distance_matrix calculates a distance matrix for each residue 
	pair i and j in a sequence
	'''
	n = len(x)
	D = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if x[i]!=math.inf and x[j]!=math.inf:
				D[i][j] = distance(x[i],x[j],y[i],y[j],z[i],z[j])
			else:
				D[i][j] = 150.0
	return D

def tm_score():
    names,sequences = load_sequences_from_file(ALIGNMENT_PATH) #load aligned seqs of targ and temp
    aligned_target  = sequences[0] #get the aligned sequence for target (first one)
    l_target        = 276 #fix to len(target_sequence) for any seq? ...nah.
    d_0             = (1.24*np.cbrt(l_target-15))-1.8 #d_0 is calculated like this..god knows why..
    scores          = [] #store the scores of the regions here to pick a max after
    regions         = []
    aligned_region  = [] #temporarily store regions here
    flag            = 0 #flag to tell if it's already reading a region or not
    for i in range(len(aligned_target)):
        if aligned_target[i]!='-':
            if flag==0: # new region, begin reading
                aligned_region = [i]
                flag = 1
            else: # already reading a region keep going...
                aligned_region.append(i)
        else:
            if flag==0:
                continue 
            else:
                regions.append(aligned_region) #append indexes of region
                flag = 0
    regions.append(aligned_region) #append last indexes of region

    for i in regions:
    	# d_i = distance()
 


# Main (this runs when calling 'python3 preprocessing.py')
def main():
	n,s       = load_sequences_from_file(ALIGNMENT_PATH)
	a,x,y,z,t = load_coordinates_from_file(TEMP_PDB_PATH,s[1])
	d = distance_matrix(x,y,z)
	tm_score()

if __name__ == '__main__':
	main()