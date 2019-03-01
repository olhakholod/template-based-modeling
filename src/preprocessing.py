# Libraries/imports
import os
import numpy as np

# Global variables
alignment_dir  = "../test/"
alignment_file = "T0951.ali"


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
			new_x.append(0)
			new_y.append(0)
			new_z.append(0)
			new_t.append(0)
	return new_a,new_x,new_y,new_z,new_t

def distances(x,y,z):
	d_x   = []
	for i in range(len(x)-1):
		row_x = []
		for j in range(i+1,len(x)):
			row_x.append(x[i] - x[j])
		d_x.append(row_x)

	print(d_x)
	print(len(x),len(d_x))


# Main (this runs when calling 'python3 preprocessing.py')
def main():
	# print("Hello world! This is main.")
	n,s       = load_sequences_from_file(alignment_dir+alignment_file)
	a,x,y,z,t = load_coordinates_from_file('../data/target_native_structure_5z82.pdb',s[0])
	distances(x,y,z)

if __name__ == '__main__':
	main()