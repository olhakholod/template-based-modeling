###############################################################################
# Libraries/imports
###############################################################################
import os
import numpy as np
import math

###############################################################################
# Global variables
###############################################################################
ALIGNMENT_PATH = "../data/alignment/aligned_sequences.ali"
TEMP_PDB_PATH  = "../data/template_4i1a.pdb"
TARG_PDB_PATH  = "../test/target_5z82.pdb"

d_mean = 2.1939
d_stdd = 0.3829

###############################################################################
# Functions
###############################################################################
def check_data(template_sequence,c_temp):
	case1,case2,case3,case4,case5 = 0,0,0,0,0
	i1,i2,i4 = [],[],[]
	for i in range(len(template_sequence)):
		if target_seq[i]==template_sequence[i] and template_sequence[i]!='-':
			case1 += 1
			i1.append(i)
			continue
		if template_sequence[i]!='-':
			if target_seq[i] == '-':
				case2 += 1
				i2.append(i)
			else:
				case4 += 1
				i4.append(i)
		else:
			if target_seq[i]!='-':
				case3 += 1
			else:
				case5 += 1
	print(case1,case2,case3,case4,case5)

	copy = 0
	for i in range(len(i4)):
		if c_temp[i4[i]][0]!=math.inf:
			copy += 1
	print(copy)

	l = i1 + i2 + i4
	corr = 0
	incorr = 0
	for i in l:
		if c_temp[i][0]!=math.inf:
			corr += 1
		else:
			incorr += 1
	print(corr,incorr)

def load_sequences_from_file(alignment_file):
	'''
	Load a list of fasta sequences from a given file with FASTA format.
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
	'''
	Load carbon-alpha coordinates in the specified PDB file, and store them in
	an array equal in length to the aligned sequences. If no amino acid is 
	specified in the sequence (loop) store math.inf as coordinate and 'UNK' as
	amino-acid label. Return the arrays of amino acids labels, the coordinates 
	c (x,y,z) and the temperature factor.
	'''
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
			new_x.append(math.inf)
			new_y.append(math.inf)
			new_z.append(math.inf)
			new_t.append(0)
	c = np.array([[new_x[i],new_y[i],new_z[i]] for i in range(len(new_x))])
	return c

def distance(point1, point2):
	'''
	Calculate euclidean distance using numpy vectors 1 and 2 for
	x,y,z coordinates.
	'''
	return np.sqrt(np.sum((point1 - point2)**2, axis=0))

def distance_explicit(x1,x2,y1,y2,z1,z2):
	'''
	Calculate euclidian distance explicitly.
	'''
	return math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))

def distance_norm(point1,point2):
	'''
	Calculate euclidean distance using numpy vector norm.
	'''
	return np.linalg.norm(point1-point2)

def distance_matrix(x,y,z): #this is a silly way to do it, but it works...
	'''
	distance_matrix calculates a distance matrix for each residue 
	pair i and j in a sequence using three arrays x,y,z containing the 
	coordinates.
	'''
	n = len(x)
	D = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if x[i]!=math.inf and x[j]!=math.inf:
				D[i][j] = distance_explicit(x[i],x[j],y[i],y[j],z[i],z[j])
			else:
				D[i][j] = 150.0
	return D

def distance_matrix_vectorized(c):
	'''
	calculates a distance matrix for each residue pair i and j in a 
	sequence using a single array c of 3-dimensional (x,y,z coordinates) 
	vectors.
	'''
	n = len(c)
	D = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if c[i][0]!=math.inf and c[j][0]!=math.inf:
				D[i][j] = distance(c[i],c[j])
			else:
				D[i][j] = 150.0
	return D	

def next_coordinate(coordinate_i):
	count = 0
	while True:
		r = np.random.normal(0,1,3)
		a = np.sqrt(np.sum(np.power(r,2)))
		c = coordinate_i + np.random.normal(3.8,0.65,1)*(r/a)
		d = distance(coordinate_i,c)
		if d > 2.5 and d < 5.1:
			break
		count += 1
		if count > 100:
			print("100 iterations, returning any random number.")
			break
	return c

def noisy_coordinate(coordinate_i):
	for i in range(10):
		r = np.random.normal(0,1,3)
		a = np.sqrt(np.sum(np.power(r,2)))
		c = coordinate_i + np.random.normal(0.25,0.1,1)*(r/a)
	return c

def tm_score(p_c,n_c):
    '''
    Calculate TM-Score (template-modeling score) between predicted 
    target structure and real target structure. Takes as arguments 
    the coordinates of two structures p_c (predicted target coordinates) and 
    n_c (native/real target coordinates), iteratively calculates scores for 
    each aligned region, and returns the maximum score.
    '''
    if len(p_c) != len(n_c):
        print("\nERROR: Predicted and native coordinate arrays have different lengths.")
        return None        
    # ----- COLLECT INDEXES OF ALIGNED REGIONS FROM SEQ -----------------------
    names,sequences = load_sequences_from_file(ALIGNMENT_PATH) #aligned seqs of targ and temp
    aligned_target  = sequences[0] #get the aligned sequence for target (first one)
    regions         = []
    aligned_region  = [] #temporarily store regions here
    flag            = 0 #flag to tell if it's currently reading a region or not
    for i in range(len(aligned_target)): # ----- COLLECT INDEXES OF ALIGNED REGIONS FROM SEQ
        if aligned_target[i]!='-': # NOT A LOOP
            if flag==0: # new region, begin reading and change flag
                aligned_region = [i]
                flag = 1
            else: # already reading a region, keep going...
                aligned_region.append(i)
        else: # LOOP
            if flag==0: # not currently reading, skip
                continue 
            else: # currently reading, append indexes and change flag to skip loop
                regions.append(aligned_region) #append indexes of region
                flag = 0
    regions.append(aligned_region) #append indexes of last region
    # ----- CALCULATE SCORES OF COORDINATES IN ALIGNED REGIONS ----------------
    l_target = 276 #this is terrible...fix to len(target_sequence) for any seq? ...nah.
    d_0      = (1.24*np.cbrt(l_target-15))-1.8 #d_0 is calculated like this..god knows why..
    scores   = [] #store the scores of the regions here to pick a max after
    for i in regions: # ----- CALCULATE SCORES OF COORDINATES IN ALIGNED REGIONS
        l_aligned = []
        for j in i:
            d_i = distance(p_c[j],n_c[j])
            l_aligned.append((1.0/(1.0+(d_i/d_0)**2)))
        this_score = (1.0/l_target) * np.sum(l_aligned)
        scores.append(this_score)
    # ----- RETURN MAX SCORE --------------------------------------------------
    return max(scores)

def rmsd(p_c,n_c):
    '''
    Calculate the root-mean-square-deviation between the predicted target 
    structure (p_c) and the real/native target structure (n_c). Square norm is 
    distance, so the distance() function is used for the distance between the 
    ith-residue in predicted and native, returns the square-root of the mean 
    of these distances.
    '''
    if len(p_c) != len(n_c):
        print("\nERROR: Predicted and native coordinate arrays have different lengths.")
        return None
    return np.sqrt(np.mean([distance(i,j) for i,j in zip(p_c,n_c)]))

def rmsd_without_loops(p_c,n_c):
    '''
    Calculate the root-mean-square-deviation between the predicted target 
    structure (p_c) and the real/native target structure (n_c) ONLY FOR THE ALIGNED 
    POSITIONS (i.e.: skipping loops). Square norm is distance, so the 
    distance() function is used for the distance between the ith-residue in predicted 
    and native, returns the square-root of the mean of these distances.
    '''
    if len(p_c) != len(n_c):
        print("\nERROR: Predicted and native coordinate arrays have different lengths.")
        return None
    names,sequences = load_sequences_from_file(ALIGNMENT_PATH) #aligned seqs of targ and temp
    aligned_target  = sequences[0] #get the aligned sequence for target (first one)
    if(len(p_c) != len(aligned_target)):
        print("\nERROR: Coordinates and sequence arrays have different lengths.")
        return None
    # All good up to here, now gather distances for all non-loop coordinates
    d = [distance(p_c[i],n_c[i]) for i in range(len(aligned_target)) if aligned_target[i]!='-']
    # return the square root of the mean distance
    return np.sqrt(np.mean(d))

def initialize_target_coordinates(s_temp,c_temp):
	c_targ = np.zeros((len(c_temp),3))
	for i in range(c_temp.shape[0]):
		if c_temp[i][0]!=math.inf:
			c_targ[i] = c_temp[i]
	c_targ = initialize_loops(s_temp,c_temp,c_targ)
	return c_targ

def initialize_loops(s_temp,c_temp,c_targ):
	out   = c_targ
	loops = []
	loop  = []
	flag  = 0 #flag 
	for i in range(len(s_temp)): # ----- COLLECT INDEXES LOOPS FROM SEQ
		if s_temp[i]=='-': # LOOP
			if flag==0: # new loop begin reading and change flag
				loop = [i]
				flag = 1
			else: # already reading, keep going...
				loop.append(i)
		else:
			if flag==0: # not currently reading, skip
				continue 
			else: # currently reading, append indexes and change flag
				loops.append(loop) #append indexes of region
				flag = 0
	loops.append(loop)
	for l in loops:
		if l[0]==0:
			for i in range(len(l)-1,-1,-1): #reverse
				out[l[i]] = next_coordinate(c_targ[l[i]+1])
		else:
			for i in l: #fwd derivation of next coordinates at ~3.8A
				out[i] = next_coordinate(c_targ[i-1])
	return out

def gradient(i,j,d_hat,std_dev):
	d = distance(i,j) 
	return -(1/std_dev)(d-d_hat)*((i - j)/d)

def gradient_descent(i,j,template_coords,D):
	x = x + LEARNING_RATE * gradient(template_coords[i],template_coords[j])

def neighbor_distances(D):
	d_1 = []
	for i in range(len(D)-1):
		d_1.append(D[i][i+1])
	return np.array(d_1)

###############################################################################
# Main ('python3 preprocessing.py')
###############################################################################
def main():
	np.random.seed(123)
	names,sequences   = load_sequences_from_file(ALIGNMENT_PATH)
	target_sequence   = sequences[0]
	template_sequence = sequences[1]
	template_coords   = load_coordinates_from_file(TEMP_PDB_PATH,template_sequence)

	# print(len([i for i in c_temp if i[0]!=math.inf]))
	template_distance_matrix = distance_matrix_vectorized(template_coords)
	template_distance_ca1    = neighbor_distances(template_distance_matrix)

	print(template_distance_matrix[0:5][0:5])
	print(neighbor_distances(template_distance_matrix[0:5][0:5]))

	# o = initialize_target_coordinates(template_sequence,template_coords)
	# print(o[-10])
	# print(template_coords[-10])

if __name__ == '__main__':
	main()