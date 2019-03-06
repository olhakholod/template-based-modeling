import numpy as np
import math
import Bio.PDB as pdb
import matplotlib as mpl
mpl.use("TkAgg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import preprocessing as pp

def distance(point1, point2):
    return np.sqrt(np.sum((point1 - point2)**2, axis=0))

def pdf(tarDistance, temDistance, sigma):
    prob = 1/(sigma * np.sqrt(2*math.pi)) * np.exp(-1/2*((tarDistance - temDistance)/sigma)**2)
    return prob

def sample_normal(mean,std_dev):
    return np.random.normal(mean,std_dev,30)

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
    # ----- COLLECT INDEXES OF ALIGNED REGIONS FROM SEQ -------------------------------------------
    names,sequences = load_sequences_from_file(ALIGNMENT_PATH) #load aligned seqs of targ and temp
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
    # ----- CALCULATE SCORES OF COORDINATES IN ALIGNED REGIONS ------------------------------------
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
    # ----- RETURN MAX SCORE ----------------------------------------------------------------------
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
    names,sequences = load_sequences_from_file(ALIGNMENT_PATH) #load aligned seqs of targ and temp
    aligned_target  = sequences[0] #get the aligned sequence for target (first one)
    if(len(p_c) != len(aligned_target)):
        print("\nERROR: Coordinates and sequence arrays have different lengths.")
        return None
    # All good up to here, now gather distances for all non-loop coordinates
    d = [distance(p_c[i],n_c[i]) for i in range(len(aligned_target)) if aligned_target[i]!='-']
    # return the square root of the mean distance
    return np.sqrt(np.mean(d))

def objProb(tarPoints, temPoints, temW, sigma):
    # initialization
    nTarPts    = tarPoints.shape[0]
    nTemplates = temPoints.shape[0]
    nTemPts    = temPoints.shape[1]

    tarDisArray = np.zeros((nTarPts, nTarPts))
    temDisArray = np.zeros((nTemplates, nTemPts, nTemPts))
    cacheArray  = np.zeros((nTemplates, nTemPts, nTemPts))

    pdfArray = np.zeros((nTemplates, nTarPts, nTarPts))

    # compute distances for target and template atoms
    # compute distances for target atoms 
    for i in range(0,nTarPts-1):
        for j in range(i+1,nTarPts):
            tarDisArray[i,j] = distance(tarPoints[i,:], tarPoints[j,:])
    
    # compute distances for template atoms         
    for k in range(0,nTemplates):
        for i in range(0,nTemPts-1):
            for j in range(i+1,nTemPts):
                temDisArray[k,i,j] = distance(temPoints[k,i,:], temPoints[k,j,:])
                
                
    # compute pdf
    for k in range(0,nTemplates):
        for i in range(0,nTarPts-1):
            for j in range(i+1,nTarPts):
                pdfArray[k,i,j] = temW[k] * pdf(tarDisArray[i,j], temDisArray[k,i,j], sigma)
                cacheArray[k,i,j] = pdfArray[k,i,j] * (tarDisArray[i,j] - temDisArray[k,i,j]) / (sigma**2) 

    # sum pdf for each distance:
    sumPdf = np.sum(pdfArray, axis=0)
    sumCache = np.sum(cacheArray, axis=0) 

    # multiply over half of sumPdf (not include diagonal):
    upperPdf_noDiag = np.triu(sumPdf, k=1)
    productPdf = np.prod(upperPdf_noDiag[upperPdf_noDiag > 0])

    # maximize value with log function
    logArray = -np.log(upperPdf_noDiag[upperPdf_noDiag > 0])
    logPdf = np.sum(logArray)
    print('logPdf: ', logPdf)

    return tarDisArray, temDisArray, pdfArray, sumPdf, sumCache, productPdf, logPdf

def gradient(tarPoints, tarDisArray, sumPdf, sumCache):
    nTar = tarDisArray.shape[0]

    gradFd = np.zeros(tarDisArray.shape)
    # gradFx = np.zeros((nTar,1))
    # gradFy = np.zeros((nTar,1))
    # gradFz = np.zeros((nTar,1))
    gradFpoints = np.zeros((nTar,3))

    gradFd = (1/sumPdf) * sumCache
    # print('gradFd: ')
    # print(gradFd)

    for i in range(0,nTar):
        if (i < nTar-1):
            # gradFx[i] = gradFd[i,i+1] * (tarPoints[i,0] - tarPoints[i+1,0])/tarDisArray[i,i+1]
            # gradFy[i] = gradFd[i,i+1] * (tarPoints[i,1] - tarPoints[i+1,1])/tarDisArray[i,i+1]
            # gradFz[i] = gradFd[i,i+1] * (tarPoints[i,2] - tarPoints[i+1,2])/tarDisArray[i,i+1]
            gradFpoints[i] = gradFd[i,i+1] * (tarPoints[i] - tarPoints[i+1])/tarDisArray[i,i+1]
        else:
            # gradFx[i] = gradFd[1,i] * (tarPoints[i,0] - tarPoints[1,0])/tarDisArray[1,i]
            # gradFy[i] = gradFd[1,i] * (tarPoints[i,1] - tarPoints[1,1])/tarDisArray[1,i]
            # gradFz[i] = gradFd[1,i] * (tarPoints[i,2] - tarPoints[1,2])/tarDisArray[1,i]
            gradFpoints[i] = gradFd[1,i] * (tarPoints[i] - tarPoints[1])/tarDisArray[1,i]
    # gradFpoints = np.hstack((gradFx, gradFy, gradFz))
    return gradFd, gradFpoints


def gradDescent(tarPoints0, temPoints, temW, sigma, alpha=0.05, tolerance=10**(-5), maxiter=200):

    tarPoints = tarPoints0
    tarDisArray, temDisArray, pdfArray, sumPdf, sumCache, productPdf, logPdf = objProb(tarPoints, temPoints, temW, sigma)
    iter = 0
    error = 0.5
    while (iter < maxiter) and (error > tolerance):

        gradFd, gradFpoints = gradient(tarPoints, tarDisArray, sumPdf, sumCache)
        tarPoints += - alpha*gradFpoints

        new_tarDisArray, new_temDisArray, new_pdfArray, new_sumPdf, new_sumCache, new_productPdf, new_logPdf = objProb(tarPoints, temPoints, temW, sigma)

        error = abs(logPdf - new_logPdf)

        if  error < tolerance:
            return tarPoints
        else:
            tarDisArray = new_tarDisArray
            sumPdf = new_sumPdf
            sumCache = new_sumCache
            logPdf = new_logPdf

        iter += 1
        print('iteration: ', iter, 'error: ', error)

    return tarPoints

def getTemplate(temFile):
    parser = pdb.PDBParser()
    struct = parser.get_structure('template', temFile)
    temPoints = list()
    for model in struct:
        for chain in model:
            for res in chain:
                for atom in res:
                   if (atom.name == 'CA'):
                        vector = atom.get_vector()
                        temPoints.append(list(vector))
    return np.array(temPoints)

def writePoints(inFile, outFile, optimalPoints):
    parser = pdb.PDBParser()
    struct = parser.get_structure('target', inFile)

    optimalPoints = np.array(optimalPoints)
    i = 0

    for model in struct:
        for chain in model:
            for res in chain:
                for atom in res:
                    if (atom.name == 'CA'):
                        atom.coord = optimalPoints[i,:]
                        i += 1

    io = pdb.PDBIO()
    io.set_structure(struct)
    io.save(outFile)
    return

def showProtein():

    return

def showFunction():

    return

def main():
    np.random.seed(1)
    tarFile = '../data/target_sequence.txt'
    temFile = '../data/template_4i1a.pdb'
    
    # select points from template that align with target (need to check)
    temPoints0 = getTemplate(temFile)
    tem1 = temPoints0[0:8, :]
    tem2 = temPoints0[10:28, :]
    tem3 = temPoints0[30:58, :]
    temPoints = np.vstack((tem1, tem2, tem3))

    tarPoints = temPoints - np.random.rand(temPoints.shape[0], temPoints.shape[1])*10
    temPoints = np.reshape(temPoints, (1,temPoints.shape[0], temPoints.shape[1]))
    # tarPoints = np.random.randint(15, 20, size=(temPoints.shape[1],3))/1.0

    # tarPoints = tarPoints.astype(float)
    # temPoints = temPoints.astype(float)

    # temW = [0.4, 0.6]
    temW = [1.0]

    sigma = 0.5

    optimalTarget = gradDescent(tarPoints, temPoints, temW, sigma, alpha=0.01, tolerance=10**(-5), maxiter=1000)
    print('optimalTarget: \n', optimalTarget)

    # visualize both target and template protein chains
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(temPoints[0,:,0], temPoints[0,:,1], temPoints[0,:,2], label='template')
    ax.plot(optimalTarget[:,0], optimalTarget[:,1], optimalTarget[:,2], label='target approximation')
    ax.legend(['Template','Target'])
    fig.savefig('target_T0951.png')
    plt.show()


if __name__ == '__main__':
    main()

