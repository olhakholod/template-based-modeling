import numpy as np
import math
import Bio.PDB as pdb

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def distance(point1, point2):
    squareDistance = np.sum((point1 - point2)**2, axis=0)
    sqrtDistance = np.sqrt(squareDistance)
    return sqrtDistance

def pdf(tarDistance, temDistance, sigma):
    prob = 1/(sigma * np.sqrt(2*math.pi)) * np.exp(-1/2*((tarDistance - temDistance)/sigma)**2)
    return prob

def objProb(tarPoints, temPoints, temW, sigma):
    nTarPts = tarPoints.shape[0]
    nTemplates = temPoints.shape[0]
    nTemPts = temPoints.shape[1]

    tarDisArray = np.zeros((nTarPts, nTarPts))
    temDisArray = np.zeros((nTemplates, nTemPts, nTemPts))
    cacheArray = np.zeros((nTemplates, nTemPts, nTemPts))

    pdfArray = np.zeros((nTemplates, nTarPts, nTarPts))

    # compute distances between target and template atoms
    for i in range(0,nTarPts-1):
        for j in range(i+1,nTarPts):
            tarDisArray[i,j] = distance(tarPoints[i,:], tarPoints[j,:])
            
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
    gradFpoints = np.zeros((nTar,3))

    gradFd = (1/sumPdf) * sumCache

    for i in range(0,nTar):
        if (i < nTar-1):
            gradFpoints[i] = gradFd[i,i+1] * (tarPoints[i] - tarPoints[i+1])/tarDisArray[i,i+1]
        else:
            gradFpoints[i] = gradFd[1,i] * (tarPoints[i] - tarPoints[1])/tarDisArray[1,i]
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
    nseed = 1
    np.random.seed(nseed)
    tarFile = 'target_T0951.fasta'
    temFile = '4i1a.pdb'
    
    temPoints0 = getTemplate(temFile)
    tem1 = temPoints0[0:8, :]
    tem2 = temPoints0[10:28, :]
    tem3 = temPoints0[30:58, :]
    temPoints = np.vstack((tem1, tem2, tem3))

    tarPoints = temPoints - np.random.rand(temPoints.shape[0], temPoints.shape[1])*10
    temPoints = np.reshape(temPoints, (1,temPoints.shape[0], temPoints.shape[1]))

    temW = [1.0]

    sigma = 0.5

    optimalTarget = gradDescent(tarPoints, temPoints, temW, sigma, alpha=0.01, tolerance=10**(-5), maxiter=1000)
    print('optimalTarget: \n', optimalTarget)

if __name__ == '__main__':
    main()

