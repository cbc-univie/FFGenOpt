import string
import os
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment

qmfactor = 0.957

#source: from afmm_read_gaupy.py:
def read_charmm(filename):
    """Takes the name of the CHARMM output file that contains
    normal modes"""

    freq = []
    X = []
    Y = []
    Z = []

    vibmod = "  VIBRATION MODE"
    eigvec = "   EIGENVECTOR:"

    try:
        input = open(filename, 'r')
    except IOError:
        print("Cannot read file " + filename)
        sys.exit(-1)
    x = input.readline()
    while len(x) != 0:
        if x[:len(vibmod)] == vibmod:
            #l = string.split(x)
            l = x.split()
            if len(l[3]) > 12:
                if isinstance(l[3][11:], float):
                    freq.append(float(l[3][11:]))
                else:
                    freq.append(0)
            else:
                freq.append(float(l[4]))
            x = input.readline()
            while x[:len(eigvec)] != eigvec:
                x = input.readline()
            # get next line
            l = input.readline().split()
            while len(l) == 7:
                X.append(float(l[4]))
                Y.append(float(l[5]))
                Z.append(float(l[6]))
                l = input.readline().split()
        x = input.readline()
    input.close()
    
    # number of atoms
    N = int(len(freq) / 3)

    # return results, skipping the zero frequencies
    # this may cause problems with lonepairs as they are also zero and come before the vib freqs
    return freq[6:], X[N * 6:], Y[N * 6:], Z[N * 6:] #original one

#from afmm_read_gaupy.py
def read_gaussian(filename):
    """Takes the name of the Gaussian output file that contains
    normal modes"""
    freq = []
    X = []
    Y = []
    Z = []
    try:
        input = open(filename, 'r')
    except IOError:
        print("Cannot read file " + filename)
        sys.exit(-1)
    x = input.readline()
    while len(x) != 0:
        if x.find(" Frequencies") == 0:
            l = x.split()
            for i in range(3):
                freq.append(float(l[i + 2]))
            l = input.readline().split()
            while l[0][0] not in string.digits:
                l = input.readline().split()
            bx1 = []
            by1 = []
            bz1 = []
            bx2 = []
            by2 = []
            bz2 = []
            bx3 = []
            by3 = []
            bz3 = []
            while len(l) == 11:
                bx1.append(float(l[2]))
                by1.append(float(l[3]))
                bz1.append(float(l[4]))
                bx2.append(float(l[5]))
                by2.append(float(l[6]))
                bz2.append(float(l[7]))
                bx3.append(float(l[8]))
                by3.append(float(l[9]))
                bz3.append(float(l[10]))
                l = input.readline().split()
            X = X + bx1 + bx2 + bx3
            Y = Y + by1 + by2 + by3
            Z = Z + bz1 + bz2 + bz3
        x = input.readline()
    input.close()
    
    return freq, X, Y, Z

#source: from afmm_read_gaupy.py:
def Compute(mdfreq, mdX, mdY, mdZ, qmfreq, qmX, qmY, qmZ):
    """Computes the Merit function"""
    N = int((len(mdfreq)/3) + 2)
    if (N != ((len(qmfreq)/3) + 2)) or (N == 2):
        print("Atoms from MD:", N, ". Atoms from QM:", int((len(qmfreq)/3) + 2))
        print("Different or zero number of atoms. Exiting.")
        sys.exit(-1)
    weight = []
    maxprojidxs = []
    nfreq = len(mdfreq)
    # cost matrix, saves qm and mm dot products, row: MD, col: QM
    costmatrix = np.zeros([nfreq, nfreq]) 
    for mdidx in range(nfreq):
        maxproj = 0.0
        maxprojidx = 0
        for qmidx in range(nfreq):
            # gets the correct values and does the dot product
            mdstart = mdidx * N
            qmstart = qmidx * N
            proj = DotProduct(mdX[mdstart:mdstart+N], mdY[mdstart:mdstart+N], mdZ[mdstart:mdstart+N], qmX[qmstart:qmstart+N], qmY[qmstart:qmstart+N], qmZ[qmstart:qmstart+N])
            #build the costmatrix 
            costmatrix[mdidx][qmidx] = proj    #1-proj if not maximizue=True
    #hungarian method, to ensure the best 1:1 mapping of QM and MM freqs
    row_ind, col_ind = linear_sum_assignment(costmatrix, maximize=True)
    maxprojidxs = col_ind

    for pos, maxprojidx in enumerate(maxprojidxs):
        maxproj = costmatrix[pos][maxprojidx]
        #print(maxproj, max(costmatrix[pos]), (maxproj/max(costmatrix[pos]))*100 ,"%")
        mdidx = pos
        # compute the different weights
        # Attention: "projection" does not work, as maxproj often 0, "frequency" also resuletd in error sometimes
        #if self.weighting == "projection":
        #    if self.TooLow(maxproj):
        #        print("Maximum projection very close to zero. Exiting.")
        #        sys.exit(-1)
        #    weight.append(1/maxproj)
        #elif self.weighting == "frequency":
        #    if self.TooLow(mdfreq[mdidx]):
        #        print("MD frequency very close to zero. Exiting.")
        #        sys.exit(-1)
        #    weight.append(1/mdfreq[mdidx])
        #elif self.weighting == "freqmatch":
        #    md = mdfreq[mdidx]
        #    qm = qmfreq[maxprojidx]*self.qmfactor
        #    weight.append(abs((md - qm)/max(md, qm)))
        #elif self.weighting == "none":
        weight.append(1.0)
    sum = 0.0
    nwsum = 0.0
    mdqmfreq = []
    for i in range(nfreq):
        # compute intermediate result to avoid doing it twice
        t = (mdfreq[i] - qmfactor*qmfreq[maxprojidxs[i]])**2
        # weighted sum of squares
        sum = sum + (weight[i]**2) * t
        # non-weighted sum of squares needed for final output
        nwsum = nwsum + t
        mdqmfreq.append(qmfreq[maxprojidxs[i]])
    return (sum/nfreq)**0.5, (nwsum/nfreq)**0.5, mdfreq, mdqmfreq, maxprojidxs

def DotProduct(mdX, mdY, mdZ, qmX, qmY, qmZ):
        sum = 0.0
        for i in range(len(mdX)):
            sum = sum + mdX[i]*qmX[i] + mdY[i]*qmY[i] + mdZ[i]*qmZ[i]
        return abs(sum)

def RunMD(mdexec, mdinp, mdout):
    """Run MD program, checking for normal termination"""
    cmd = mdexec + " < " + mdinp + " 2> mderror.log 1> " + mdout
    status = os.system(cmd)
    if os.WIFEXITED(status):
        if os.WEXITSTATUS(status) != 0:
            print("MD program returned an error! Not exiting.")
            #print("MD program returned an error! Exiting.")
            #sys.exit(-1)
    else:
        print("MD program signalled! Not exiting.")
        #sys.exit(-1)
