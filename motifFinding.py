"""
These functions deal with cluster finding
"""
from clusterScoring import *
"""
This function searches for conserved motiffs of at least length l
"""
def getConservedMotifs(X,Y,minLength):
    locationsLength = {}
    motifsLength = {}
    for i in range(minLength,50):
        locationsLength[i] = [None]*len(X)
        motifsLength[i] = {}
    for i in range(-len(Y),len(X)):
        start = -i
        if start < 0: start = 0
        end = len(X)-i
        if end > len(Y): end = len(Y)
        length = 0
        mismatch1 = start
        mismatch2 = start
        for j in range(start,end):
            if X[i+j] == Y[j]:
                length += 1
            else:
                if length >= minLength: #add to clusters
                    motif = X[i+mismatch2+1:i+mismatch2+length+1]
                    locationsLength[length][i+mismatch2+1] = motif
                    if motifsLength[length].has_key(motif):
                        if i+mismatch2+1 not in motifsLength[length][motif][0]:
                            motifsLength[length][motif][0] += [i+mismatch2+1]
                        if mismatch2+1 not in motifsLength[length][motif][1]:
                            motifsLength[length][motif][1] += [mismatch2+1]
                    else:
                        motifsLength[length][motif] = [[i+mismatch2+1],[mismatch2+1]]
                mismatch2 = mismatch1 #update mismatch
                length = j - mismatch1 #update length
                mismatch1 = j #update mismatch
    return [locationsLength,motifsLength]
"""
This function searches for conserved motiffs of at least length l
"""
def getKnownMotifs(X,Y,pwmsLength,cutoff):
    locationsLength = {}
    motifsLength = {}
    for length in pwmsLength.keys():
        locationsLength[length] = [None]*len(X)
        motifsLength[length] = {}
        for i in range(0,len(X)-length):
            motif = X[i:i+length]
            tf = MaxProbability(motif,pwmsLength)
            if tf[1] > cutoff:
                locationsLength[length][i] = tf[0]
                if motifsLength[length].has_key(tf[0]):
                    motifsLength[length][tf[0]][0] += [i]
                else:
                    motifsLength[length][tf[0]] = [[i],[]]
        for j in range(0,len(Y)-length):
            motif = Y[j:j+length]
            tf = MaxProbability(motif,pwmsLength)
            if tf[1] > cutoff:
                if motifsLength[length].has_key(tf[0]):
                    motifsLength[length][tf[0]][1] += [j]
                else:
                    motifsLength[length][tf[0]] = [[],[j]]
    return [locationsLength,motifsLength]
