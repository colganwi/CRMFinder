"""
finds program finds cis regulatory modules by searching for clusters of conserved motiffs.
use GenomeCRMFinder sequence1.fasta sequence2.fasta occupancy.bed
"""
from sys import *
from fileLoading import *
from motifFinding import *
from clusterScoring import *
from clusterFinding import *

"""
This is the main method
"""
def CRMFinder():
    print('')
    if len(argv) != 4:
        print("Error, incorrect number of arguments")
        return
    try:
        file = open(argv[1], 'r')
        sequence1 = loadFasta(file)
        file.close()
    except IOError:
        print("Invalid sequence 1")
        return
    try:
        file = open(argv[2], 'r')
        sequence2 = loadFasta(file)
        file.close()
    except IOError:
        print("Invalid sequence 2")
        return
    try:
        file = open(argv[3], 'r')
        occupancy = loadOccupancy(file)
        file.close()
    except IOError:
        print("Invalid occupancy 1")
        return
    try:
        file = open('parameters.txt', 'r')
        parameters = loadParameters(file) #load parameters file
        file.close()
    except IOError:
        print("Invalid parameters file")
        return
    print('Running CRMFinder with parameters:')
    print(parameters)
    print(' ')
    [locationsLength,motifsLength] = getConservedMotifs(sequence1,sequence2,parameters['length'])
    clusters = findClusters(sequence1,sequence2,parameters,locationsLength,motifsLength,occupancy)
    for cluster in clusters:
        if cluster['score'] > parameters['score']:
            printCluster(sequence1,sequence2,cluster,parameters,0,0)

CRMFinder()
