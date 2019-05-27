"""
These functions deal with cluster finding
"""
from clusterScoring import *
"""
This function finds and scores all clusters of given length
"""
def findClusters(X,Y,parameters,locationsLength,motifsLength,occupancy):
    size = parameters['region']
    clusters = []
    for i in range(0,len(X)-size,size/4):
        for j in range(0,len(Y)-size,size/4):
            xFound = ['-']*size
            yFound = ['-']*size
            matching = 0
            matchingScore = 0
            notMatching = 0
            tfs = {}
            etsSites = 0
            for length in locationsLength.keys()[::-1]:

                for k in range(0,size-length):
                    motif = locationsLength[length][i+k]
                    if motif != None:
                        for yLocation in motifsLength[length][motif][1]:
                            if yLocation >= j and yLocation < j+size-length and xFound[k] == '-' and xFound[k+length] == '-' and yFound[yLocation-j] == '-' and yFound[yLocation-j+length] == '-':
                                matching += 1
                                matchingScore += length
                                for m in range(0,length):
                                    xFound[k+m] = matching
                                    yFound[yLocation-j+m] = matching
                            else:
                                notMatching += 1
            region = X[i:i+size+1]
            if 'GGAA' in region:
                etsSites += 1
            if 'GGAT' in region:
                etsSites += 1
            if 'ATCC' in region:
                etsSites += 1
            if 'TTCC' in region:
                etsSites += 1
            if 'ATTA' in region:
                etsSites += 1
            if 'TAAT' in region:
                etsSites += 1
            if occupancy != None:
                occupancyScore = getOccupancyScore(i,length,occupancy)
            else:
                occupancyScore = 0;
            xFound = ''.join(str(e) for e in xFound)
            yFound = ''.join(str(e) for e in yFound)
            orderScore = scorecalculator(xFound,yFound)
            matchingScore = float(matchingScore)/size
            #distance = abs((float(i)-float(j))/len(X))
            distance = abs((float(i)-float(939))/len(X))
            score = parameters['a']*matchingScore-parameters['b']*distance-parameters['c']*orderScore-parameters['d']*occupancyScore+parameters['e']*etsSites
            clusters += [{'score':score,'xStart':i,'yStart':j,'xFound':xFound,'yFound':yFound,'matching':matchingScore,'tfs':tfs,'occupancyScore':occupancyScore,'orderScore':orderScore,'distance':distance,'etsSites':etsSites}]
    return clusters

"""
This function prints out a give cluster
"""
def printCluster(X,Y,cluster,parameters,xLocation,yLocation):
    print('Cluster Score: '+str(cluster['score']))
    print(str(cluster['xStart']+xLocation)+' '+X[cluster['xStart']:cluster['xStart']+parameters['region']])
    print(str(cluster['xStart']+xLocation)+' '+cluster['xFound'])
    print(str(cluster['yStart']+yLocation)+' '+Y[cluster['yStart']:cluster['yStart']+parameters['region']])
    print(str(cluster['yStart']+yLocation)+' '+cluster['yFound'])
    print('Matching Score: '+str(round(cluster['matching'],3))+' Distance Score: '+str(round(cluster['distance'],3))+' Order Score: '+str(round(cluster['orderScore'],3))+' Occupancy score: '+str(round(cluster['occupancyScore'],3))+' Ets sites: '+str(cluster['etsSites']))
    print(' ')
