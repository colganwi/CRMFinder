"""
These functions handle file loading
"""
"""
This function loads in a fasta file
"""
def loadFasta(file):
    lines = file.readlines()
    sequence = str(lines[1].strip('\n').strip('\r'))
    return sequence
"""
This function loads in a gff3 file
"""
def loadGff3(file):
    genes = {}
    lines = file.readlines()
    for i in range(4,len(lines)):
         line = lines[i].strip('\n').strip('\r').split('\t')
         if line[2] == 'gene':
             data = line[8]
             genes[int(line[3])] = [int(data.split(';')[1].split('=')[1]),int(line[4])]
    return genes
"""
This function loads in the parameters file.
"""
def loadParameters(file):
    parameters = {}
    lines = file.readlines()
    line = lines[0].strip('\n').strip('\r').split(':')
    parameters[line[0]] = line[1]
    line = lines[1].strip('\n').strip('\r').split(':')
    parameters[line[0]] = float(line[1])
    for i in range(2,5):
        line = lines[i].strip('\n').strip('\r').split(':')
        parameters[line[0]] = int(line[1])
    for i in range(5,10):
        line = lines[i].strip('\n').strip('\r').split(':')
        parameters[line[0]] = float(line[1])
    return parameters
"""
This function loads in the occupancy file.
"""
def loadOccupancy(file):
    occupancy = []
    lines = file.readlines()
    for i in range(1,len(lines)):
        line = lines[i].strip('\n').strip('\r')
        line = line.split('\t')
        occupancy += [float(line[3])]
    return occupancy
"""
This function reads in the pwm file
"""
def loadPwm(file):
    pwms_length = {}
    lines = file.readlines()
    for i in range(0, len(lines), 6):
        tf = lines[i].split()[1]
        stringA = lines[i+1][lines[i+1].index("[")+1:lines[i+1].index("]")].split()
        A = [int(j) for j in stringA]
        length = len(A)
        stringC = lines[i+2][lines[i+2].index("[")+1:lines[i+2].index("]")].split()
        C = [int(j) for j in stringC]
        stringG = lines[i+3][lines[i+3].index("[")+1:lines[i+3].index("]")].split()
        G = [int(j) for j in stringG]
        stringT = lines[i+4][lines[i+4].index("[")+1:lines[i+4].index("]")].split()
        T = [int(j) for j in stringT]
        columnsum = A[0] + C[0] + G[0] + T[0]
        if pwms_length.has_key(length):
            pwms_length[length] += [[tf, [A, C, G, T], columnsum]]
        else:
            pwms_length[length] = [[tf, [A, C, G, T], columnsum]]
    return pwms_length
