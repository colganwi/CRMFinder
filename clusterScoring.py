"""
These function deal with scoring cluster
"""
import numpy as np
"""
This function does insertionsort
"""
def insertionsort(sequence, length):
    score = 0.0
    for index in range(1,len(sequence)):
        currentvalue = sequence[index]
        position = index
        while position>0 and sequence[position-1]>currentvalue:
             sequence[position]=sequence[position-1]
             position = position-1
        if sequence[position] != currentvalue:
            sequence[position] = currentvalue
            score +=1
    return (score/length)
"""
This function calculates the order score
"""
def scorecalculator(x, y):
    a = []
    b = []
    c = [] # corrected
    # creating a better way of picturing our sequences
    for i in range(len(x)):
        if (i == 0 or x[i] != x[i-1]) and x[i] != "-":
            a += x[i]
    for j in range(len(y)):
        if (j == 0 or y[j] != y[j-1]) and y[j] != "-":
            b += y[j]
    if len(a) == 0:
        return 1
    #setting the values of the characters in sequence one to be an arbitrary increasing sequence eg: 1,2,3,4...
    temp = 1
    #creating a dictionary to store the new mapping of conserved regions to arbitray numbers
    dictionary = {}
    for ch in a:
        dictionary[ch] = temp
        temp += 1
    for ch in b:
        value = dictionary.get(ch)
        c.append(value)
    length = len(c)
    score = insertionsort(c, length)
    return score
"""
This function return a subsequence
"""
def getSubsequence(sequence,location):
    return sequence[0][location[0]-sequence[1][0]:location[1]-sequence[1][0]]
"""
This function calculates the occupancyScore
"""
def getOccupancyScore(i,length,occupancy):
    return float(sum(occupancy[i:i+length]))/length
"""
finds the maximum porbability score for a given input sequence
"""
def MaxProbability(motif, pwms_length):
    if not pwms_length.has_key(len(motif)):
        return -1000
    pwms = pwms_length[len(motif)]
    best_likelihood = -1000
    best_tf = None
    for binding_site in pwms:
        prob = 1.0
        pwm = binding_site[1]
        columnsum = binding_site[2]
        for i in range(0,len(motif)):
            if motif[i] == 'A':
                prob *= float(pwm[0][i])/columnsum
            elif motif[i] == 'C':
                prob *= float(pwm[1][i])/columnsum
            elif motif[i] == 'G':
                prob *= float(pwm[2][i])/columnsum
            else:
                prob *= float(pwm[3][i])/columnsum
        if prob == 0.0:
            likelihood = -1000
        else:
            likelihood = np.log(prob/(.25*len(motif)))
        if likelihood > best_likelihood:
            best_likelihood = likelihood
            best_tf = binding_site[0]
    return [best_tf,best_likelihood]
