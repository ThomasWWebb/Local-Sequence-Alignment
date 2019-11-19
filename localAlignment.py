import numpy
import time
import sys

def dynprog(uniqueChars, subMatrix, seq1, seq2):
    start = time.time()
    seq1 = " " + seq1
    seq2 = " " + seq2
    switched = False
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seq2 = temp
        switched = True
        
    scoringMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)
    backtrackingMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)

    #directionDict = {'L': 1, 'U': 2, 'D': 3}
    for stepper in range(0, len(seq1)):
        backtrackingMatrix[stepper][0] = 2
    for stepper in range(0, len(seq2)):
        backtrackingMatrix[0][stepper] = 1

    charPosDict = {}
    for stepper in range (0, len(uniqueChars)):
        charPosDict[uniqueChars[stepper]] = stepper
    charPosDict["_"] = len(uniqueChars)

    maxScoreCoords = [0,0]
    maxScore = 0

    for rowStepper in range(1, len(seq1)):
        for columnStepper in range(1, len(seq2)):
            seq1SubPos = charPosDict[seq1[rowStepper]]
            seq2SubPos = charPosDict[seq2[columnStepper]]
            
            score = subMatrix[seq1SubPos][seq2SubPos] + scoringMatrix[rowStepper - 1][columnStepper - 1]
            scoringMatrix[rowStepper][columnStepper] = score
            backtrackingMatrix[rowStepper][columnStepper] = 3
            
            if score < 0:
                scoringMatrix[rowStepper][columnStepper] = 0
                backtrackingMatrix[rowStepper][columnStepper] = 0
                score = 0 
            upScore = subMatrix[seq1SubPos][charPosDict["_"]] + scoringMatrix[rowStepper - 1][columnStepper]
            
            if upScore > score:
                scoringMatrix[rowStepper][columnStepper] = upScore
                backtrackingMatrix[rowStepper][columnStepper] = 2
                score = upScore
            leftScore = subMatrix[charPosDict["_"]][seq2SubPos] + scoringMatrix[rowStepper][columnStepper - 1]
            
            if leftScore > score:
                scoringMatrix[rowStepper][columnStepper] = leftScore
                backtrackingMatrix[rowStepper][columnStepper] = 1

            if score > maxScore:
                maxScore = score
                maxScoreCoords = [columnStepper, rowStepper]
    
    print(maxScore)

    direction = ""
    x = maxScoreCoords[1]
    y = maxScoreCoords[0]
    best_score = scoringMatrix[x][y]
    score = best_score
    best_alignment = ["",""]

    seq1Indices = []
    seq2Indices = []
    
    direction = backtrackingMatrix[x][y]
    
    while not score == 0: 
        if direction == 1:
            best_alignment[0] = seq2[x] + best_alignment[0]
            best_alignment[1] = "-" + best_alignment[1]
            y -= 1
        elif direction == 2:
            best_alignment[0] =  "-" + best_alignment[0]
            best_alignment[1] = seq1[y] + best_alignment[1]
            x -=1
        elif direction == 3:
            best_alignment[0] = seq1[x] + best_alignment[0]
            seq1Indices.insert(0, x - 1) 
            best_alignment[1] = seq2[y] + best_alignment[1]
            seq2Indices.insert(0, y - 1) 
            x -=1
            y -=1
        direction = backtrackingMatrix[x][y]
        score = scoringMatrix[x][y]

        
    stop = time.time()
    time_taken=stop-start
    print('Time taken: '+str(time_taken))
    displayAlignment(best_alignment)
    
    if switched:
        return seq2Indices, seq1Indices
    return seq1Indices, seq2Indices
    
    

def dynproglin(uniqueChars, subMatrix, seq1, seq2):
    start = time.time()
    seq1 = " " + seq1
    seq2 = " " + seq2
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seg2 = temp

    charPosDict = {}
    for stepper in range (0, len(uniqueChars)):
        charPosDict[uniqueChars[stepper]] = stepper
    charPosDict["_"] = len(uniqueChars)

    scoringColumns = numpy.zeros([len(seq1), 2], dtype = int)

    maxScoreCoords = [0,0]
    maxScore = 0

    maxScoreCoords = forwardPass(seq1, seq2,charPosDict, subMatrix )

    reversedSeq1 = reverseSeq(seq1, maxScoreCoords[1])
    reversedSeq2 = reverseSeq(seq2, maxScoreCoords[0])

    maxScoreCoords = forwardPass(reversedSeq1, reversedSeq2,charPosDict, subMatrix)

    localSeq1 = reverseSeq(reversedSeq1, maxScoreCoords[1])
    localSeq2 = reverseSeq(reversedSeq2, maxScoreCoords[0])
    
    print(localSeq1)
    print(localSeq2)

    #best_alignment = globalAligmentLin(cutSeq1, cutSeq2, charPosDict, subMatrix)

def reverseSeq(seq, pos):
    return (' ' + ''.join(reversed(seq[1:pos + 1])))

def forwardPass(seq1, seq2,charPosDict, subMatrix):
    
    scoringColumns = numpy.zeros([2, len(seq1)], dtype = int)
    maxScoreCoords = [0,0]
    maxScore = 0
    
    for seq2Stepper in range(1, len(seq2)):
        
        for seq1Stepper in range(1, len(seq1)):
            seq1SubPos = charPosDict[seq1[seq1Stepper]]
            seq2SubPos = charPosDict[seq2[seq2Stepper]]
            
            score = subMatrix[seq1SubPos][seq2SubPos] + scoringColumns[0][seq1Stepper - 1]
            scoringColumns[1][seq1Stepper] = score
            
            if score < 0:
                scoringColumns[1][seq1Stepper] = 0
                score = 0

            upScore = subMatrix[seq1SubPos][charPosDict["_"]] + scoringColumns[1][seq1Stepper - 1]
            if upScore > score:
                scoringColumns[1][seq1Stepper] = upScore
                score = upScore
            
            leftScore = subMatrix[charPosDict["_"]][seq2SubPos] + scoringColumns[0][seq1Stepper]
            if leftScore > score:
                scoringColumns[1][seq1Stepper] = leftScore

            if score > maxScore:
                maxScore = score
                maxScoreCoords = [seq2Stepper, seq1Stepper]
                
        scoringColumns[0] = scoringColumns[1]

    return maxScoreCoords

#def globalAlignmentLin(seq1, seq2, charPosDict, subMatrix):

def displayAlignment(alignment):
    string1 = alignment[0]
    string2 = alignment[1]
    string3 = ''
    for i in range(min(len(string1),len(string2))):
        if string1[i]==string2[i]:
            string3=string3+"|"
        else:
            string3=string3+" "
    print('Alignment ')
    print('String1: '+string1)
    print('         '+string3)
    print('String2: '+string2+'\n\n')
