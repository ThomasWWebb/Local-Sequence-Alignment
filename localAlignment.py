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
            best_alignment[0] = seq2[y] + best_alignment[0]
            best_alignment[1] = "-" + best_alignment[1]
            y -= 1
        elif direction == 2:
            best_alignment[0] =  "-" + best_alignment[0]
            best_alignment[1] = seq1[x] + best_alignment[1]
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
    
    if switched:
        return best_score, seq2Indices, seq1Indices
    return best_score, seq1Indices, seq2Indices
    
    

def dynproglin(uniqueChars, subMatrix, seq1, seq2):
    start = time.time()
    seq1 = " " + seq1
    seq2 = " " + seq2
    switched = False
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seq2 = temp
        switched = True
        

    charPosDict = {}
    for stepper in range (0, len(uniqueChars)):
        charPosDict[uniqueChars[stepper]] = stepper
    charPosDict["_"] = len(uniqueChars)

    scoringColumns = numpy.zeros([len(seq1), 2], dtype = int)

    maxScoreCoords = [0,0]
    maxScore = 0

    maxScoreCoords, maxScore = forwardPass(seq1, seq2,charPosDict, subMatrix, False )
    reversedSeq1 = reverseSeqByPos(seq1, maxScoreCoords[1])
    reversedSeq2 = reverseSeqByPos(seq2, maxScoreCoords[0])

    maxScoreCoords, _ = forwardPass(reversedSeq1, reversedSeq2,charPosDict, subMatrix, False)
    localSeq1 = reverseSeqByPos(reversedSeq1, maxScoreCoords[1])
    localSeq2 = reverseSeqByPos(reversedSeq2, maxScoreCoords[0])

    best_alignment = globalAlignmentLin(localSeq1, localSeq2, charPosDict, subMatrix)

    seq1Indices = []
    seq2Indices = []
    stepper = len(best_alignment[0]) - 1
    maxScoreCoords[0] -= 1
    maxScoreCoords[1] -= 1
    

    while stepper > -1:
        if (best_alignment[0])[stepper] != "-" and (best_alignment[1])[stepper] != '-':
            seq1Indices.insert(0, maxScoreCoords[1])
            maxScoreCoords[1] -= 1
            seq2Indices.insert(0, maxScoreCoords[0])
            maxScoreCoords[0] -= 1
        elif (best_alignment[0])[stepper] == "-":
            maxScoreCoords[0] -= 1
        else:
            maxScoreCoords[1] -= 1
        stepper -= 1

    if switched:
        return (maxScore, seq2Indices, seq1Indices)
    else:
        return (maxScore, seq1Indices, seq2Indices)

def reverseSeqByPos(seq, pos):
    return (' ' + ''.join(reversed(seq[1:pos + 1])))

def reverseSeq(seq):
    return ''.join(reversed(seq))

def forwardPass(seq1, seq2,charPosDict, subMatrix, returnColumn):
    seq1, seq2 = addSpace(seq1, seq2)
    
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

    if returnColumn:
        return scoringColumns[0]
    else:
        return maxScoreCoords, maxScore

def globalAlignmentLin(seq1, seq2, charPosDict, subMatrix):
    if len(seq1) > 0:
        if seq1[0] == " ":
            seq1 = seq1[1:]
    if len(seq2) > 0:
        if seq2[0] == " ":
            seq2 = seq2[1:]

    z = ""
    w = ""
    if len(seq1) == 0:
        for stepper in range(0, len(seq2)):
            z = z + "-"
            w = w + seq2[stepper]
    elif len(seq2) == 0:
        for stepper in range(0, len(seq1)):
            z = z + seq1[stepper]
            w = w + "-"
    elif len(seq1) == 1 or len(seq2) == 1:
        (z, w) = needlemanWunsch(seq1, seq2, charPosDict, subMatrix)
    else:
        seq2Mid = len(seq2) // 2
        columnL = forwardPass(seq1, seq2[0:seq2Mid], charPosDict, subMatrix, True)
        columnL = columnL
        columnR = forwardPass(reverseSeq(seq1),reverseSeq(seq2[seq2Mid:]), charPosDict, subMatrix, True)
        columnR = columnR[::-1]
        res_list = [columnL[i] + columnR[i] for i in range(len(columnL))]     
        seq1Mid = numpy.argmax(res_list)
        z1, w1 = globalAlignmentLin(seq1[:seq1Mid], seq2[0:seq2Mid], charPosDict, subMatrix)
        z2, w2 = globalAlignmentLin(seq1[seq1Mid:], seq2[seq2Mid:], charPosDict, subMatrix)
        z = z1 + z2
        w = w1 + w2
    return z,w

def addSpace(seq1, seq2):
    if seq1[0] != " ":
        seq1 = " " + seq1
    if seq2[0] != " ":
        seq2 = " " + seq2
    return seq1, seq2
        

def needlemanWunsch(seq1, seq2, charPosDict, subMatrix):
    seq1, seq2 = addSpace(seq1, seq2)

    switched = False
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seq2 = temp
        switched = True

    charPos = charPosDict[seq2[1]]
    
    scoringMatrix = numpy.zeros([len(seq2), len(seq1)], dtype = int)
    backtrackingMatrix = numpy.zeros([len(seq2), len(seq1)], dtype = int)
    scoringMatrix[1, 0] = subMatrix[charPos][charPosDict["_"]]
    backtrackingMatrix[1, 0] = 2
    
    for stepper in range(1, len(seq1)):
        scoringMatrix[0, stepper] = subMatrix[charPosDict[seq1[stepper]]][charPosDict["_"]] + scoringMatrix[0][stepper - 1]
        backtrackingMatrix[0, stepper] = 1

    for columnStepper in range(1, len(seq1)):
        seq1SubPos = charPosDict[seq1[columnStepper]]
        
        score = subMatrix[seq1SubPos][charPos] + scoringMatrix[0][columnStepper - 1]
        scoringMatrix[1][columnStepper] = score
        backtrackingMatrix[1][columnStepper] = 3
        
        upScore = subMatrix[charPos][charPosDict["_"]] + scoringMatrix[0][columnStepper]
        
        if upScore > score:
            scoringMatrix[1][columnStepper] = upScore
            backtrackingMatrix[1][columnStepper] = 2
            score = upScore
        leftScore = subMatrix[charPosDict["_"]][seq1SubPos] + scoringMatrix[1][columnStepper - 1]
        
        if leftScore > score:
            scoringMatrix[1][columnStepper] = leftScore
            backtrackingMatrix[1][columnStepper] = 1


    best_alignment = ["",""]
    x = 1
    y = len(seq1) - 1
    direction = backtrackingMatrix[x][y]
    
    while not direction == 0:
        if direction == 1:
            best_alignment[0] = seq1[y] + best_alignment[0]
            best_alignment[1] = "-" + best_alignment[1]
            y -= 1
        elif direction == 2:
            best_alignment[0] =  "-" + best_alignment[0]
            best_alignment[1] = seq2[x] + best_alignment[1]
            x -=1
        elif direction == 3:
            best_alignment[0] = seq1[y] + best_alignment[0]
            best_alignment[1] = seq2[x] + best_alignment[1]
            x -=1
            y -=1
        direction = backtrackingMatrix[x][y]

    if switched:
        return best_alignment[1], best_alignment[0]
    return best_alignment

    
    
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
