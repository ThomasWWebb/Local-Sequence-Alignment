import numpy

def dynprog(uniqueChars, subMatrix, seq1, seq2):
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
        if stepper < len(seq2):
            backtrackingMatrix[0][stepper] = 1
        
    charPosDict = createCharPosDict(uniqueChars)
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

    seq1Indices = []
    seq2Indices = []
    
    direction = backtrackingMatrix[x][y]
    
    while not score == 0:
        if direction == 1:
            y -= 1
        elif direction == 2:
            x -=1
        elif direction == 3:
            seq1Indices.insert(0, x - 1) 
            seq2Indices.insert(0, y - 1) 
            x -=1
            y -=1
        direction = backtrackingMatrix[x][y]
        score = scoringMatrix[x][y]
    
    if switched:
        return best_score, seq2Indices, seq1Indices
    return best_score, seq1Indices, seq2Indices

def createCharPosDict(uniqueChars):
    charPosDict = {}
    for stepper in range (0, len(uniqueChars)):
        charPosDict[uniqueChars[stepper]] = stepper
    charPosDict["_"] = len(uniqueChars)
    return charPosDict
    
##########################################################################
def dynproglin(uniqueChars, subMatrix, seq1, seq2):
    seq1 = " " + seq1
    seq2 = " " + seq2
    switched = False
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seq2 = temp
        switched = True
        

    charPosDict = createCharPosDict(uniqueChars)

    scoringColumns = numpy.zeros([len(seq1), 2], dtype = int)

    maxScoreCoords = [0,0]
    maxScore = 0

    maxScoreCoords, maxScore = forwardPass(seq1, seq2,charPosDict, subMatrix, False )
    originalMaxScoreCoords = maxScoreCoords
    reversedSeq1 = reverseSeqByPos(seq1, maxScoreCoords[1])
    reversedSeq2 = reverseSeqByPos(seq2, maxScoreCoords[0])

    maxScoreCoords, _ = forwardPass(reversedSeq1, reversedSeq2,charPosDict, subMatrix, False)
    localSeq1 = reverseSeqByPos(reversedSeq1, maxScoreCoords[1])
    localSeq2 = reverseSeqByPos(reversedSeq2, maxScoreCoords[0])

    best_alignment = globalAlignmentLin(localSeq1, localSeq2, charPosDict, subMatrix)

    seq1Indices = []
    seq2Indices = []
    stepper = len(best_alignment[0]) - 1
    maxScoreCoords = originalMaxScoreCoords
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

##########################################################################

def heuralign(uniqueChars, subMatrix, seq1, seq2):

    switched = False
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seq2 = temp
        switched = True

    charPosDict = createCharPosDict(uniqueChars)  
    ktup = 2
    bandWidth = 7
    
    seq1KtupDict = getKtupSeqDict(seq1, ktup)
    diagDict = findMatches(seq2, seq1KtupDict, ktup)
    if diagDict == {}:
        ktup = 1
        seq1KtupDict = getKtupSeqDict(seq1, ktup)
        diagDict = findMatches(seq2, seq1KtupDict, ktup)
    diagScoreArray = []
    for key in diagDict:
        diagScore = 0
        score, seq1End, seq2End = -1,-1,-1
        for value in diagDict[key]:
            if seq1End < value[1] and seq2End < value[2]:
                score, seq1End, seq2End = extendDiag(value, seq1, seq2, subMatrix, charPosDict, ktup)
                diagScore += score
        diagScoreArray.append([key, diagScore])
    diagScoreArray.sort(reverse = True, key = lambda x: x[1])

    seq1 = "_" + seq1
    seq2 = "_" + seq2

    scoringMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)
    backtrackingMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)

    for stepper in range(0, len(seq1)):
        backtrackingMatrix[stepper][0] = 2
        if stepper < len(seq2):
            backtrackingMatrix[0][stepper] = 1

    topDiags = []
    bestDiagScore = 0
    bestIndexes = []

    if len(diagScoreArray) < 10:
        topDiags = diagScoreArray[0:len(diagScoreArray)]
    else:
        topDiags = diagScoreArray[0:10]

    for arrayStepper in range(0, len(topDiags)):

        topDiag = diagScoreArray[arrayStepper]
        startCoords = getDiagStartCoords(topDiag[0])
        maxScoreCoords = [0,0]
        maxScore = 0

        seq1Len = len(seq1)
        seq2Len = len(seq2)

        
        rowStepper = startCoords[0] + 1 
        columnStepper = startCoords[1] + 1

        if rowStepper > 1:
            rowStepper = rowStepper - bandWidth
            columnStepper = columnStepper - bandWidth
            if rowStepper < 1:
                while rowStepper < 1:
                    rowStepper += 1
                    columnStepper += 1

        bandedMatrix = numpy.zeros([2, bandWidth * 2 + 2], dtype = int)
        rowArray = numpy.zeros([1,  bandWidth * 2 + 2], dtype = int)
        bandedBacktrackMatrix = numpy.ones([1, bandWidth * 2 + 2], dtype = int)
        bandedBacktrackMatrix = numpy.append(bandedBacktrackMatrix, rowArray, 0)
        bandedBacktrackMatrix[1][0] = 2
        maxSubArrayStepper = len(bandedMatrix[0])

        subRowStepper = 1
        while rowStepper < seq1Len and columnStepper < seq2Len + bandWidth:
            subArrayStepper = 1
            for subColumnStepper in range(columnStepper - bandWidth, columnStepper + bandWidth + 1):
                if subColumnStepper > 0 and subColumnStepper < seq2Len:
                    seq1SubPos = charPosDict[seq1[rowStepper]]
                    seq2SubPos = charPosDict[seq2[subColumnStepper]]
                    
                    score = subMatrix[seq1SubPos][seq2SubPos] + bandedMatrix[subRowStepper - 1][subArrayStepper]
                    bandedMatrix[subRowStepper][subArrayStepper] = score
                    bandedBacktrackMatrix[subRowStepper][subArrayStepper] = 3
                    
                    if score < 0:
                        bandedMatrix[subRowStepper][subArrayStepper] = 0
                        bandedBacktrackMatrix[subRowStepper][subArrayStepper] = 0
                        score = 0

                    if subColumnStepper > columnStepper - bandWidth and subArrayStepper < maxSubArrayStepper -1:
                        upScore = subMatrix[seq1SubPos][charPosDict["_"]] + bandedMatrix[subRowStepper - 1][subArrayStepper + 1]
                        if upScore > score:
                            bandedMatrix[subRowStepper][subArrayStepper] = upScore
                            bandedBacktrackMatrix[subRowStepper][subArrayStepper] = 2
                            score = upScore
                        
                    if subColumnStepper <= columnStepper + bandWidth:
                        leftScore = subMatrix[charPosDict["_"]][seq2SubPos] + bandedMatrix[subRowStepper][subArrayStepper - 1]
                        if leftScore > score:
                            bandedMatrix[subRowStepper][subArrayStepper] = leftScore
                            bandedBacktrackMatrix[subRowStepper][subArrayStepper] = 1

                    if score > maxScore:
                        maxScore = score
                        maxScoreCoords = [subArrayStepper, subRowStepper]
                        maxScoreSeqPos = [rowStepper, subColumnStepper]
                subArrayStepper += 1
            bandedMatrix = numpy.append(bandedMatrix, rowArray, 0)
            bandedBacktrackMatrix = numpy.append(bandedBacktrackMatrix, rowArray, 0)
            
            subRowStepper +=1
            rowStepper += 1
            columnStepper += 1

        if maxScore > bestDiagScore:
            bestDiagScore = maxScore

            direction = ""
            x = maxScoreCoords[1]
            y = maxScoreCoords[0]

            if len(seq1) > 0:
                if seq1[0] == "_":
                    backtrackSeq1 = seq1[1:]
                else:
                    backtrackSeq1 = seq1
            if len(seq2) > 0:
                if seq2[0] == "_":
                    backtrackSeq2 = seq2[1:]
                else:
                    backtrackSeq2 = seq2
            
            seq1X = maxScoreSeqPos[0] - 1
            seq2Y = maxScoreSeqPos[1]
            if seq2Y >= seq2Len:
                seq2Y = seq2Len - 1
            seq2Y -= 1
            best_score = bandedMatrix[x][y]
            score = best_score
            best_alignment = ["",""]

            seq1Indices = []
            seq2Indices = []

            direction = bandedBacktrackMatrix[x][y]
            
            while not score == 0:
                if direction == 1:
                    best_alignment[0] = backtrackSeq2[seq2Y] + best_alignment[0]
                    best_alignment[1] = "-" + best_alignment[1]
                    y -= 1
                    seq2Y -= 1
                elif direction == 2:
                    best_alignment[0] =  "-" + best_alignment[0]
                    best_alignment[1] = backtrackSeq1[seq1X] + best_alignment[1]
                    x -=1
                    y += 1
                    seq1X -= 1
                elif direction == 3:
                    best_alignment[1] = backtrackSeq1[seq1X] + best_alignment[1]
                    seq1Indices.insert(0, seq1X)
                    best_alignment[0] = backtrackSeq2[seq2Y] + best_alignment[0]
                    seq2Indices.insert(0, seq2Y) 
                    x -=1
                    seq1X -= 1
                    seq2Y -= 1
                direction = bandedBacktrackMatrix[x][y]
                score = bandedMatrix[x][y]
            
    if switched:
        return (bestDiagScore, seq2Indices, seq1Indices)
    else:
        return (bestDiagScore, seq1Indices, seq2Indices)
    
def getKtupSeqDict(seq, ktup):
    seqDict = {}
    index = []
    for stepper in range(0, len(seq) + 1 - ktup):
        key = seq[stepper:stepper + ktup]
        if key in seqDict:
            indexArray = seqDict[key]
            indexArray.append(stepper)
            seqDict[key] = indexArray
        else:
            seqDict[key] = [stepper]
    return seqDict

def findMatches(seq, seq1KtupDict, ktup):
    diagDict = {}
    for stepper in range(0, len(seq) + 1 - ktup):
        key = seq[stepper:stepper + ktup]
        if key in seq1KtupDict:
            for index in seq1KtupDict[key]:
                diag = index - stepper
                if diag in diagDict:
                    diagArray = diagDict[diag]
                    diagArray.append([key, index, stepper])
                    diagDict[diag] = diagArray
                else:
                    diagDict[diag] = [[key, index, stepper]]
    return diagDict   

def extendDiag(diagPair, seq1, seq2, subMatrix, charPosDict, ktup):
    score = scoreSeq(diagPair[0], diagPair[0], subMatrix, charPosDict)
    bestScore = score

    startIndex = [diagPair[1], diagPair[2]]
    endIndex = [diagPair[1] + ktup - 1, diagPair[2] + ktup - 1]
    threshold = -1
    savedCoords = [startIndex[0], startIndex[1], endIndex[0], endIndex[1]]

    seq1Len = len(seq1)
    seq2Len = len(seq2)

    while score >= threshold and ((startIndex[0] and startIndex[1]) or (endIndex[0] < seq1Len - 1 and endIndex[1] < seq2Len - 1)):

        if score >= threshold and startIndex[0] and startIndex[1]:
            startIndex[0] -= 1
            startIndex[1] -= 1
            seq1Pos = charPosDict[seq1[startIndex[0]]]
            seq2Pos = charPosDict[seq2[startIndex[1]]]
            score += subMatrix[seq1Pos][seq2Pos]
            if score > bestScore:
                savedCoords = [startIndex[0], startIndex[1], savedCoords[2], savedCoords[3]]
                bestScore = score
            seq1Sub = seq1[startIndex[0]:endIndex[0] +1]
            seq2Sub = seq2[startIndex[1]:endIndex[1] +1]
        if score >= threshold and endIndex[0] < seq1Len - 1 and endIndex[1] < seq2Len - 1:
            endIndex[0] += 1
            endIndex[1] += 1
            seq1Pos = charPosDict[seq1[endIndex[0]]]
            seq2Pos = charPosDict[seq2[endIndex[1]]]
            score += subMatrix[seq1Pos][seq2Pos]
            if score > bestScore:
                savedCoords = [savedCoords[0], savedCoords[1], endIndex[0], endIndex[1]]
                bestScore = score
            seq1Sub = seq1[startIndex[0]:endIndex[0] +1]
            seq2Sub = seq2[startIndex[1]:endIndex[1] +1]
    
    seq1Sub = seq1[savedCoords[0]:savedCoords[2] +1]
    seq2Sub = seq2[savedCoords[1]:savedCoords[3] +1]
    return bestScore, savedCoords[2], savedCoords[3]
        
def scoreSeq(seq1, seq2, subMatrix, charPosDict):
    stepper = 0
    score = 0
    while stepper < len(seq1) and stepper < len(seq2):
        seq1Pos = charPosDict[seq1[stepper]]
        seq2Pos = charPosDict[seq2[stepper]]
        score = score + subMatrix[seq1Pos][seq2Pos]
        stepper +=1
    return score

def getDiagStartCoords(diag):
    if diag == 0:
        return (0,0)
    elif diag < 0:
        return(0, -diag)
    return (diag, 0) 
