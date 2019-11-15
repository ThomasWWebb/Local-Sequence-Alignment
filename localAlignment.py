import numpy

def dynprog(uniqueChars, subMatrix, seq1, seq2):
    seq1 = " " + seq1
    seq2 = " " + seq2
    
    if len(seq1) < len(seq2):
        temp = seq1
        seq1 = seq2
        seg2 = temp
        

    scoringMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)
    backtrackingMatrix = numpy.zeros([len(seq1), len(seq2)], dtype = int)


    scoreList = []
    #directionDict = {'L': 1, 'U': 2, 'D': 3}
    for stepper in range(0, len(seq1)):
        backtrackingMatrix[stepper][0] = 2
    for stepper in range(0, len(seq2)):
        backtrackingMatrix[0][stepper] = 1

    charPosDict = {}
    for stepper in range (0, len(uniqueChars)):
        charPosDict[uniqueChars[stepper]] = stepper
    charPosDict["_"] = len(uniqueChars)
    print(charPosDict)

    for rowStepper in range(1, len(seq1)):
        for columnStepper in range(1, len(seq2)):
            score = subMatrix[charPosDict[seq1[rowStepper]]][charPosDict[seq2[columnStepper]]] + scoringMatrix[rowStepper - 1][columnStepper - 1]
            scoringMatrix[rowStepper][columnStepper] = score
            if score < 0:
                scoringMatrix[rowStepper][columnStepper] = 0
                backtrackingMatrix[rowStepper][columnStepper] = 3
                score = 0 
            upScore = subMatrix[charPosDict[seq1[rowStepper]]][charPosDict["_"]] + scoringMatrix[rowStepper - 1][columnStepper]
            if upScore > score:
                scoringMatrix[rowStepper][columnStepper] = upScore
                backtrackingMatrix[rowStepper][columnStepper] = 2
                score = upScore
            leftScore = subMatrix[charPosDict["_"]][charPosDict[seq2[columnStepper]]] + scoringMatrix[rowStepper][columnStepper - 1]
            if leftScore > score:
                scoringMatrix[rowStepper][columnStepper] = leftScore
                backtrackingMatrix[rowStepper][columnStepper] = 1
                

    '''
    for rowStepper in range(1, maxLen -1):
        basicMatch = 0       
        for columnStepper in range(1, maxLen -1):                
            if seq1[columnStepper] == seq2[rowStepper]:
                if seq1[columnStepper] == 'A' or seq1[columnStepper] == 'C':
                    basicMatch = 3
                if seq1[columnStepper] == 'G' or seq1[columnStepper] == 'T':
                    basicMatch = 2
            elif seq1[columnStepper] == ' ' or seq2[rowStepper] == ' ':
                basicMatch = -2 
            elif not (seq1[columnStepper] == seq2[rowStepper]):
                basicMatch = -1
                
            scoreList = []
            
            tempScore = scoringMatrix[rowStepper + 1][columnStepper] -2
            scoreList.append(tempScore)
              
            tempScore = scoringMatrix[rowStepper][columnStepper + 1] -2
            scoreList.append(tempScore)

            tempScore = scoringMatrix[rowStepper][columnStepper] + basicMatch
            scoreList.append(tempScore)
            scoringMatrix[rowStepper + 1][columnStepper + 1]  = max(scoreList)
            index = scoreList.index(max(scoreList))     
            backtrackingMatrix[rowStepper + 1][columnStepper + 1] = (index + 1 )
'''
    print(scoringMatrix)
    print(backtrackingMatrix)
