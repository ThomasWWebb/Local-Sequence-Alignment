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

    maxScoreCoords = [0,0]
    maxScore = 0

    for rowStepper in range(1, len(seq1)):
        for columnStepper in range(1, len(seq2)):
            score = subMatrix[charPosDict[seq1[rowStepper]]][charPosDict[seq2[columnStepper]]] + scoringMatrix[rowStepper - 1][columnStepper - 1]
            scoringMatrix[rowStepper][columnStepper] = score
            backtrackingMatrix[rowStepper][columnStepper] = 3
            if score < 0:
                scoringMatrix[rowStepper][columnStepper] = 0
                backtrackingMatrix[rowStepper][columnStepper] = 0
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

            if score > maxScore:
                maxScore = score
                maxScoreCoords = [columnStepper, rowStepper]

    print(subMatrix)
    print(scoringMatrix)
    print(backtrackingMatrix)
    print(maxScoreCoords)
    
    direction = ""
    x = maxScoreCoords[1]
    y = maxScoreCoords[0]
    best_score = scoringMatrix[x][y]
    score = best_score
    best_alignment = ["",""]
    
    direction = backtrackingMatrix[x][y]
    
    while not score == 0: 
        print(direction)
        print("score = {}".format(scoringMatrix[x][y]))
        if direction == 1:
            best_alignment[0] = seq2[x] + best_alignment[0]
            best_alignment[1] = "-" + best_alignment[1]
            y -= 1
        elif direction == 2:
            best_alignment[0] =  "-" + best_alignment[0]
            best_alignment[1] = seq1[y] + best_alignment[1]
            x -=1
            print("new score up = {}".format(scoringMatrix[x][y]))
        elif direction == 3:
            best_alignment[0] = seq1[x] + best_alignment[0]
            best_alignment[1] = seq2[y] + best_alignment[1]
            x -=1
            y -=1
        direction = backtrackingMatrix[x][y]
        score = scoringMatrix[x][y]
        print(best_alignment)
    displayAlignment(best_alignment)

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
