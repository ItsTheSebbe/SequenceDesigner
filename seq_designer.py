import json
import sys
import numpy as np


def ParseJson():
    """
    Parse cadnano json file
    """

    if len(sys.argv) < 3:
        sys.exit(
            "No file provided. Please use \"Python3 seq_designer.py <filename.json> <scaffoldseq.txt>\" ")
    else:
        inputJson = sys.argv[1]
        inputScaffold = sys.argv[2]

    with open(inputScaffold, 'r') as file:
        scaffold_seq = file.read().replace('\n', '')

    # outputFile = sys.argv[2]

    with open(inputJson, 'r') as json_data:
        cadnano_data = json.load(json_data)

    # filename = cadnano_data['name']

    strand_data = cadnano_data['vstrands']
    num_strands = len(strand_data)
    length_strands = len(strand_data[0]['scaf'])

    # Initialize arrays
    scaffold = np.empty(num_strands, dtype=object)
    staples = np.empty(num_strands, dtype=object)
    row = np.empty(num_strands, dtype=object)
    col = np.empty(num_strands, dtype=object)
    num = np.empty(num_strands, dtype=object)

    for i in range(num_strands):
        scaffold[i] = strand_data[i]['scaf']
        staples[i] = strand_data[i]['stap']
        row[i] = strand_data[i]['row']
        col[i] = strand_data[i]['col']
        num[i] = strand_data[i]['num']

    return num_strands, length_strands, scaffold, staples, row, col, num, strand_data, scaffold_seq


def StapleSearch(staples, num_strands, length_strands):
    """
    Returns all startbases of staple strands
    """
    stapleStartBases = []
    # Check all strands
    for i in range(num_strands):
        # Check all bases in each strand
        for j in range(length_strands):

            startSearchBase = [i, j]
            startBase = TraverseEntireReverse(staples, startSearchBase)

            if startBase != [-1, -1]:
                stapleStartBases.append(startBase)

    # Extracts only unique start bases
    stapleStartBases = [list(x) for x in set(tuple(x)
                                             for x in stapleStartBases)]

    # Find corresponding end bases
    stapleEndBases = [None]*len(stapleStartBases)
    for i in range(len(stapleStartBases)):
        stapleEndBases[i] = TraverseEntire(staples, stapleStartBases[i])

    return stapleStartBases, stapleEndBases


def ScaffoldSearch(scaffold, num_strands, length_strands):
    """
    Looks for a start base and returns if found.
    If not found, exits program
    """

    currentBase = [0, 0]  # Start searching at strand 0, base 0
    currentBlock = scaffold[currentBase[0]][currentBase[1]]

    # Go through first empty bases until nonempty base is found
    while currentBlock == [-1, -1, -1, -1]:

        # Shift one base to the right
        nextBase = [currentBase[0], currentBase[1]+1]
        nextBlock = scaffold[nextBase[0]][nextBase[1]]

        # If no scaffold found in all strands
        if(nextBase[1] == length_strands-1 and nextBase[0] == num_strands - 1):
            sys.exit("Scaffold does not exist")

        # If entire strand is empty, start searching next strand at base 0
        if(nextBase[1] == length_strands-1):
            nextBase = [currentBase[0]+1, 0]
            nextBlock = scaffold[nextBase[0]][nextBase[1]]

        currentBase = nextBase
        currentBlock = nextBlock

    firstNonEmpty = currentBase

    # From first nonempty base, traverse until end is reached
    while currentBase != [-1, -1]:
        nextBase, nextBlock = Traverse(scaffold, currentBase)

        # Break if end base is found
        if nextBase == [-1, -1]:
            break

        # If a loop is reached, there is no breakpoint
        if nextBase == firstNonEmpty:
            sys.exit("Scaffold does not have breakpoint")

        currentBlock = nextBlock
        currentBase = nextBase

    endBase = currentBase

    # For even strands start base is on the right of end base
    if endBase[0] % 2 == 0:
        scaffoldStartBase = [endBase[0], endBase[1]+1]
    # For odd strands start base is on the left of end base
    else:
        scaffoldStartBase = [endBase[0], endBase[1]-1]

    # Find corresponding end base
    scaffoldEndBase = TraverseEntire(scaffold, scaffoldStartBase)

    # Check if start and end align
    ValidateScaffold(scaffoldStartBase, scaffoldEndBase)

    return scaffoldStartBase, scaffoldEndBase


def ValidateScaffold(startBase, endBase):
    """
    Check if start and end base connect correctly
    """

    # check if end and start base are next to each other!
    if startBase[0] != endBase[0]:
        sys.exit("Start and end of scaffold do not connect!\nMake sure the start and end bases connect and that you don't have multiple scaffolds")
    else:
        if endBase[0] % 2 == 0 and endBase[1]+1 != startBase[1]:
            sys.exit("Start and end of scaffold do not connect!\nMake sure the start and end bases connect and that you don't have multiple scaffolds")
        elif endBase[0] % 2 == 1 and endBase[1]-1 != startBase[1]:
            sys.exit("Start and end of scaffold do not connect!\nMake sure the start and end bases connect and that you don't have multiple scaffolds")


def Traverse(strand, startBase):
    """
    Traverse strand by a single base
    """

    currentBase = startBase
    currentBlock = strand[currentBase[0]][currentBase[1]]
    nextBase = [currentBlock[2], currentBlock[3]]
    nextBlock = strand[nextBase[0]][nextBase[1]]
    return nextBase, nextBlock


def ReverseTraverse(strand, startBase):
    """
    Traverse strand by a single base in reverse
    """

    currentBase = startBase
    currentBlock = strand[currentBase[0]][currentBase[1]]
    prevBase = [currentBlock[0], currentBlock[1]]
    prevBlock = strand[prevBase[0]][prevBase[1]]
    return prevBase, prevBlock


def TraverseEntire(strand, startSearchBase):
    """
    Traverse strand, returns end base
    """

    currentBase = startSearchBase
    currentBlock = strand[currentBase[0]][currentBase[1]]

    # If current base is empty, return empty base
    if currentBlock == [-1, -1, -1, -1]:
        return [-1, -1]

    nextBase, nextBlock = Traverse(strand, currentBase)

    # Traverse scaffold until prevBase is [-1,-1]
    while nextBase != [-1, -1]:
        currentBlock = nextBlock
        currentBase = nextBase

        nextBase, nextBlock = Traverse(strand, currentBase)

    endBase = currentBase

    return endBase


def TraverseEntireReverse(strand, startSearchBase):
    """
    Traverse strand in reverse, returns start base
    """

    currentBase = startSearchBase
    currentBlock = strand[currentBase[0]][currentBase[1]]

    # If current base is empty, return empty base
    if currentBlock == [-1, -1, -1, -1]:
        return [-1, -1]

    prevBase, prevBlock = ReverseTraverse(strand, currentBase)

    # Traverse scaffold until prevBase is [-1,-1]
    while prevBase != [-1, -1]:
        currentBlock = prevBlock
        currentBase = prevBase

        prevBase, prevBlock = ReverseTraverse(strand, currentBase)

    startBase = currentBase

    return startBase


def PrintScaffold(sequence, view=0):
    """
    Prints scaffold, 0 = detailed view, 1 = cadnano view
    """
    print("Scaffold: ")

    if view == 0:
        for seq in sequence:
            print(seq)
    elif view == 1:
        print("Start,End,Sequence,Length")
        print(str(sequence[0][0]) + "[" + str(sequence[0][1]) + "]," +
              str(sequence[-1][0]) + "[" + str(sequence[-1][1]) + "],", end='')
        cnt = 0
        for j in range(len(sequence)):
            print(str(sequence[j][2]), end='')
            cnt += 1
        print("," + str(cnt))

    print("")


def PrintStaples(sequence, view=0):
    """
    Prints staples, 0 = detailed view, 1 = cadnano view
    """
    print("Staples:")

    if view == 0:
        for i in range(len(sequence)):
            print("Staple " + str(i) + ":")
            for seq in sequence[i]:
                print(seq)
            print("")

    elif view == 1:
        print("Start,End,Sequence,Length")
        for i in range(len(sequence)):
            currentSequence = sequence[i]
            print(str(currentSequence[0][0]) + "[" + str(currentSequence[0][1]) + "]," +
                  str(currentSequence[-1][0]) + "[" + str(currentSequence[-1][1]) + "],", end='')
            cnt = 0
            for j in range(len(currentSequence)):
                print(str(currentSequence[j][2]), end='')
                cnt += 1
            print("," + str(cnt))

    print("")


def FindScaffoldSequence(scaffold, startBase, rawScaffoldSequence):
    """
    Finds and returns scaffold sequence
    """

    finalSequence = []

    currentBase = startBase
    currentBlock = scaffold[currentBase[0]][currentBase[1]]

    currentBase.append(rawScaffoldSequence[0])
    finalSequence.append(currentBase)

    nextBase, nextBlock = Traverse(scaffold, currentBase)
    cnt = 1

    # Traverse scaffold until nextBase is [-1,-1]
    while nextBase != [-1, -1]:

        currentBase = nextBase
        currentBlock = nextBlock

        currentBase.append(rawScaffoldSequence[cnt])
        finalSequence.append(currentBase)

        nextBase, nextBlock = Traverse(scaffold, currentBase)
        cnt += 1

    return finalSequence


def Complement(base):
    """
    Returns the complementary base
    """

    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    elif base == 'C':
        return 'G'
    else:
        sys.exit("Not a valid base")


def FindScaffoldBase(base, scaffoldSequence):
    """
    Runs through scaffold to find the letter at position base. 
    If found, returns complement of that letter
    """

    baseLetter = -1
    for i in range(len(scaffoldSequence)):
        if base[0] == scaffoldSequence[i][0] and base[1] == scaffoldSequence[i][1]:
            baseLetter = Complement(scaffoldSequence[i][2])

    if baseLetter == -1:
        sys.exit("Staple not connected to scaffold")
    return baseLetter


def FindStapleSequences(staples, stapleStartBases, scaffoldSequence):
    """
    Finds staple sequences, returns them in a list of lists
    """
    
    finalSequence = [None] * len(stapleStartBases)

    for i in range(len(stapleStartBases)):

        currentBase = stapleStartBases[i]
        currentBlock = staples[currentBase[0]][currentBase[1]]

        baseLetter = FindScaffoldBase(currentBase, scaffoldSequence)
        currentBase.append(baseLetter)
        finalSequence[i] = [currentBase]

        nextBase, nextBlock = Traverse(staples, currentBase)

        # Traverse scaffold until nextBase is [-1,-1]
        while nextBase != [-1, -1]:

            currentBase = nextBase
            currentBlock = nextBlock

            baseLetter = FindScaffoldBase(currentBase, scaffoldSequence)
            currentBase.append(baseLetter)
            finalSequence[i].append(currentBase)

            nextBase, nextBlock = Traverse(staples, currentBase)

    return finalSequence

def main():
    """
    Main program loop
    """

    # load data
    num_strands, length_strands, scaffold, staples, row, col, num, strand_data, rawScaffoldSequence = ParseJson()

    # staples
    stapleStartBases, stapleEndBases = StapleSearch(
        staples, num_strands, length_strands)

    # scaffold
    scaffoldStartBase, scaffoldEndBase = ScaffoldSearch(
        scaffold, num_strands, length_strands)

    # Returns scaffold sequence
    scaffoldSequence = FindScaffoldSequence(
        scaffold, scaffoldStartBase, rawScaffoldSequence)

    # Returns staple sequences
    stapleSequence = FindStapleSequences(
        staples, stapleStartBases, scaffoldSequence)

    # IO
    PrintScaffold(scaffoldSequence, 1)
    PrintStaples(stapleSequence, 1)

main()
