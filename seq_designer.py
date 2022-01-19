from ast import For, If
import json
import sys
import numpy as np


def ParseJson():
    """
    Parse cadnano json file
    """

    # inputFile = "testsmall.json"
    inputFile = sys.argv[1]
    # outputFile = sys.argv[2]

    with open(inputFile, 'r') as json_data:
        cadnano_data = json.load(json_data)

    filename = cadnano_data['name']
    # print("Importing " + filename + "...")

    strand_data = cadnano_data['vstrands']
    num_strands = len(strand_data)
    length_strands = len(strand_data[0]['scaf'])

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

    return num_strands, length_strands, scaffold, staples, row, col, num, strand_data


def FindStart(scaffold, length_strands):
    """
    Looks for start base and returns if found
    """

    currentBase = [0, 0]  # Start searching at strand 0, base 0
    currentBlock = scaffold[currentBase[0]][currentBase[1]]

    # Go through first empty bases until nonempty base is found
    while currentBlock == [-1, -1, -1, -1]:

        nextBase = [currentBase[0], currentBase[1]+1]
        nextBlock = scaffold[nextBase[0]][nextBase[1]]

        # If entire strand is empty, start searching next strand at base 0
        if(nextBase[1] == length_strands-1):
            nextBase = [currentBase[0]+1, 0]
            nextBlock = scaffold[nextBase[0]][nextBase[1]]

        currentBase = nextBase
        currentBlock = nextBlock

    firstNonEmpty = currentBase
    # From first nonempty base, traverse until end is reached
    while currentBase != [-1, -1]:
        nextBase, nextBlock = TraverseScaffold(scaffold, currentBase)

        # Break if end base is found
        if nextBase == [-1, -1]:
            break
        
        # If a loop is reached, there is no breakpoint
        if nextBase == firstNonEmpty:
            print("Scaffold does not have breakpoint")
            return [-1,-1]

        currentBlock = nextBlock
        currentBase = nextBase

    endBase = currentBase

    # For even strands start base is on the right of end base
    if endBase[0] % 2 == 0:
        startBase = [endBase[0], endBase[1]+1]
    # For odd strands start base is on the left of end base
    else:
        startBase = [endBase[0], endBase[1]-1]

    return startBase


def TraverseScaffold(scaffold, startBase):
    """
    Traverse scaffold by a single base
    """
    currentBase = startBase
    currentBlock = scaffold[currentBase[0]][currentBase[1]]
    nextBase = [currentBlock[2], currentBlock[3]]
    nextBlock = scaffold[nextBase[0]][nextBase[1]]
    return nextBase, nextBlock


def TraversePrintScaffold(scaffold, startBase):
    """
    Traverse and print scaffold
    """
    if startBase == [-1, -1]:
        print("Invalid start base")
        return

    currentBase = startBase
    currentBlock = scaffold[currentBase[0]][currentBase[1]]
    nextBase = [currentBlock[2], currentBlock[3]]

    # Traverse scaffold until nextBase is [-1,-1]
    print("Printing scaffold...")
    while nextBase != [-1, -1]:
        print(currentBlock)
        nextBase, nextBlock = TraverseScaffold(scaffold, currentBase)
        currentBlock = nextBlock
        currentBase = nextBase


num_strands, length_strands, scaffold, staples, row, col, num, strand_data = ParseJson()
startBase = FindStart(scaffold, length_strands)
TraversePrintScaffold(scaffold, startBase)
