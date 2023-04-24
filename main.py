import sys
from Levenshtein import distance
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import numpy as np


def parseArgs(argv):
    '''parse out Command line options.'''
    try:
        
        parser = ArgumentParser(description="a program to calculate Fibonacci's number", formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-f", "--fasta_file", dest="fastafile", action="store", help="fasta file for which you want to calc GC% [default: %(default)s]")
        parser.add_argument("-m", "--max_dist", dest="maxdist", action="store", help="filter Levenshtein distances greater than this [default: %(default)s]")
    
        # Process arguments
        args = parser.parse_args()

        global fastaFile
        global maxDist
        
        fastaFile = args.fastafile
        maxDist = args.maxdist

        # check the user specified a fasta file, if not warn and and exit
        if fastaFile:
            print("fasta file is <" + fastaFile + ">")
        else:
            print("you must specify a fasta file")
            exit
            
        if maxDist:                    
            print("max Levenshtein distance to keep is <" + str(maxDist) + ">")
        
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        


def levenshtein(seq1, seq2):
    '''
    simple implementation of levenshtein distance calculation
    '''
    xdim = len(seq1) + 1
    ydim = len(seq2) + 1
    leven_matrix = np.zeros((xdim, ydim))
    x=0
    while x< xdim:
        leven_matrix[x, 0]=x
        x+=1
    y=0
    while y< ydim:
        leven_matrix[0, y]=y
        y+=1

    x=0
    while x< xdim:
        y=0
        while y < ydim:
            if seq1[x-1] == seq2[y-1]:
                leven_matrix [x,y] = min(
                    leven_matrix[x-1, y] + 1,
                    leven_matrix[x-1, y-1],
                    leven_matrix[x, y-1] + 1
                )
            else:
                leven_matrix [x,y] = min(
                    leven_matrix[x-1,y] + 1,
                    leven_matrix[x-1,y-1] + 1,
                    leven_matrix[x,y-1] + 1
                )

    return (leven_matrix[xdim - 1, ydim-1])

        

def generateDistanceMatrix():
    
    '''
    generate matrix of Levenshtein distances between all pairs
    '''
    
    ldDistMatrix = np.zeros((len(sequenceLines),len(sequenceLines)))
    
    i=0
    j=0
    
    while i < len(sequenceLines):
        j=i+1
        while j < len(sequenceLines):
            ldDistMatrix[i][j] = distance(sequenceLines[i], sequenceLines[j])
            j+=1
        i+=1
        
    return ldDistMatrix




def readFastaFile(filename):
    '''
    load specified fasta file and store header and sequence as entries in two lists
    :param self:
    :return:
    '''

    print("load sequences from fasta file <" + fastaFile + ">")
    global headerLines
    global sequenceLines

    # load the fasta lines into a list
    try:
        fFA = open(filename, 'r')
        fastaLines = fFA.readlines()
        fFA.close()
    except Exception as e:
        raise(e)

    headerLines = []
    headerLine = ""
    sequenceLines = []
    sequence = ""

    s = 0
    for fastaLine in fastaLines:
        if fastaLine[0] == '>':
            if s > 0:
                headerLines.append(headerLine)
                sequenceLines.append(sequence)
                sequence = ""
            headerLine = fastaLine[1:].strip()
            sequence = ""
            
        else:
            sequence = sequence + fastaLine.strip()
        s += 1

    headerLines.append(headerLine)
    sequenceLines.append(sequence)   
         
    print("loaded <" + str(s-1) + "> sequences")
    
    return len(headerLines)




def testPlot(DistMatrix):
    
    import networkx as nx
    import pandas as pd
    import matplotlib.pyplot as plt
    
    G=nx.from_numpy_matrix(DistMatrix)
    weights = pd.DataFrame(DistMatrix).reset_index().melt('index').values.tolist()
    tuples = [tuple(w) for w in weights]
    
    G.add_weighted_edges_from(tuples)
    pos=nx.spring_layout(G)

    edge_weight = nx.get_edge_attributes(G,'weight')
    d = nx.degree(G)
    node_adjacencies = []
    node_text = []
    for adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))


    nx.draw(G, pos, with_labels=True, font_weight='bold', node_size=[v * 100 for v in dict(d).values()], node_color=node_adjacencies)
    
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_weight)

    plt.show()


def main(argv=None): 

    if argv is None:
        argv = sys.argv
        
    
    # parse_args to get filename
    parseArgs(argv)
    
    # load fasta file
    readFastaFile(fastaFile)
    
    # generate Levenshtein distance matrix
    ldDistMat = generateDistanceMatrix()
    
    # make a pretty plot
    testPlot(ldDistMat)
    
if __name__ == '__main__':

    sys.exit(main())