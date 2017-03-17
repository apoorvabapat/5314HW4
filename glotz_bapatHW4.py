#!/usr/bin/env python
from sys import *
import numpy as np
import sys, getopt, math, os.path, argparse

#TAKES IN FASTA AND COMPARES EVERY SEQUENCE AGAINST ALL OTHERS
#TODO: ADD SCORES INTO MATRIX (HIGHEST SCORES INDICATE GREATEST MATCH) to use with UPGMA
#cant figure out how to take in scoring matrix, code is formatted to expect scoring1.txt

def usage():
 "python hw4.py -f <filename> -g <gap> -s <scoring>

gap_penalty = -1

def matcher(seq1, seq2):
 
    if seq1 == seq2:
	return 3
    if seq1 == '-' or seq2 == '-':
	return gap_penalty
    elif seq1 != seq2:
	return -1


def finalize(align1, align2):
    align1 = align1[::-1]    #reverse sequence 1
    align2 = align2[::-1]    #reverse sequence 2
    
    i,j = 0,0
    
    #score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:                
            symbol = symbol + align1[i]
            score += matcher(align1[i], align2[i])
    
        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += matcher(align1[i], align2[i])
            symbol += ' '
            found = 0
    
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':          
            symbol += ' '
            score += gap_penalty
    
    
    
    #Highest scores indicate best match
    print 'Score =', score, '\n'
    print align1, '\n'
    print align2, '\n'
	
   

def needlemanW(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = np.zeros((m+1, n+1))      # the DP table
   
    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + matcher(seq1[i-1], seq2[j-1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment 
    align1, align2 = '', ''
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: # end toching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + matcher(seq1[i-1], seq2[j-1]):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1
    
    finalize(align1, align2)

# parse the FASTA file and calls the global alignment algorithm for every combination of sequences
def pF(infile):
 sequenceList = []
 with open(infile, 'r') as f:
  for line in f:
   if line.startswith('>'):
    pass
   else:
    temp = ''
    if not line.startswith('>'):
     temp = temp + line.strip('\n')
    sequenceList.append(temp)
 sequenceList = [i.strip(' ') for i in sequenceList]
 #print sequenceList
 i = 0
 j = 1
 for i in iter(sequenceList):
  for j in iter(sequenceList):
   needlemanW(i, j)


 
#main function
def main():
 parse = argparse.ArgumentParser()
 parse.add_argument('-F', '-f', '--filename')
 parse.add_argument('-G', '-g', '--gap')
 parse.add_argument('-S', '-s', '--scoring')
 arg = parse.parse_args()
 infile = arg.filename
 gap_penalty = arg.gap
 matrix = arg.scoring

 programName = parse.prog

	# print statements
 print("\nFile: " + infile + "\n")
 pF(infile)
 

# run the main function
if __name__ == '__main__':
 main()

'''
for i in range(sequenceList[0:]):
  for j in range(sequenceList[1:]):
    needlemanW(int(i), int(j))
'''