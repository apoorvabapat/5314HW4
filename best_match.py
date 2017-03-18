#!/usr/bin/env python
from sys import *
import numpy as np
import sys, getopt, math, os.path, argparse
import os
import random

#TAKES IN FASTA AND COMPARES EVERY SEQUENCE AGAINST ALL OTHERS
#TODO: ADD SCORES INTO MATRIX (HIGHEST SCORES INDICATE GREATEST MATCH) to use with UPGMA
#cant figure out how to take in scoring matrix, code is formatted to expect scoring1.txt

def usage():
 "python hw4.py -f <filename> -g <gap> -s <scoring>"
 
gap_penalty = -1
def matcher(seq1, seq2):
 
    if seq1 == seq2:
	return 3
    if seq1 == '-' or seq2 == '-':
	return gap_penalty
    elif seq1 != seq2:
	return -1


def finalize(align1, align2,x,y,distance_matrix):
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

    # print 'Score=', score, '\n'
    # print align1, '\n'
    # print align2, '\n'
    distance_matrix[x][y]=score


	
   

def needlemanW(seq1, seq2,x,y,distance_matrix):
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
    
    finalize(align1, align2,x,y,distance_matrix)

# parse the FASTA file and calls the global alignment algorithm for every combination of sequences
def pF(infile):
 sequenceList = []
 sequence_name=[]
 seq=0
 with open(infile, 'r') as f:
  for line in f:
   if line.startswith('>'):
       line=line.split()
       sequence_name.append(line[0][1:])
       seq+=1
       pass
   else:
    temp = ''
    if not line.startswith('>'):
     temp = temp + line.strip('\n')
    #Keeps track of the sequence names so as to represent in the Newick tree
    sequenceList.append(temp)
 sequenceList = [i.strip(' ') for i in sequenceList]
 print sequence_name
 n=len(sequenceList)
 distance_matrix=np.zeros((n,n))

 x = 0
 y = 0
 for i in iter(sequenceList):
  y=0
  for j in iter(sequenceList):
   needlemanW(i, j,x,y,distance_matrix)
   y+=1
  x+=1


#Converts the matrix into lower trainagular matrix for applying UPGMA
 M=[]
 for i in range(len(distance_matrix)):
    list1=[]
    for j in range(0,i):
        list1.append(distance_matrix[i][j])
    M.append(list1)
 return UPGMA(M,sequence_name)

#Returns the minimum value in the table
def lowest_value(table):
    min=10000
    x,y=-1,-1

    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j]<min:
                min=table[i][j]
                x,y=i,j
    #Returns the lowest value row and column
    return x,y


def edit_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


def edit_table(table,a,b):
    if b<a:
        a,b=b,a

    new_row=[]
    for i in range(0,a):
        new_row.append(float(table[a][i]+table[b][i])/2)
    table[a]=new_row

    for i in range(a+1, b):
        table[i][a] = float((table[i][a]+table[b][i]))/2
        
    #   We get the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        # Remove the (repeated second index column entry
        del table[i][b]

    # Remove the (repeated) second index row
    del table[b]    


def UPGMA(table, labels):

    # Until all labels have been joined
    while len(labels) > 1:
        # Locate lowest cell in the table
        x, y = lowest_value(table)

        # Join the table on the cell co-ordinates
        edit_table(table, x, y)

        # Update the labels accordingly
        edit_labels(labels, x, y)

    # Return the final label
    return labels[0]


 
#main function
def main():
 parse = argparse.ArgumentParser()
 #print parse
 parse.add_argument('-F', '-f', '--filename')
 parse.add_argument('-G', '-g', '--gap')
 parse.add_argument('-S', '-s', '--scoring')
 parse.add_argument('-T', '-t', '--tree')
 arg = parse.parse_args()
 infile = arg.filename
 gap_penalty = arg.gap
 matrix = arg.scoring
 tree=arg.tree
 wr=open(tree,'w')
 programName = parse.prog

	# print statements
 print("\nFile: " + infile + "\n")
 label=pF(infile)
 print label
 wr.write(label)


 

# run the main function
if __name__ == '__main__':
 main()

'''
for i in range(sequenceList[0:]):
  for j in range(sequenceList[1:]):
    needlemanW(int(i), int(j))
'''


