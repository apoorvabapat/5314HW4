import 
import os, system
import random


def lowest_value(table):
	min=10000
	x,y=-1,-1

	for i in range(len(table)):
		for j in range(len(table[i])):
			if table[i][j]<min:
				min=table[i][j]
				x,y=i,j
	return x,y



def join_table(table,a,b):
	if b<a:
		a,b=b,a

	new_row=[]
	for i in range(0,a):
		new_row.appen d(float(table[a][i]+table[b][i])/2)
	table[a]=new_row

	for i in range(a+1, b):
        print a
        print table[i][a]
        table[i][a] = float((table[i][a]+table[b][i]))/2
        
    #   We get the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]
    print table
