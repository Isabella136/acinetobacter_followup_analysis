import pandas as pd
import numpy as np
import sys

blast_output = sys.argv[1:-1]

matrix_precursor = dict()
genes_list = list()
runs_list = list()

for i in range(len(blast_output)):
    with open(blast_output[i]) as curr_output:
        runs_list.append(' - '.join(
            (blast_output[i].split('/')[-3], blast_output[i].split('/')[-2])))
        try:
            tsv_handler = pd.read_csv(curr_output, sep='\t')
        except:
            continue
        for row in tsv_handler.iterrows():
            if row[1][0] in matrix_precursor:
                matrix_precursor[row[1][0]].append(runs_list[-1])
            else:
                genes_list.append(row[1][0])
                matrix_precursor[row[1][0]] = [runs_list[-1]]

matrix = np.zeros((len(runs_list), len(matrix_precursor)))

for col, (gene, runs) in enumerate(matrix_precursor.items()):
    for run in runs:
        row = runs_list.index(run)
        matrix[row][col] = 1

with open(sys.argv[-1], 'w') as matrix_file:
    matrix_file.write(',' + ','.join(genes_list) + '\n')

    for row in range(len(runs_list)):
        matrix_file.write(runs_list[row])
        for col in matrix[row] :
            matrix_file.write(',{:.0f}'.format(col))
        matrix_file.write('\n')