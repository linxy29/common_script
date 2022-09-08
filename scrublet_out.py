# Xinyi Lin, 202209
# Goal: Get the cell barcodes of doublets
# Useage: python /home/linxy29/code/R/common_script/csv2loom.py "day7"

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Get the cell barcodes of doublets.')
parser.add_argument('inputPath', help='the cellranger output folder')
parser.add_argument('outputPath', help='the output folder')
parser.add_argument('sample', help='the sample name')

args = parser.parse_args()

counts_matrix = scipy.io.mmread(inputPath + sample + '_unzip/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(inputPath + sample + '_unzip/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
d = {'score': doublet_scores, 'prediction': predicted_doublets}
df = pd.DataFrame(data=d)
df.to_csv(output_path + sample + '_doublet.csv', index=False)
print("Done.")