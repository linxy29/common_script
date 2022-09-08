# Xinyi Lin, 202209
# Goal: Get the cell barcodes of doublets
# Useage: python ~/code/common_script/scrublet_out.py ./ ./ "Murray_b01"

import argparse
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Get the cell barcodes of doublets.')
parser.add_argument('inputPath', help='the folder containing sample_unzip folders')
parser.add_argument('outputPath', help='the output folder')
parser.add_argument('sample', help='the sample name')

args = parser.parse_args()

counts_matrix = scipy.io.mmread(args.inputPath + args.sample + '_unzip/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(args.inputPath + args.sample + '_unzip/features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
d = {'score': doublet_scores, 'prediction': predicted_doublets}
df = pd.DataFrame(data=d)
df.to_csv(args.outputPath + args.sample + '_doublet.csv', index=False)
print("Done.")