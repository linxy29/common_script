# Xinyi Lin, 202206
# Goal: Change the output of get_pySCENIC.r to loom file
# Useage: python /home/linxy29/code/R/common_script/csv2loom.py "day7"

import argparse
import os, sys
import loompy as lp
import numpy as np
import scanpy as sc

parser = argparse.ArgumentParser(description='Process some csv.')
parser.add_argument('sample', help='the sample name')

args = parser.parse_args()

os.getcwd()
os.listdir(os.getcwd()) 


x=sc.read_csv(args.sample + ".csv")
print("Read the csv file.")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create(args.sample + ".loom",x.X.transpose(),row_attrs,col_attrs)
print("Done.")