# Xinyi Lin, 202310
# Goal: convert the excel file from GSA (https://ngdc.cncb.ac.cn/gsa-human/) to txt file for further download
# Usage: python GSA_xlsx2txt.py -i <input file path> -o <output file path> -s <sheet name>

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Convert the excel file from GSA (https://ngdc.cncb.ac.cn/gsa-human/) to txt file for further download')
parser.add_argument('-i', '--input', help='input file path', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)
parser.add_argument('-s', '--sheet', help='sheet name', default='Run', required=False)

args = parser.parse_args()

# Read the Excel file
file_path = args.input  
sheet_name = args.sheet
df = pd.read_excel(file_path, sheet_name=sheet_name, engine='openpyxl')

# Creating the new concatenated rows
new_rows = []
for _, row in df.iterrows():
    new_rows.append(f"{row['Accession']}\t{row['File name 1']}")
    new_rows.append(f"{row['Accession']}\t{row['File name 2']}")

# Write the results to a .txt file
output_file_path = args.output
with open(output_file_path, 'w') as file:
    file.write("\n".join(new_rows))
