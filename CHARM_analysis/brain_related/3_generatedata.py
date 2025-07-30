"""
Scripts for generate regression data.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
import multiprocessing
multiprocessing.set_start_method('fork')
import concurrent.futures

import numpy as np
import pandas as pd
import sys
import re
import tqdm
import pickle

sys.path.append('/shareb/zliu/analysis/')
sys.path.append('/shareb/zliu/analysis/CHARMtools')
from CHARMtools import Cell3Ddev as Cell3D
from CHARMtools import MultiCell3D

metadata = pd.read_csv("./all.metadata.tsv", sep="\t",header=None)
metadata.columns = ["cellname", "celltype","mc_15","mc_20","mc_25"]

with open("brain.pkl", "rb") as f:
    brain = pickle.load(f)
cells = list(brain.cells_dict.values())

#variable_genes = pd.read_csv("./all.top5k.variablegenes.tsv",header=None).values.T.astype(str)[0][:5000]
variable_genes = pd.read_table("./groud_truth/considered_genes.txt", header=None).values.T.astype(str)[0]

tss = pd.read_csv("./tss.bed",sep="\t",header=None)
tss.columns = ["chrom","start","end","gene","strand"]
tss["genome_coord"] = tss.apply(lambda x: f"{x['chrom']}:{x['start']//5000*5000-int(2e6)}-{x['start']//5000*5000+5000+int(2e6)}", axis=1)
tss = tss.query('gene in @variable_genes').query('chrom != "chrM"')[["gene","genome_coord"]].drop_duplicates().reset_index(drop=True)
tss["tss_id"] = tss["gene"] + "_" + (tss.groupby("gene").cumcount() + 1).astype(str)

with open("brain.pkl", "rb") as f:
    brain = pickle.load(f)
cells = list(brain.cells_dict.values())

def _get_data(gene):
    try:
        genome_coord = tss.query('gene == @gene')["genome_coord"].values[0]

        def _process_cell(cell):
            temp_tdgs = []
            chrom = genome_coord.split(":")[0]
            
            for suffix in ['a:', 'b:']:
                temp_tdg = cell.get_data(genome_coord.replace(":", suffix)).reset_index(drop=True).copy()
                tss_row = temp_tdg.iloc[len(temp_tdg)//2]
                tss_loc = tss_row[["x","y","z"]].values
                
                temp_tdg["distance"] = np.sqrt((temp_tdg['x'] - tss_loc[0])**2 + (temp_tdg['y'] - tss_loc[1])**2 + (temp_tdg['z'] - tss_loc[2])**2)
                temp_tdg["linear_bin"] = np.abs(temp_tdg['pos'] - tss_row["pos"]) // 5000
                temp_tdg["obsexp_distance"] = temp_tdg["distance"] / cell.expected[f"{chrom}{suffix[0]}"][temp_tdg["linear_bin"].astype(int)]
                temp_tdg["distance_to_tss"] = temp_tdg['pos'] - tss_row["pos"]
                temp_tdg["cellname"] = cell.cellname
                
                temp_tdgs.append(temp_tdg[["chrom","pos","obsexp_distance","distance","atac","ct","distance_to_tss","cellname"]])
            
            return temp_tdgs
        
        results = [_process_cell(cell) for cell in cells]
        temp_tdgs = [item for sublist in results for item in sublist]
        concat_df = pd.concat(temp_tdgs)
        concat_df['chrom'] = concat_df['chrom'].str.replace('[ab]', '', regex=True)

        concat_df_mean = concat_df.groupby(["cellname","chrom","pos"]).mean()
        concat_df_mean["obsexp_contact"] = 1 / concat_df_mean["obsexp_distance"]
        concat_df_mean.reset_index(inplace=True)
        concat_df_mean.fillna(0,inplace=True)

        concat_df_mean.to_csv(f"data/{gene}.tsv.gz",sep="\t",compression="gzip",index=False)
        return 0
    except Exception as e:
        return 1


with concurrent.futures.ProcessPoolExecutor(max_workers=200) as executor:
    results = list(tqdm.tqdm(executor.map(_get_data, variable_genes), total=len(variable_genes)))

print(f"Failed: {sum(results)}")    
print("All done!")
