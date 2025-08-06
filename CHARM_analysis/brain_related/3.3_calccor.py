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

sys.path.append('/shareb/zliu/analysis/')
sys.path.append('/shareb/zliu/analysis/CHARMtools')
from upload.CHARMtools.archieve import Cell3Ddev as Cell3D
from CHARMtools import MultiCell3D

import matplotlib.pyplot as plt
import seaborn as sns
import tqdm

from scipy import stats
import statsmodels

np.random.seed(42)
import pickle
import warnings
from sklearn.cluster import KMeans

import pybedtools
import pickle


metadata = pd.read_csv("./all.metadata.tsv", sep="\t",header=None)
metadata.columns = ["cellname", "celltype","mc_id","mc_20","mc_25"]
metadata_mc = metadata.groupby("mc_id")["celltype"].agg(lambda x: x.mode()[0]).reset_index()

rnamat = pd.read_csv("SCT_metacellmean.tsv.gz",sep="\t",index_col=0)
rnamat.index.name = "mc_id"
variable_genes = pd.read_table("./groud_truth/considered_genes.txt", header=None).values.T.astype(str)[0]

with open("brain.pkl", "rb") as f:
    brain = pickle.load(f)
cells = list(brain.cells_dict.values())


tss = pd.read_csv("tss.bed",sep="\t",header=None)
tss.columns = ["chrom","start","end","gene","strand"]
tss["genome_coord"] = tss.apply(lambda x: f"{x['chrom']}:{x['start']//5000*5000-int(2e6)}-{x['start']//5000*5000+5000+int(2e6)}", axis=1)
tss = tss.query('chrom != "chrM"')[["gene","genome_coord"]].drop_duplicates().reset_index(drop=True)
tss["tss_id"] = tss["gene"] + "_" + (tss.groupby("gene").cumcount() + 1).astype(str)

print("Finished loading data.")

def calccor(tss_id):
    # suppress warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            gene,genome_coord = tss.query('tss_id == @tss_id')[["gene","genome_coord"]].values[0]

            temp_tdgs = []
            for cell in cells:
                chrom = genome_coord.split(":")[0]
                temp_tdg = cell.get_data(genome_coord.replace(":","a:")).reset_index(drop=True).copy()
                tss_row = temp_tdg.iloc[temp_tdg.index.values[len(temp_tdg.index.values)//2],:]
                tss_location = tss_row[["x","y","z"]].values
                temp_tdg["distance"] = np.sqrt((temp_tdg['x'] - tss_location[0])**2 + (temp_tdg['y'] - tss_location[1])**2 + (temp_tdg['z'] - tss_location[2])**2)
                temp_tdg["linear_bin"] = np.abs(temp_tdg['pos'] - tss_row["pos"]) // 5000
                temp_tdg["obsexp_distance"] = temp_tdg["distance"] / cell.expected[chrom + "a"][temp_tdg["linear_bin"].values.astype('int')]
                temp_tdg["distance_to_tss"] = temp_tdg['pos'] - tss_row["pos"]
                temp_tdg["cellname"] = cell.cellname
                temp_tdgs.append(temp_tdg[["chrom","pos","obsexp_distance","distance","atac","ct","distance_to_tss","cellname"]])

                temp_tdg = cell.get_data(genome_coord.replace(":","b:")).reset_index(drop=True).copy()
                tss_row = temp_tdg.iloc[temp_tdg.index.values[len(temp_tdg.index.values)//2],:]
                tss_location = tss_row[["x","y","z"]].values
                temp_tdg["distance"] = np.sqrt((temp_tdg['x'] - tss_location[0])**2 + (temp_tdg['y'] - tss_location[1])**2 + (temp_tdg['z'] - tss_location[2])**2)
                temp_tdg["linear_bin"] = np.abs(temp_tdg['pos'] - tss_row["pos"]) // 5000
                temp_tdg["obsexp_distance"] = temp_tdg["distance"] / cell.expected[chrom + "b"][temp_tdg["linear_bin"].values.astype('int')]
                temp_tdg["distance_to_tss"] = temp_tdg['pos'] - tss_row["pos"]
                temp_tdg["cellname"] = cell.cellname
                temp_tdgs.append(temp_tdg[["chrom","pos","obsexp_distance","distance","atac","ct","distance_to_tss","cellname"]])

            concat_df = pd.concat(temp_tdgs)
            concat_df['chrom'] = concat_df['chrom'].apply(lambda x: re.sub('[ab]','',x))
            summarise_df = concat_df.merge(metadata,how="left").groupby(['mc_id','distance_to_tss'])[["obsexp_distance","distance","atac","ct"]].mean()
            df_reset = summarise_df.reset_index()
            merged_df = df_reset.merge(rnamat[[gene]], left_on='mc_id', right_index=True)
            merged_df = merged_df.merge(metadata_mc, on='mc_id', how='left')

            temp_rna = rnamat[[gene]].copy()
            kmeans = KMeans(n_clusters=2, random_state=42).fit(rnamat[[gene]].values)
            temp_rna["group"] = kmeans.labels_

            index_order = temp_rna.groupby("group").mean().sort_values(gene).index.values
            if index_order[0] == 1:
                temp_rna["group"] = temp_rna["group"].replace({0:"high",1:"low"})
            else:
                temp_rna["group"] = temp_rna["group"].replace({1:"high",0:"low"})
            high_mc_id = temp_rna.query('group == "high"').index.values
            low_mc_id = temp_rna.query('group == "low"').index.values

            # temp_rna = rnamat[[gene]].copy()
            # temp_rna = temp_rna.merge(metadata_mc, left_index=True, right_on='mc_id', how='left')
            # temp_rna = temp_rna.groupby('celltype')[gene].mean().reset_index()
            # temp_rna = temp_rna.sort_values(by=gene, ascending=False)
            # high_celltype = temp_rna.iloc[0]['celltype']
            # high_mc_id = metadata_mc.query('celltype == @high_celltype')['mc_id'].values
            # low_celltype = temp_rna.iloc[-1]['celltype']
            # low_mc_id = metadata_mc.query('celltype == @low_celltype')['mc_id'].values

            def calculate_stats(group):
                stats_dict = {}
                for col in ['obsexp_distance','distance', 'atac','ct']:
                    valid_data = group[[col, gene]].dropna()
                    x = valid_data[col]
                    y = valid_data[gene]
                    if len(x) >= 2:
                        r, p = stats.pearsonr(x, y)
                        #r, p = stats.spearmanr(x, y)
                        stats_dict[f"{col}_corr"] = r
                        stats_dict[f"{col}_pval"] = p
                    else:
                        stats_dict[f"{col}_corr"] = 0
                        stats_dict[f"{col}_pval"] = 1
                        
                return pd.Series(stats_dict)

            correlations = merged_df.groupby('distance_to_tss').apply(calculate_stats)
            correlations['obsexp_distance_corr']=correlations['obsexp_distance_corr'].fillna(0)
            correlations['distance_corr']=correlations['distance_corr'].fillna(0)
            correlations['atac_corr']=correlations['atac_corr'].fillna(0)
            correlations['ct_corr']=correlations['ct_corr'].fillna(0)

            correlations['obsexp_distance_pval']=correlations['obsexp_distance_pval'].fillna(1)
            correlations['distance_pval']=correlations['distance_pval'].fillna(1)
            correlations['atac_pval']=correlations['atac_pval'].fillna(1)
            correlations['ct_pval']=correlations['ct_pval'].fillna(1)

            #high_exp = merged_df.query('mc_id in @high_mc_id').groupby('distance_to_tss')[["obsexp_distance","distance"]].mean()
            high_exp = merged_df.query('mc_id in @high_mc_id')\
                    .groupby(['celltype', 'distance_to_tss'])[["obsexp_distance", "distance"]]\
                    .mean()\
                    .reset_index()
            high_exp = high_exp.sort_values('distance', ascending=True)
            high_exp = high_exp.drop_duplicates(subset='distance_to_tss', keep='first')
            high_exp = high_exp.set_index('distance_to_tss')[['obsexp_distance', 'distance']].sort_index()
            low_exp = merged_df.query('mc_id in @low_mc_id').groupby('distance_to_tss')[["obsexp_distance","distance"]].mean()
            high_exp.columns = ['obsexp_distance_highexp','distance_highexp']
            low_exp.columns = ['obsexp_distance_lowexp','distance_lowexp']

            result = pd.concat([correlations, high_exp, low_exp], axis=1)

            result.to_csv(f"./corres/{tss_id}.tsv.gz",sep="\t",index=False)
            return None
    except Exception as e:
        print(f"Error processing {tss_id}: {e}")
        return None

tss_ids = variable_genes
tss_ids = [tss_id + "_1" for tss_id in tss_ids]

# with concurrent.futures.ProcessPoolExecutor(200,mp_context=multiprocessing.get_context("fork")) as executor:
#     list(tqdm.tqdm(executor.map(calccor, tss["tss_id"].sample(frac=1)), total=tss.shape[0]))

with concurrent.futures.ProcessPoolExecutor(200,mp_context=multiprocessing.get_context("fork")) as executor:
   list(tqdm.tqdm(executor.map(calccor, tss_ids), total=len(tss_ids)))
