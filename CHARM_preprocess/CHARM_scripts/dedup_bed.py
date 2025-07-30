import pandas as pd
import numpy as np
import gzip
import os
from sklearn.cluster import DBSCAN
import tqdm
import warnings
warnings.filterwarnings("ignore")

def DBSCAN_wrapper(column,eps):
    try:
        return DBSCAN(eps=eps, min_samples=1).fit(column.values.reshape(-1, 1)).labels_
    except:
        pass

def dedup_bed(rawbed:pd.DataFrame,eps) -> pd.DataFrame:
    """
    dedup pairs by DBSCAN clustering
    """
    raw_bed_num = rawbed.shape[0]
    # try chr1 or chrom1
    rawbed["cluster"] = rawbed.groupby("chrom")["pos"].transform(DBSCAN_wrapper,eps=eps)

    rawbed = rawbed.groupby(["chrom","cluster"]).head(n=1)
    rawbed = rawbed.drop(["cluster","pos"], axis=1)
    dedup_bed_num = rawbed.shape[0]
    rate = 100*(raw_bed_num-dedup_bed_num)/raw_bed_num

    rawbed["start"] = rawbed["start"].astype(int)
    rawbed["end"] = rawbed["end"].astype(int)
    rawbed["score"] = rawbed["score"].astype(int)

    return rawbed, rate

def dedup_wrapper(bed,eps,type="normal"):
    #type in ["normal","spilitpool"]
    if type == "normal":
        bed,rate = dedup_bed(bed,eps)
        print("Duplication rate is %.2f%%" % rate)
    elif type == "splitpool":
        total = bed.shape[0]
        # split dataframe to list of dataframe by column cell
        bed_list = [group[1] for group in bed.groupby("cell")]
        # dedup each dataframe
        dedup_bed_list = [dedup_bed(bed,eps)[0] for bed in tqdm.tqdm(bed_list)]
        # combine dataframe
        bed = pd.concat(dedup_bed_list)
        bed = bed.reset_index(drop=True)
        print("Mean duplication rate is %.2f%%" % ((1 - bed.shape[0]/total)*100))
    else:
        raise ValueError("type should be in ['normal','splitpool']")
    return bed


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", help="input bed file", required=True)
    parser.add_argument("-o","--output", help="output bed file", required=True)
    parser.add_argument("-t","--type", help="dedup type", default="normal", choices=["normal","splitpool"])

    parser.add_argument("-e","--eps",help="eps for DBSCAN", default=1, type=int)
    parser.add_argument("-q","--mapq",help="threshold for mapping quality", default=0, type=int)
    parser.add_argument("--R2_3end",help="dedup R2 3' end", default=True, type=bool)

    args = parser.parse_args()

    atac_bed = pd.read_csv(args.input,sep="\t",header=None)
    atac_bed.columns = ["chrom","start","end","readID","score","strand"]
    # we use 5' end of R2 for dedup
    atac_bed = atac_bed.assign(pos = np.where(atac_bed['strand'] == '+', atac_bed['start'], atac_bed['end']-1)).sort_values(['chrom','pos'])
    print("Running in " + args.type + " mode")
    if args.type == "splitpool":
        atac_bed['cell'] = atac_bed['readID'].str.split('_').str[1]
        atac_bed['readID'] = atac_bed['readID'].str.split('_').str[0]

        atac_bed = atac_bed = atac_bed.query('cell != "unmatched"')
        print("Demultiplexed ATAC reads: " + str(atac_bed.shape[0]))

    eps = args.eps
    dedup_atac_bed = dedup_wrapper(atac_bed,eps,type = args.type)
    if args.R2_3end:
        dedup_atac_bed = dedup_atac_bed.assign(pos = np.where(dedup_atac_bed['strand'] == '-', dedup_atac_bed['start'], dedup_atac_bed['end']-1)).sort_values(['chrom','pos'])
        dedup_atac_bed = dedup_wrapper(dedup_atac_bed,eps,type = args.type)
    dedup_atac_bed = dedup_atac_bed.query('score > @args.mapq')
    print("Deduped ATAC reads after filtering: " + str(dedup_atac_bed.shape[0]))

    dedup_atac_bed.to_csv(args.output,sep="\t",header=None,index=False)
