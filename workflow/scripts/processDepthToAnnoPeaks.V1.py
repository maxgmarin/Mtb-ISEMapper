#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import bioframe as bf


def parse_RvDepth_To_NP(in_Depth_TSV, genome_length=4411532):
    RvDepth_data = pd.read_csv(in_Depth_TSV, sep ="\t", header=None)
    RvDepth_data.columns = ["chrom", "position", "depth"]
    
    # Step 1: Initialize a NumPy array the length of the into genome
    depth_array = np.zeros(genome_length)
    
    # Step 2: Update depth values in the depth_array to match the input depth data
    pos_to_update = RvDepth_data['position'] - 1  # Convert positions to zero-based indexing
    depth_to_update = RvDepth_data['depth']
    np.put(depth_array, pos_to_update, depth_to_update)

    return depth_array

def Anno_ISPeaks(i_Peaks_DF, H37Rv_Anno_DF):

    i_ISPeaks_Anno_DF = bf.overlap(i_Peaks_DF, H37Rv_Anno_DF)

    i_ISPeaks_Anno_DF = i_ISPeaks_Anno_DF.drop(["end", "chrom_", "start_", "end_", "strand_", "feature_", "functional_category_"], axis = 1)
    i_ISPeaks_Anno_DF.columns = ["chrom", "Pos", "Rv_GeneID", "symbol"]
    i_ISPeaks_Anno_DF['symbol'] = i_ISPeaks_Anno_DF['symbol'].fillna(i_ISPeaks_Anno_DF['Rv_GeneID'])
    return i_ISPeaks_Anno_DF

def parse_H37RvAnnoTSV(i_Rv_Regions_TSVGZ):
    RvAnno_DF = pd.read_csv(i_Rv_Regions_TSVGZ, sep = "\t")
    RvAnno_DF.columns = ['chrom', 'start', 'end', 'strand', 'H37rv_GeneID', 'symbol', 'feature',
           'functional_category']
    return RvAnno_DF

def DepthToAnnoPeaks(i_depth_array, H37Rv_Anno_DF):

    # Step 1: Find peaks of depth, indicating insertion points of IS6110
    peaks, _ = find_peaks(i_depth_array, distance = 1500, prominence = 20)

    print ("Number of IS6110 Insertions detected of IS6110: ", len(peaks))
    
    # Step 2: Convert peaks to a DF (BED-like)
    peaks_DF = pd.DataFrame(peaks, columns=['position'])
    peaks_DF['start'] = peaks_DF['position'].astype(int) - 1
    peaks_DF['chrom'] = 'NC_000962.3'
    peaks_DF['end'] = peaks_DF['start'] + 1
    peaks_DF = peaks_DF[['chrom', 'start', 'end']]

    # Step 3: Annotate peaks DF by overlapping w/ H37Rv annotations DF
    i_ISPeaks_Anno_DF = Anno_ISPeaks(i_Peaks_DF, H37Rv_Anno_DF)
    
    return i_ISPeaks_Anno_DF



if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description='Identify IS6110 insertion points based on peaks of depth')

    # Add arguments
    parser.add_argument('-i', '--in_depth_tsv', type=str, help='Path to the IS6110 alignment depth tsv (from samtools depth)')
    parser.add_argument('-a', '--h37rv_anno_tsv', type=str, help='Path to a TSV containing H37Rv annotations (Genes + intergenic regions)')     
    parser.add_argument('-o', '--out_peaks', type=str, help='The path to the annotated output peaks file')

    # Parse the arguments
    args = parser.parse_args()

    IS6110_Depth_NP = parse_RvDepth_To_NP(args.in_depth_tsv, 4411532)

    H37Rv_Anno_DF = parse_H37RvAnnoTSV(args.h37rv_anno_tsv)

    Peaks_Anno_DF = DepthToAnnoPeaks(IS6110_Depth_NP, H37Rv_Anno_DF)

    Peaks_Anno_DF.to_csv(args.out_peaks, sep="\t", index=False)


