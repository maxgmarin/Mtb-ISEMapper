#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

def Depth_TSV_To_Peaks(in_Depth_TSV, out_Peaks_TXT, genome_length):

    RvDepth_data = pd.read_csv(in_Depth_TSV, sep ="\t", header=None)
    RvDepth_data.columns = ["chrom", "position", "depth"]

    

    # Step 1: Initialize a NumPy array the length of the into genome
    depth_array = np.zeros(genome_length)

    # Step 2: Update depth values in the depth_array to match the input depth data
    pos_to_update = RvDepth_data['position'] - 1  # Convert positions to zero-based indexing
    depth_to_update = RvDepth_data['depth']
    np.put(depth_array, pos_to_update, depth_to_update)

    # Step 3: Find peaks of depth, indicating insertion points of IS6110
    peaks, _ = find_peaks(depth_array, distance = 1500, prominence = 20)

    print ("Number of IS6110 Insertions detected of IS6110: ", len(peaks))
    print (peaks)

    peaks_list = list(peaks)

    fileout1 = open(out_Peaks_TXT, 'w')
    for p in peaks_list:
        fileout1.write(str(p) + "\n")
    fileout1.close()  



if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description='Identify IS6110 insertion points based on peaks of depth')

    # Add arguments
    parser.add_argument('-i', '--in_depth_tsv', type=str, help='The path to the input depth file')
    parser.add_argument('-o', '--out_peaks', type=str, help='The path to the output peaks file')

    # Parse the arguments
    args = parser.parse_args()

    # Run the main function    
    Depth_TSV_To_Peaks(args.in_depth_tsv, args.out_peaks, 4411532)

