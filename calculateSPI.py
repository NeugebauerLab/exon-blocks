##################################################################################
## Import required packages 
##################################################################################

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, combine_pvalues, false_discovery_control
import seaborn as sns
import argparse


##########################################################################
## Define a function to format the readthrough cutoff ranges input
##########################################################################

def list_of_strings(arg):
    return list(map(str, arg.split(',')))


##################################################################################
## Main code
##################################################################################

def run(output_prefix, comparisons_condition1, comparisons_condition2, p_mask_cutoff, path_to_files, spliceq_files, samples, treatments):

	# Set the index for the first sample 
	idx = 0

	# Set the sample name for the first sample 
	sample_name = samples[idx]

	# Load the first sample file 
	sample_df = pd.read_csv(spliceq_files[idx], delimiter = '\t', low_memory=False)

	# Calculate the spliced counts for each intron entry
	sample_df[sample_name] = sample_df['sj5_cov_split'] + sample_df['sj3_cov_split']

	# Keep only a subset of columns 
	sample_df_subset = sample_df[['transcript_ID', 'intron_ID', sample_name]]

	# Set the final_spliced_df equal to the subset of columns from the first sample 
	final_spliced_df = sample_df_subset

	# Check if there are remaining samples 
	if len(samples) > 1:
	    # Loop through remaining samples 
	    for idx in range(0, len(samples)):
	        # Set the sample name for the current sample 
	        sample_name = samples[idx]
	        
	        # Load the current sample file 
	        sample_df = pd.read_csv(spliceq_files[idx], delimiter = '\t', low_memory=False)
	        
	        # Calculate the spliced counts for each intron entry 
	        sample_df[sample_name] = sample_df['sj5_cov_split'] + sample_df['sj3_cov_split']

	        # Keep only a subset of columns 
	        sample_df_subset = sample_df[['transcript_ID', 'intron_ID', sample_name]]

	        # Add the subset df for the current sample to the final df 
	        final_spliced_df = final_spliced_df.merge(sample_df_subset)


	final_spliced_df.to_csv('final_spliced_df.tsv', sep = '\t', header = True, index = False)

	# Set the index for the first sample 
	idx = 0

	# Set the sample name for the first sample 
	sample_name = samples[idx]

	# Load the first sample file 
	sample_df = pd.read_csv(spliceq_files[idx], delimiter = '\t', low_memory=False)

	# Calculate the total counts for each intron entry
	sample_df[sample_name] = sample_df['sj5_cov_split'] + sample_df['sj3_cov_split'] + sample_df['sj5_cov_nonsplit'] + sample_df['sj3_cov_nonsplit']

	# Keep only a subset of columns 
	sample_df_subset = sample_df[['transcript_ID', 'intron_ID', sample_name]]

	# Set the final_total_df equal to the subset of columns from the first sample 
	final_total_df = sample_df_subset

	# Check if there are remaining samples 
	if len(samples) > 1:
	    # Loop through remaining samples 
	    for idx in range(0, len(samples)):
	        
	        # Set the sample name for the current sample 
	        sample_name = samples[idx]
	        
	        # Load the current sample file 
	        sample_df = pd.read_csv(spliceq_files[idx], delimiter = '\t', low_memory=False)
	        
	        # Calculate the total counts for each intron entry 
	        sample_df[sample_name] = sample_df['sj5_cov_split'] + sample_df['sj3_cov_split'] + sample_df['sj5_cov_nonsplit'] + sample_df['sj3_cov_nonsplit']

	        # Keep only a subset of columns 
	        sample_df_subset = sample_df[['transcript_ID', 'intron_ID', sample_name]]

	        # Add the subset df for the current sample to the final df 
	        final_total_df = final_total_df.merge(sample_df_subset)


	final_total_df.to_csv('final_total_df.tsv', sep = '\t', header = True, index = False)

	for idx_spi in range(0, len(comparisons_condition1)):

	    # Define condition 1
	    condition1 = comparisons_condition1[idx_spi]

	    # Define condition 2
	    condition2 = comparisons_condition2[idx_spi]

	    # Find the indices of treatments that match condition 1
	    condition1_idx = [idx for idx, treatment in enumerate(treatments) if treatment == condition1]

	    # Find the indices of treatments that match condition 2
	    condition2_idx = [idx for idx, treatment in enumerate(treatments) if treatment == condition2]

	    # Convert the list of samples to a numpy array 
	    samples_numpy = np.array(samples)

	    # # Select the samples of interest
	    target_samples = list(np.concatenate((samples_numpy[condition1_idx], samples_numpy[condition2_idx])))

	    # Create a matrix of total counts 
	    total_mat = final_total_df[target_samples].to_numpy()

	    # Create a matrix of spliced counts 
	    spliced_mat = final_spliced_df[target_samples].to_numpy()

	    # Create a matrix of unspliced counts 
	    unspliced_mat = total_mat - spliced_mat

	    # Create a matrix of SPI values 
	    # se_mat = spliced_mat / total_mat
	    se_mat = spliced_mat / total_mat

	    # Initialize an empty list to store p-values 
	    pvals = [] 
	    
	    number_of_replicates = len(target_samples) / 2

	    if number_of_replicates == 3:
	        
	        # Calculate p-values 
	        for i in range(len(total_mat)):
	            cont_table_1 = [[unspliced_mat[i,0], spliced_mat[i,0]], [unspliced_mat[i,3], spliced_mat[i,3]]]
	            cont_table_2 = [[unspliced_mat[i,1], spliced_mat[i,1]], [unspliced_mat[i,4], spliced_mat[i,4]]]
	            cont_table_3 = [[unspliced_mat[i,2], spliced_mat[i,2]], [unspliced_mat[i,5], spliced_mat[i,5]]]

	            try:
	                p1 = chi2_contingency(cont_table_1).pvalue
	            except:
	                p1 = np.nan

	            try:
	                p2 = chi2_contingency(cont_table_2).pvalue
	            except:
	                p2 = np.nan

	            try:
	                p3 = chi2_contingency(cont_table_3).pvalue
	            except:
	                p3 = np.nan

	            if all(np.isfinite([p1, p2, p3])):
	                p = combine_pvalues([p1, p2, p3]).pvalue
	            else:
	                p = np.nan

	            pvals.append(p)
	            
	    elif number_of_replicates == 2:
	        
	        # Calculate p-values
	        for i in range(len(total_mat)):
		        cont_table_1 = [[unspliced_mat[i,0], spliced_mat[i,0]], [unspliced_mat[i,2], spliced_mat[i,2]]]
		        cont_table_2 = [[unspliced_mat[i,1], spliced_mat[i,1]], [unspliced_mat[i,3], spliced_mat[i,3]]]

		        try:
		            p1 = chi2_contingency(cont_table_1).pvalue
		        except:
		            p1 = np.nan

		        try:
		            p2 = chi2_contingency(cont_table_2).pvalue
		        except:
		            p2 = np.nan

		        if all(np.isfinite([p1, p2])):
		            p = combine_pvalues([p1, p2]).pvalue
		        else:
		            p = np.nan

		        pvals.append(p)

	    # Convert pvals to numpy array 
	    pvals = np.array(pvals)

	    if number_of_replicates == 3:

	        # Create a mask to distinguish np.nan values from actual values 
	        p_mask = np.isfinite(pvals) & (np.sum(total_mat[:, :3], axis = 1) >= p_mask_cutoff) & (np.sum(total_mat[:, 3:6], axis = 1) >= p_mask_cutoff) # mask is a list of falses and trues 

	    elif number_of_replicates == 2:

	    	# Create a mask to distinguish np.nan values from actual values 
	        p_mask = np.isfinite(pvals) & (np.sum(total_mat[:, :2], axis = 1) >= p_mask_cutoff) & (np.sum(total_mat[:, 2:4], axis = 1) >= p_mask_cutoff) # mask is a list of falses and trues 

	    # Calculate adjusted p-values 
	    pvals_adj = false_discovery_control(pvals[p_mask])

	    if number_of_replicates == 3:

	        # Calculate SPI for condition 1
	        se_C = np.sum(spliced_mat[:, :3], axis = 1) / np.sum(total_mat[:, :3], axis = 1)

	        # Calculate SPI for condition 2
	        se_T = np.sum(spliced_mat[:, 3:6], axis = 1) / np.sum(total_mat[:, 3:6], axis = 1)

	    elif number_of_replicates == 2:

	        # Calculate SPI for condition 1
	        se_C = np.sum(spliced_mat[:, :2], axis = 1) / np.sum(total_mat[:, :2], axis = 1)

	        # Calculate SPI for condition 2
	        se_T = np.sum(spliced_mat[:, 2:4], axis = 1) / np.sum(total_mat[:, 2:4], axis = 1)

	    # Calculate log10 of adjusted p-values 
	    y = -np.log10(pvals_adj)

	    # Filter the SPI values to remove any np.nan values 
	    # x = (se_C - se_T)[p_mask]
	    x = (se_T - se_C)[p_mask]

	    # Classify as significant or not significant based on adjusted p-value 
	    z = ['Not Significant' if pval > 0.01 else 'Significant' for pval in pvals_adj]

	    # Generate a volcano plot 
	    fig = plt.figure(figsize = (3,3))
	    plt.scatter(x[pvals_adj > 0.01], y[pvals_adj > 0.01], 20, color = "gray", alpha = 0.2, edgecolors = "none")
	    plt.scatter(x[pvals_adj <= 0.01], y[pvals_adj <= 0.01], 20, color = "red", alpha = 0.2, edgecolors = "none")
	    plt.xlabel('Splicing Difference')
	    plt.ylabel('-log10(p-value)')
	    plt.title(f'Significant Introns: {np.sum(pvals_adj <= 0.01)}')
	    plt.show() 
	    fig.savefig(f'{output_prefix}_{condition1}_vs_{condition2}_SPI_volcanoPlot.png', bbox_inches='tight', dpi = 300)

	    # Generate histogram of adjusted p-values 
	    fig = plt.figure(figsize=(4, 3))
	    bins = np.arange(0, 1, 0.05)
	    plt.hist(pvals_adj, bins=bins)
	    plt.xlabel('Adjusted p-value')
	    plt.ylabel('Number of Introns')
	    plt.title(f'{condition1} vs {condition2}')
	    plt.show()
	    fig.savefig(f'{output_prefix}_{condition1}_vs_{condition2}_histogram_adjustedPValues_for_SPI.png', bbox_inches='tight', dpi = 300)

	    # Generate output table for R 
	    pd.DataFrame({'transcript_ID':final_total_df['transcript_ID'][p_mask], 'intron_ID':final_total_df['intron_ID'][p_mask], 'SPI_diff':x, 'adjusted_pvalue':y, 'Significance':z}).to_csv(f'{output_prefix}_{condition1}_vs_{condition2}_CoSE_adjustedpval.tsv', sep = '\t', header = True, index = False)

	    # Generate output table of raw SPI scores
	    se_df_out = pd.DataFrame(se_mat)
	    se_df_out.columns = final_spliced_df[target_samples].columns
	    se_df_out['transcript_ID'] = final_spliced_df['transcript_ID']
	    se_df_out['intron_ID'] = final_spliced_df['intron_ID']
	    se_df_out.to_csv(f'{output_prefix}_{condition1}_vs_{condition2}_RAW_SPI.tsv', sep = '\t', header = True, index = False)

	    # Generate output table of total counts
	    total_df_out = pd.DataFrame(total_mat)
	    total_df_out.columns = final_total_df[target_samples].columns
	    total_df_out['transcript_ID'] = final_total_df['transcript_ID']
	    total_df_out['intron_ID'] = final_total_df['intron_ID']
	    total_df_out.to_csv(f'{output_prefix}_{condition1}_vs_{condition2}_totalCount_RAW_SPI.tsv', sep = '\t', header = True, index = False)

	    # Generate output table with p-values 
	    output_df = pd.DataFrame(se_mat)
	    
	    if number_of_replicates == 2:
	        
	        output_df = output_df.rename(columns = {0:target_samples[0], 1:target_samples[1], 2:target_samples[2], 3:target_samples[3]})
	        
	    elif number_of_replicates == 3:
	        
	        output_df = output_df.rename(columns = {0:target_samples[0], 1:target_samples[1], 2:target_samples[2], 3:target_samples[3], 4:target_samples[4], 5:target_samples[5]})
	    
	    output_df = output_df[p_mask]
	    output_df['SPI_diff'] = x
	    output_df['pvals_adj'] = pvals_adj
	    output_df['transcript_ID'] = final_spliced_df['transcript_ID'][p_mask]
	    output_df['intron_ID'] = final_spliced_df['intron_ID'][p_mask]
	    output_df_sort = output_df.sort_values('pvals_adj')
	    output_df_sort.to_csv(f'{output_prefix}_SPI_for_{condition1}_vs_{condition2}.tsv', sep = '\t', header = True, index = False)


if __name__ == '__main__':

    ## Define inputs

    # Create a parser
    parser = argparse.ArgumentParser(description = "Calculate Splicing Per Intron (SPI)")

    # Add arguments to parser
    parser.add_argument("-out", required = True, help = "output prefix")
    parser.add_argument("-condition1", type = list_of_strings, required = True, default = False, help = "a list of condition 1 for each comparison")
    parser.add_argument("-condition2", type = list_of_strings, required = True, default = False, help = "a list of condition 2 for each comparison")
    parser.add_argument("-pMaskCutoff", type = int, required = False, default = 0, help = 'minimum number of total reads for a p-value to be calculated for an intron CoSE value')
    parser.add_argument("-m", required = True, help = "tab-separated manifest file containing filtered reads file names, samples, and treatments (each in their own column, in that order)")
    parser.add_argument("-filePath", type = str, required = False, default = '', help = 'path to directory containing the sample and manifest file')

    # Create arguments from parser
    args = parser.parse_args()

    output_prefix = args.out
    comparisons_condition1 = args.condition1
    comparisons_condition2 = args.condition2
    p_mask_cutoff = args.pMaskCutoff

    path_to_files = args.filePath

    input_file = args.m

    input_file_withPath = path_to_files + input_file

    manifest = pd.read_csv(input_file_withPath, sep = '\t', header = None)
    spliceq_files = list(manifest[0])

    spliceq_files = [path_to_files + file for file in spliceq_files]

    samples = list(manifest[1])
    treatments = list(manifest[2])

    run(output_prefix, comparisons_condition1, comparisons_condition2, p_mask_cutoff, path_to_files, spliceq_files, samples, treatments)

#python calculateSPI_25June2024.py \
# -out eif4a3_nSmash_SPI \
# -condition1 ctr \
# -condition2 kd \
# -m spliceq_nSMASH_manifest.tsv \
# -filePath /home/ab3525/ycga_work/SRS_nsmash/spliceq/




















    

