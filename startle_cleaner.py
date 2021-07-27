#!/bin/python
'''
Author: David Gruskin (david.gruskin@columbia.edu)
Function of this script: Clean startle box output files
Inputs:
    - input_file [path to file (or directory of files) to be cleaned]
    - output_file [path to where output should be saved]

Outputs:
    - output [cleaned csv file]
Last updated: 07/2021
'''

'''
When "and" for every trial, then odds are always unpaired, evens always paired
When not "and" for every trial, if and, paired


'''




# Import some stuff
import numpy as np
import pandas as pd
import string
import argparse
import re
import glob
import os

def main():

    # Initialize variables here:
    parser = argparse.ArgumentParser(description='Let''s clean some files!')
    parser.add_argument("-input_file", "--input_file", required=True, type=str,default = 'None')
    parser.add_argument("-output_file", "--output_file", required=True, type=str,default = 'None')
    args = vars(parser.parse_args())

    # Get the filename from the first argument
    filename1 = args['input_file']

    # Set the column IDs
    columns_list     = ['chamber','null','white_noise_block_1','stim2_unpaired_latency','stim2_unpaired_peak_time','stim2_unpaired_peak_value','stim2_unpaired_duration','stim2_unpaired_total','stim2_paired_latency','stim2_paired_peak_time','stim2_paired_peak_value','stim2_paired_duration','stim2_paired_total']
    big_columns_list = ['dam','group','litter','sex','mouse_number'] + columns_list
    
    # Initialize dataframes
    output = pd.DataFrame(columns = big_columns_list) 
    all_trials = pd.DataFrame(columns = big_columns_list) 

    # If the filename refers to a directory of files, grab all filenames. Otherwise, just grab input.
    filename_list = []
    if os.path.isdir(filename1):
        multiple_files = 1
        os.chdir(filename1)
        for filename in glob.glob("*Computed.txt"):
            filename_list.append(filename)
    else:
        filename_list.append(filename1)

    # Find all filenames that end in "Computed.txt"
    for filename in filename_list:
        print(filename)

        # Open the file and read the text
        file1    = open(filename,"r")
        Lines    = file1.readlines()
        num_dams = re.findall(r'\d{3}(?!\d)', args['input_file'])

        # Initialize lists
        group_list   = []
        litter_list  = []
        chamber_list = []
        sex_list     = []
        mouse_list   = []
        dam_list     = []

        # Break down the filename to find import information and store that information in lists
        counter = 0
        for x in filename.split('/')[-1].split('_')[:-2]:
            if len(re.findall(r'\d{3}(?!\d)', x)) > 0:
                dam_id = int(x)
                group_id = filename.split('/')[-1].split('_')[counter+1]
                litter_id = filename.split('/')[-1].split('_')[counter+2]
            if "Chamber" in x:
                dam_list.append(dam_id)
                group_list.append(group_id)
                litter_list.append(litter_id)
                chamber_list.append(int(x[-1]))
                sex_list.append(filename.split('/')[-1].split('_')[counter-1][0])
                mouse_list.append(filename.split('/')[-1].split('_')[counter-1][1])
            counter +=1

        # Convert text into csv by separating at new lines
        Lines    = list(filter(lambda x: x != '\n', Lines)) 
        df_input = pd.DataFrame([x.split('\t') for x in Lines])

        #df_input.to_csv(args['output_file'] + '_input'+ '.csv',index=False)
        
        # Loop through chambers
        counter = 0
        for chamber in chamber_list:

            # Get locations of block 1 and block 2 data for current chamber
            block1_start = df_input.index[df_input[0].str.contains(f"Chamber: {str(chamber)}  Block: 1")].to_list()[0]
            block2_start = df_input.index[df_input[0].str.contains(f"Chamber: {str(chamber)}  Block: 2")].to_list()[0]
            try:
                chamber_end = df_input.index[df_input[0].str.contains(f"Chamber: {str(chamber+1)}  Block: 1")].to_list()[0]
            except:
                chamber_end = df_input.shape[0]

            # Calculate means after finding row indices for variables of interest
            null_ids      = df_input.index[df_input[1].str.contains("Null Period:", na=False)]
            null_ids      = [i for i in null_ids if block1_start<i<block2_start]
            baseline_null = ([float(i) for i in df_input.iloc[null_ids,6].to_list()])[-1]
            nulls         = df_input.iloc[null_ids,6].to_list()[:-1]
            baseline_ids  = [x+2 for x in null_ids]
            baselines     = df_input.iloc[baseline_ids[:-1],6].to_list()

            stim2_ids        = df_input.index[df_input[1].str.contains("Stimulus 2 Period:", na=False)]
            average_id       = df_input.index[df_input[0].str.contains("Average All Trials:", na=False)]
            average_id       = [i for i in average_id if block1_start<i<block2_start]

            # Get row indices for Stimulus 2 info, adjusting for discrepency in indexing for chamber 4
            if chamber == 4:
                stim2_ids_block2 = [i for i in stim2_ids if chamber_end>i>block2_start][:-1]
            else:
                stim2_ids_block2 = [i for i in stim2_ids[:-1] if chamber_end>i>block2_start][:-1]
            and_ids = [x-3 for x in stim2_ids_block2]
            and_count = 0
            counter1 = 0
            and_list = []
            not_and_list = []
            paired_ids = []
            unpaired_ids = []

            # Determine if file reports paired vs unpaired directly or if it must be inferred 
            for string in df_input.iloc[and_ids,0].to_list():
                if "and" in string:
                    and_count +=1
                    and_list.append(counter1)
                else:
                    not_and_list.append(counter1)
                counter1+=1
            if and_count > 11: # only 9 trials should have stimulus 1 ~~and~~ stimulus 2, use 11 in case of edge cases
                file_type = 'every_and'
            else:
                file_type = 'not_every_and'
            paired_ids = [and_ids[i] for i in and_list]
            unpaired_ids = [and_ids[i] for i in not_and_list]

            stim2_ids_block1     = [i for i in stim2_ids[:-1] if block2_start>i>block1_start][:-1]
            baseline_white_noise = float(df_input.iloc[average_id[0]+3,6])
            stim2s               = df_input.iloc[stim2_ids_block1,6].to_list()[:-1]
            print(chamber)

            # Find block 2 stim 2 data unpaired/paired data based on file type
            if file_type == 'every_and':
                unpaired_white_noise = df_input.iloc[stim2_ids_block2[0::2],2:7].astype(float).mean().to_list()
                paired_white_noise   = df_input.iloc[stim2_ids_block2[1::2],2:7].astype(float).mean().to_list()
                all_trials_data = np.concatenate((df_input.iloc[stim2_ids_block2[0::2],2:7].astype(float).values,df_input.iloc[stim2_ids_block2[1::2],2:7].astype(float).values),axis=1)
                all_trials_data = np.column_stack(((np.full((9, 1), dam_list[counter])),(np.full((9, 1), group_list[counter])),(np.full((9, 1), litter_list[counter])),(np.full((9, 1), sex_list[counter])),(np.full((9, 1), mouse_list[counter])),(np.full((9, 1),chamber)),np.array(nulls),np.array(baselines),all_trials_data))

            elif file_type == 'not_every_and':
                unpaired_white_noise = df_input.iloc[[x +3 for x in unpaired_ids],2:7].astype(float).mean().to_list()
                paired_white_noise   = df_input.iloc[[x +3 for x in paired_ids],2:7].astype(float).mean().to_list()
                all_trials_data = np.concatenate((df_input.iloc[[x +3 for x in unpaired_ids],2:7].astype(float).values,df_input.iloc[[x +3 for x in paired_ids],2:7].astype(float).values),axis=1)
                all_trials_data = np.column_stack(((np.full((9, 1), dam_list[counter])),(np.full((9, 1), group_list[counter])),(np.full((9, 1), litter_list[counter])),(np.full((9, 1), sex_list[counter])),(np.full((9, 1), mouse_list[counter])),(np.full((9, 1),chamber)),np.array(nulls),np.array(baselines),all_trials_data))



            # Concatenate variables of interest into numpy array
            # Get the average of the numerical data, along with descriptive data
            all_trials_data_avg = np.column_stack(([all_trials_data[1,0:6]],[all_trials_data[:,6:].astype(np.float64).mean(axis=0)]))
            all_trials.drop_duplicates(inplace=True)
            output.drop_duplicates(inplace=True)

            # Append output from this file to main dataframe
            all_trials = pd.concat([all_trials, pd.DataFrame(data = all_trials_data, columns = big_columns_list)], axis=0)
            data_list  = [chamber]+[baseline_null]+[baseline_white_noise] + unpaired_white_noise +paired_white_noise
            output     = pd.concat([output, pd.DataFrame(data = all_trials_data_avg, columns = big_columns_list)], axis=0)
            counter +=1

    # Save output
    all_trials.to_csv(args['output_file'] + '_trials'+ '.csv',index=False)
    output.to_csv(args['output_file'] + '_average'+ '.csv',index=False)

    print("output saved to: " + args['output_file'] + '_average'+ '.csv' )

if __name__ == "__main__":
    main()
