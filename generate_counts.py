#!/bin/bash
#############################
# File Name : get_motif_counts.py
#
# Purpose : Parse fastq SQUISH output to get counts of codes and motifs
#
# Creation Date : 10-05-2018
#
# Last Modified : Tue 19 Nov 2019 12:03:44 AM PST
#
# Created By : Rob Bierman, Julia Olivieri, and Darren Liu
#
##############################
import argparse
import collections
import numpy as np
import os
import pandas as pd
import sq_utils
import time


def get_tc_and_code(r, rev, motif_len, motif, code_seqs, run):
    """Find target and code within read"""
    left_wind = 0 
    right_wind = len(r) - motif_len
    min_dist = len(motif)
    tc_ind = left_wind
    orig_r = r
    if rev:
        r = sq_utils.rev_comp(r)
    for i in range(left_wind, right_wind + 1):
        if degen_hamming(r[i:i + motif_len], motif) < min_dist:
            min_dist = degen_hamming(r[i:i + motif_len], motif) 
            tc_ind = i
    tc_seq = r[tc_ind:tc_ind + motif_len]
    
    r = orig_r
    if rev:
        code_area = r[:len(r) - tc_ind - motif_len]
    else:
        code_area = sq_utils.rev_comp(r[tc_ind + motif_len:])
        
    max_ind = -1
    code_len = len(code_seqs[0]) - 4
    for code_seq in code_seqs:
        code_seq = code_seq[4:]
        code_ind = code_area.find(code_seq)
        if code_ind != -1 and code_ind > max_ind:
            max_ind = code_ind
            code_len = len(code_seq)
    if max_ind == -1:

        code_seq_r = code_area[-code_len - 4:]
    else:
        code_seq_r = code_area[max_ind - 4:max_ind + code_len]
    return tc_seq, code_seq_r, max_ind, tc_ind


def degen_hamming(s, degen):
    """Find distance between sequences where the second can be degenerate"""
    if len(s) != len(degen):
        raise ValueError("Undefined for sequences of unequal length")
    dist = 0
    for i in range(len(s)):
        same = False
        if s[i] == degen[i]:
            same = True
        elif s[i] == "A" and degen[i] in {"R","M","W","H","D","V","N"}:
            same = True
        elif s[i] == "T" and degen[i] in {"Y","K","W","H","B","D","N"}:
            same = True
        elif s[i] == "G" and degen[i] in {"R","K","S","B","D","V","N"}:
            same = True
        elif s[i] == "C" and degen[i] in {"Y","M","S","H","B","V","N"}:
            same = True
        if not same:
            dist += 1
    return dist


def count_reads(lib,code_dict,spike_dict,prob_dict, code_seqs, make_align):
    """Get target and code count for the given library"""
    parse_style = lib["ParseStyle"]
    counts = collections.defaultdict(int)
    r1_path = os.path.join(data_path,run,lib['R1FileName'])
    r2_path = os.path.join(data_path,run,lib['R2FileName'])
    motif = lib['Motif']
    motif_len = len(motif)
    code_lens = list(set(len(k) for k in code_dict.keys()))

    new_motif = "SASSASASASASAASASASSASASASAASSASASA"
    old_motif = "SASSASASASSAASASASSASASSASSASSA"
    bases = ["A", "T", "C", "G"]
    count = 0
    print_out = False
    for r1,r2 in zip(sq_utils.read_fastq(r1_path),sq_utils.read_fastq(r2_path)):
        count += 1
        if count % 5000 == 0:
            print_out = True
            print(count)
        else:
            print_out = False

        if parse_style == 1:
            if make_align:
                tc_seq_r1, code_seq_r1, code_ind_r1, tc_ind_r1 = get_tc_and_code(r1, False, 
                                         motif_len, motif, code_seqs, run)
                tc_seq_r2, code_seq_r2, code_ind_r2, tc_ind_r2 = get_tc_and_code(r2, True, 
                                         motif_len, motif, code_seqs, run)
            else:
                for code_len in code_lens:
                    code_seq_r1 = r1[motif_len:motif_len+code_len]
                    code_seq_r1 = sq_utils.rev_comp(code_seq_r1)
                    tc_seq_r1 = r1[:motif_len]
       
                    code_seq_r2 = r2[:code_len]
                    tc_seq_r2 = r2[code_len:code_len+motif_len]
                    tc_seq_r2 = sq_utils.rev_comp(tc_seq_r2)
                    if code_seq_r1 == code_seq_r2 and code_seq_r1 in code_dict:
                        break
        elif parse_style == 2:
            if make_align:
                tc_seq_r1, code_seq_r1, code_ind, tc_ind = get_tc_and_code(r1, True, 
                                         motif_len, motif, code_seqs, run)
                tc_seq_r2, code_seq_r2, code_ind, tc_ind = get_tc_and_code(r2, False, 
                                         motif_len, motif, code_seqs, run)
            else:
                tc_seq_r2 = r2[:motif_len]
                motif_ind = r1.find("AGATCGGAAG")
                tc_seq_r1 = sq_utils.rev_comp(r1[motif_ind - motif_len:motif_ind])
                for code_len in code_lens:
                    code_seq_r2 = r2[motif_len:motif_len+code_len]
                    code_seq_r2 = sq_utils.rev_comp(code_seq_r2)
                    if motif_ind - motif_len - code_len < 0:
                        code_seq_r1 = r1[:code_len]
                    else:
                        code_seq_r1 = r1[motif_ind - motif_len - code_len:motif_ind - motif_len]
                    if code_seq_r1 == code_seq_r2 and code_seq_r1 in code_dict:
                        break

        if code_seq_r1 == code_seq_r2 and tc_seq_r1 == tc_seq_r2:

            key = tc_seq_r1+'_'+code_seq_r1
            counts[key] += 1

    rows = {'TargetCode':[],'Count':[],'Spike':[],'Prob':[],
            'TargetPattern':[],'Code':[],'CodeSeq':[]}

    for key,count in counts.items():
        tc_seq,code_seq = key.split('_')
        rows['TargetCode'].append(tc_seq)
        rows['Count'].append(count)
        rows['Spike'].append(spike_dict[tc_seq])
        rows['Prob'].append(prob_dict[tc_seq])
        rows['TargetPattern'].append(sq_utils.contains_motif(tc_seq,motif))
        rows['Code'].append(code_dict[code_seq])
        rows['CodeSeq'].append(code_seq)

    for spike in spike_dict.keys():
        if spike not in rows["TargetCode"]:
            for code_seq in code_seqs:
                rows["TargetCode"].append(spike)
                rows["Count"].append(0)
                rows["Spike"].append(spike_dict[spike])
                rows["Prob"].append(spike_dict[spike])
                rows["TargetPattern"].append(True)
                rows["Code"].append(code_dict[code_seq])
                rows["CodeSeq"].append(code_seq)
    count_df = pd.DataFrame(rows)
    count_df['RunName'] = lib['RunName']
    count_df['LibName'] = lib['LibName']
    count_df['LibNum'] = lib['LibNum']
    count_df = count_df.sort_values(by=['Code','Spike'],ascending=[True,False])
    return count_df
    

def collapse_counts_by_round(count_df):
    """Make only one row for each target per library and run"""
    codes = sorted(count_df['Code'].unique())
    rows = []
    groups = count_df.groupby(['RunName','LibName','TargetCode'])
    for (run,lib,tc_seq),r in groups:
        row = {}
        row['RunName'] = run
        row['LibName'] = lib
        row['LibNum'] = r['LibNum'].iloc[0]
        row['TargetCode'] = tc_seq
        row['TargetPattern'] = r['TargetPattern'].iloc[0]
        row['Spike'] = r['Spike'].iloc[0]
        row['Prob'] = r['Prob'].iloc[0]
        for code in codes:
            row[code] = 0
        for code,code_g in r.groupby(['Code']):
            row[code] = code_g['Count'].sum()
            
        rows.append(row)

    c_df = pd.DataFrame(rows)
    return c_df


def collapse_sub_codes(count_df):
    """Collapse codes with decimals"""
    count_df.head()
    codes = [c for c in count_df.columns if c[:4] == 'Code']
    code_prefix = []
    for c in codes:
        code,sub_code = c.split('.')
        code,num = code.split('Code')
        prefix = 'code'+str(int(num))
        code_prefix.append(prefix)
        count_df[prefix] = 0

    for prefix,code in zip(code_prefix,codes):
        count_df[prefix] += count_df[code]
        del count_df[code]

    return count_df


def add_target_base_counts(count_df):
    """Add counts for each base"""
    for base in ['A','T','C','G']:
        count_df[base+'_count'] = count_df['TargetCode'].str.count(base)
    return count_df


def get_code_dict(lib, run, one_tube):
    """Get concentrations of each code"""
    code_series = all_codes[all_codes['Name'] == lib['CodingSeries']]
    seq_cols = sorted(n for n in code_series.columns if '_Seq' in n)
    code_dict = collections.defaultdict(lambda:'NoCode')
    code_seqs = []
    for seq_col in seq_cols:
        code_seq = code_series.iloc[0][seq_col]
        # in case a space was added in the csv (for readability)
        if not pd.isnull(code_seq):
            code_seq = code_seq.replace(" ", "")

            code_num,_ = seq_col.split('_')

            #This for loop is just for multi-pot codes like "one_tube"
            for i,code_seq in enumerate(code_seq.split('-')):
                code_seqs.append(code_seq)
                for seq in sq_utils.expand_degens(code_seq):
                    if one_tube:
                        code_dict[seq] = code_num+'.'+str(i+1)
                    else:
                        code_dict[seq] = "Code0" + str(max(i + 1,int(code_num[-2:]))) + ".1"
                    if i > 0:
                        print("one tube version: {}\nother version: {}\n".format(code_num+'.'+str(i+1), "Code0" + str(max(i + 1,int(code_num[-2:]))) + ".1"))
        
    return code_dict, code_seqs


def get_spike_dict(lib):
    """Get concentration of each spike"""
    spikes = all_spikes[all_spikes.Name == lib.SpikeSeries]
    seq_fold = collections.defaultdict(lambda:1)
 
    degen_mags = {k:len(v) for k,v in sq_utils.degen_seqs().items()}
    tc_degen = np.prod([degen_mags[b] for b in lib.Motif])
    tc_conc = lib.TargetLibConc

    seq_cols = [n for n in spikes.columns if '_Seq' in n]
    conc_cols = [n for n in spikes.columns if '_Conc' in n]

    min_conc = 0

    if tc_conc == 0:
        for conc_col in conc_cols:
            if min_conc > spikes.iloc[0][conc_col] or min_conc == 0:
                print(spikes.iloc[0][conc_col])
                if spikes.iloc[0][conc_col] != 0:
                    min_conc = spikes.iloc[0][conc_col]
        print("tc_conc: {}".format(tc_conc))
        min_conc = min_conc/10.

    base = tc_conc/tc_degen


    for seq_col,conc_col in zip(seq_cols,conc_cols):
        seq = spikes.iloc[0][seq_col]
        conc = spikes.iloc[0][conc_col]
        if not pd.isnull(seq) and not pd.isnull(conc):
            seq = seq.replace(" ","")
            seq = seq.upper()
            if seq.startswith("GTGCTCTTCCGATCT"):
                seq = seq[15:]

            if tc_conc == 0:
                fold = conc/float(min_conc)
            else:
                fold = (conc/base)+1

            seq_fold[seq] = fold
            print("seq: {}\nfold: {}\n".format(seq,fold))

    return seq_fold
 

def get_prob_dict(lib):
    """Get effective probability of each spike"""
    seq_fold = get_spike_dict(lib)

    degen_mags = {k:len(v) for k,v in sq_utils.degen_seqs().items()}
    tc_degen = np.prod([degen_mags[b] for b in lib.Motif])
    tc_conc = lib.TargetLibConc

    species_base = tc_conc/tc_degen
    base = tc_conc+sum(seq_fold.values())*species_base
    species_prob = species_base/base

    eff_spikes = collections.defaultdict(lambda: species_prob)
    for seq,spike in seq_fold.items():
        eff_spikes[seq] = (spike*species_base)/(base)

    return eff_spikes


#######################################################################
#                                                                     #
#                         Main                                        #
#                                                                     #
####################################################################### 
if __name__ == '__main__': # Run to look at (normally just want a list of 1 run)
    t0 = time.time()
    parser = argparse.ArgumentParser(description="Get code counts from fastq files for SQUICH")
    parser.add_argument("-o","--one_tube",action="store_true",help="If codes are 'one tube' collapse all code counts into Code 1")
    args = parser.parse_args()
    runs = ["20180829_SQUISH14"]

    runs = runs[::-1]

    make_align = True 
    one_tube = args.one_tube 
    if make_align:
        make_align_label = "_align"
    else:
        make_align_label = ""
    if one_tube:
        one_tube_label = ""
    else: 
        one_tube_label = "_newcodelabel"

    # Metadata files
    data_path = '.' 
    all_libs = pd.read_csv('all_libraries.csv')
    all_codes = pd.read_csv('all_codes.csv')
    all_degs = pd.read_csv('all_degs.csv')
    all_spikes = pd.read_csv('all_spikes.csv')


    # Loop through all the runs
    for run in runs:
        libs = all_libs[all_libs.RunName == run]

        # Generate the expanded counts format
        all_counts_path = run+'_expanded_counts.csv'
        all_counts = pd.DataFrame()

        for i,lib in libs.iterrows():

            print('Working on: '+lib['LibName']+'\n')
            code_dict, code_seqs = get_code_dict(lib, run, one_tube)
            spike_dict = get_spike_dict(lib)
            prob_dict = get_prob_dict(lib)
            lib_counts  = count_reads(lib,code_dict,spike_dict,prob_dict, code_seqs, make_align)
            all_counts = pd.concat([all_counts,lib_counts])
            print("LibName: {}".format(lib.LibName))

        rounds = sorted(all_counts['Code'].unique())
        all_counts = add_target_base_counts(all_counts)

        print('Collapsing counts')
        collapsed = collapse_counts_by_round(all_counts)
        print("\tcollapsed by round")
        collapsed = add_target_base_counts(collapsed)
        print("\tadded target base count")
        collapsed = collapsed.sort_values(by=['RunName','LibName','Spike','TargetPattern'],ascending=[True,True,False,True])
        print("\tvalues sorted")

        print("\twritten to csv")

        # Generate the collapsed codes format
        collapsed_codes_path = run+'_collapsed_codes.csv'
        print('Collapsing codes')

        collapsed = collapse_sub_codes(collapsed)
        collapsed = add_target_base_counts(collapsed)
        collapsed.to_csv('{}_collapsed_codes{}{}.csv'.format(run, make_align_label, one_tube_label),index=False)
        print('Done with',run)
    print("time: {}".format(time.time() - t0))



