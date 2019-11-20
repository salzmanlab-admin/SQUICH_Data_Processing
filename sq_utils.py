#!/bin/bash
#############################
# File Name : squish_utils.py
#
# Purpose : Functions manipulating sequencing reads
#
# Creation Date : 10-05-2018
#
# Last Modified : Wed 23 May 2018 12:03:17 PM PDT
#
# Created By : Rob Bierman
#
##############################

from Bio import SeqIO

def motif_count(fastq_name,motif):
    total_reads = 0
    motif_reads = 0
    for seq in read_fastq(fastq_name):
        total_reads += 1
        if contains_motif(seq,motif):
            motif_reads += 1

    return motif_reads,total_reads


def read_fastq(fastq_name,seqs_only=True):
    print("reading fastq")
    in_loop = False
    for r in SeqIO.parse(fastq_name, 'fastq'):
        if not in_loop:
            in_loop = True
            print("in loop {}".format(r)) 
        
        if seqs_only:
            yield str(r.seq)
        else:
            yield r


def contains_motif(seq,motif,max_mismatches=0):
    matches = {'A':['A'],'T':['T'],'C':['C'],'G':['G'],
               'S':['C','G'],'W':['A','T']}

    for i in range(len(seq)-len(motif)+1):
        sub_seq = seq[i:i+len(motif)]
        found = True
        mismatches = 0
        for i,b in enumerate(sub_seq):
            mismatches += 0 if b in matches[motif[i]] else 1
            if mismatches > max_mismatches:
                found = False
                break

        if found:
            return True
            
    return False

def degen_seqs():
    degens = {'A':['A'],'C':['C'],'G':['G'],'T':['T'],
              'R':['A','G'],'Y':['C','T'],'S':['C','G'],'W':['A','T'],'K':['G','T'],'M':['A','C'],
              'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],
              'N':['A','C','G','T']}

    return degens

def expand_degens(seq):
    degens = degen_seqs()

    def expand_degens_recur(seq,builds=['']):
        if not seq:
            return builds

        h,t = seq[0],seq[1:]
        builds = [build+degen for degen in degens[h] for build in builds]
        return expand_degens_recur(t,builds)

    return expand_degens_recur(seq)


def score_v1(collapsed):
    collapsed['ScoreV1'] = 0
    for i,r in enumerate(rounds[:-1]):
        weighted = collapsed[r]*(10**i)
        collapsed['ScoreV1'] += np.where(i == 0 or collapsed['ScoreV1']+10 > weighted, weighted, 0)

    return collapsed

def rev_comp(seq):
    return comp_seq(seq)[::-1]

def comp_seq(seq):
    comps = {'A':'T','T':'A','G':'C','C':'G','S':'S', 'N':'N'}
    return ''.join(comps[b] for b in seq)

