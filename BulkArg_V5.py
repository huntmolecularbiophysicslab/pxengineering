# -*- coding: utf-8 -*-

#Author: Nooriel Banayan
#Contributors: John Hunt, Henry Hunt, Quangi Lu

import os
import sys
import time
import cupy as cp
import numpy as np
from numpy import inf
from sys import exit
import pandas as pd
from Bio import AlignIO
import matplotlib as mpl
from collections import Counter
from itertools import combinations
import matplotlib.pyplot as plt
from Bio.Align import MultipleSeqAlignment
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.dpi'] = 300


def find_index(string, character):
    return [i for i, ltr in enumerate(string) if ltr == character]

def percent_identity(seq1, seq2):
    
    #length of the sequence
    seq_length1 = len(seq1)
    seq_length2 = len(seq2)
    if seq_length1 != seq_length2:
        raise RuntimeError("Sequences in alignment are not the same length")
    match = 0
    #for each site between two sequences
    for site1, site2 in zip(seq1, seq2):
        #check to see if they are the same
        if site1 == site2:
            match += 1
        else:
            pass
    
    mismatch = seq_length1 - match 
    perc_id = float(match/seq_length1)
    
    return perc_id, mismatch
        

def get_refseq_index(alignment, refseq_id):
    
    #Get the row index of the reference sequence
    for idx, seq in enumerate(alignment):
        if seq.id == refseq_id:
            break
    return idx
    

def refseq_columns(alignment):
    
    #Get the column positions of the sites in the reference sequence that
    #do not contain a dash
    
    #create an empty list
    refseq_idx = []
    
    #iterate through sites in the reference sequence
    for idx in range(len(alignment)):
        
        #if the site is not a dash, thats a hit and append the index
        #to the list
        if alignment[idx] != "-":
            refseq_idx.append(idx)
            
    return refseq_idx

def clean_alignment(aln, refseq_idx):
    n = len(aln[0])
    
    #create a new MultipleSequenceAlignment Object Starting with record ids
    new_alignment = aln[:, :0]
    
    #iterate through column indices of reference sequence with no dashes
    for idx in refseq_idx:
        
        if idx < n:
            #Add columns with no dashes in the reference sequence to the new 
            #alignment object
            new_alignment = new_alignment + aln[:, idx:idx+1]
        
        elif idx == n-1:
            #if the index is the last column of the alignment
            #put the semicolon after the index in the brackets
            #all this to say, Saturday Night Live is NOT funny
            new_alignment = new_alignment + aln[:, idx:]
            
    return new_alignment

def drop_sequences(alignment, cutoff):
    
    #create a new MSA object to store sequences in
    new_alignment_file = MultipleSeqAlignment([])
    seq_length = len(alignment[0])
    dropped = 0
    for i in range(len(alignment)):
        #if the sequence has less dash frequency than the cutoff, keep that sequence
        if alignment[i].seq.count("-")/seq_length <= cutoff:
            #append it to the MSA object
            new_alignment_file.append(alignment[i])
            #count the number of sequences ignored
        else:
            dropped += 1
    print("Dropped {0} sequences with >{1}% dashes".format(dropped, cutoff*100))
    
    return new_alignment_file


def sequence_percent_identities(alignment, ref_index, outpath):
    
    #will be storing everything into a dataframe 
    df = pd.DataFrame(columns = ["SEQID", "PERCENT_IDENTITY", "MISMATCHES", "SEQ_LENGTH"])
    
    #sequence length
    seq_length = len(alignment[0])
    
    #total number of homologs
    number_of_seqs = len(alignment)
    
    #sequence id of all sequences
    seq_id = [alignment[i].id for i in range(number_of_seqs)]
    
    #lists to store percent ID and number of mismatches between homologs and the reference sequnece
    perc_i = []
    mis_m = []
    
    #lopp through all sequences, calculate percent ID between homologs and reference
    #store percent IDs and number of mismatches in lists
    for i in range(number_of_seqs):
        pi, mismatch = percent_identity(alignment[ref_index], alignment[i])
        perc_i.append(pi)
        #add the percent identitty to the dbxrefs attribute of the seqrecord object (i really hate BioPython)
        alignment[i].dbxrefs = pi
        mis_m.append(mismatch)
    
    #add lists as columns to dataframes and save dataframe as csv
    df["SEQID"] = seq_id
    df["PERCENT_IDENTITY"] = perc_i
    df["MISMATCHES"] = mis_m
    df["SEQ_LENGTH"] = seq_length

    df.to_csv(outpath + "/{}_sequence_summaries.csv".format(protein_id))
    
    return df

def bin_percent_identities(alignment):
    
    binned_seqs = {"99" : [], 
                   "90" : [],
                   "80" : [], 
                   "70" : [], 
                   "60" : [], 
                   "50" : [], 
                   "40" : [], 
                   "30" : [], 
                   "20" : [], 
                   "10" : [], 
                   "1" : [], }
    #loop through the alignment file and compare sequence ids with their percent id
    #append sequences to the appropriate bin by its percent id
    for i in alignment:
        if i.dbxrefs >= 0.99:
            binned_seqs["99"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.90 and i.dbxrefs < 0.99:
            binned_seqs["90"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.80 and i.dbxrefs < 0.90:
            binned_seqs["80"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.70 and i.dbxrefs < 0.80:
            binned_seqs["70"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.60 and i.dbxrefs < 0.70:
            binned_seqs["60"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.50 and i.dbxrefs < 0.60:
            binned_seqs["50"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.40 and i.dbxrefs < 0.50:
            binned_seqs["40"].append(list(str(i.seq).upper()))
        elif i.dbxrefs >= 0.30 and i.dbxrefs < 0.40:
            binned_seqs["30"].append(list(str(i.seq).upper()))         
        elif i.dbxrefs >= 0.20 and i.dbxrefs < 0.30:
            binned_seqs["20"].append(list(str(i.seq).upper()))          
        elif i.dbxrefs >= 0.10 and i.dbxrefs < 0.20:
           binned_seqs["10"].append(list(str(i.seq).upper()))
        elif i.dbxrefs < 0.10:
            binned_seqs["1"].append(list(str(i.seq).upper()))
     
    for key in binned_seqs.keys():
        binned_seqs[key] = np.array(binned_seqs[key])
        
    return binned_seqs

def shannon_entropy(AAsequence):

         #Length of sequence
         length = len(AAsequence)
         #Count the occurence of each amino acid in the sequence
         counts = Counter(AAsequence)
         #list to store entropies
         shannon_es = []
         #loop through counter
         for aa, countaa in counts.items():
                 #shannon entropy = sum(-1 *P(x)*log_10(P(x)))
                 freq = countaa/length
                 shannon_es.append( -1*freq*np.log(freq))

         return np.sum(np.array(shannon_es))   

def AA_counting(binned_seqs, refseq_index, old_aa, new_aa, reference_sequence, outpath):
    
    df_dict = {}
    total_homologs_dict = {}
    
    #loop through alignment and percent identities
    #we are only looking at reference sites at this point, not alignment sites
    for pi, seq in binned_seqs.items():
        
        df = pd.DataFrame(columns = ["REFPOS","REFSEQ","MOST_FREQUENT_AA", "SHANNON_ENTROPY",
                                     "A","C","D","E","F","G","H","I","K","L",
                                     "M","N","P","Q","R","S","T","V","W","Y","X","-"])
        
        amino_acids = ["A","C","D","E","F","G","H","I","K","L",
                      "M","N","P","Q","R","S","T","V","W","Y",
                      "X","-"]
        
        df["REFPOS"] = [*range(1,len(binned_seqs["99"][0])+1)]

        #if there is an alignment
        if len(seq) > 0:
            for i in df.index:
                #retrieve total occurences of amino acids at a given reference site
                cnt = Counter(seq[:,i])
                #calculate the shannon entropy at that site
                shannon_ent = shannon_entropy(seq[:,i])
                #store entropy at that site
                df.at[i,"SHANNON_ENTROPY"] = shannon_ent
                #store the reference amino acid at that site
                df.at[i,"REFSEQ"] = reference_sequence[i]
                #store the most frequently occuring amino acid at that site
                df.at[i,"MOST_FREQUENT_AA"] = max(cnt, key = cnt.get)  
                #store the number of times every amino acid occurs at that site
                #store that number
                for aa in list(df.columns[4:]):
                    if aa in cnt.keys():
                        df.at[i, aa] = cnt[aa]
                    elif aa not in cnt.keys():
                        df.at[i, aa] = 0
                        
            df["COUNT_SUM"] = df[amino_acids].sum(axis=1)
            df["FRACTION_{}".format(old_aa)] = df[str(old_aa)].astype("float").divide(df["COUNT_SUM"].astype("float"))
            df["FRACTION_{}".format(new_aa)] = df[str(new_aa)].astype("float").divide(df["COUNT_SUM"].astype("float"))
            df["FRACTION_DASH"] = df["-"].astype("float").divide(df["COUNT_SUM"].astype("float"))
            df["{0}/{1}".format(new_aa, old_aa)] = df["FRACTION_{}".format(new_aa)].astype("float").divide(df["FRACTION_{}".format(old_aa)].astype("float"))     

            #save the csv   
            df.to_csv(outpath + "/{0}_{1}_AAcounting.csv".format(protein_id, pi))
        
            df_dict[pi] = df
            
            total_homologs_dict[pi] = df["COUNT_SUM"][0]
      
        
    fraction_df = pd.DataFrame()
        
    for pi, aa_counts in df_dict.items():
        tmp = aa_counts[aa_counts["REFSEQ"] == str(old_aa)]
        fraction_df["{}%_ID".format(pi)] = tmp["{}".format(new_aa)].astype("float").div(tmp["{}".format(old_aa)].astype("float"))
    
    fraction_df.to_csv(outpath + "/{0}_Ranking.csv".format(protein_id))
        
    return df_dict, fraction_df, total_homologs_dict

def GPU_PairwiseIdentity(seqs):
    
    NUMBER_OF_SEQ = seqs.shape[0]
    LENGTH_OF_SEQ = seqs.shape[1]
    seqsCP = cp.asarray(seqs.view(np.dtype('int32')))
    
    if type(seqs) is not np.ndarray:
        print("pairwiseIdentity was passed a {} requires a 2D numpy array".format(
            type(seqs)))
        exit(1)
    if len(seqs.shape) != 2:
        print("pairwiseIdentity was passed a {}D numpy array requires a 2D numpy array".format(
            len(seqs.shape)))
        exit(1)
    if seqs.dtype != np.dtype('<U1'):
        print("pairwiseIdentity was passed an numpy array with dtype {} requires dtype <U1 ie. a character array.".format(seqs.dtype))
        exit(1)


    identityKernel = cp.ReductionKernel(
    'T x, T y',
    'int32 z',
    'x == y && x != (short)\'-\'',
    'a + b',
    'z = a',
    '0',
    'perIdent'
    )

    dashKernel = cp.ReductionKernel(
    'T x, T y',
    'int32 z',
    'x == (short)\'-\' || y == (short)\'-\'',
    'a + b',
    'z = a',
    '0',
    'perIdent'
    )


    Out_arr = cp.zeros((NUMBER_OF_SEQ, NUMBER_OF_SEQ))
    Out_dash_arr = cp.zeros((NUMBER_OF_SEQ, NUMBER_OF_SEQ))

    for i in range(NUMBER_OF_SEQ):
        Out_arr[:, i] = identityKernel(
            seqsCP[i], seqsCP, axis=1)
        Out_dash_arr[:, i] = dashKernel(seqsCP[i], seqsCP, axis=1)
    
    ADJUSTED_LENGTH_OF_SEQS = LENGTH_OF_SEQ - Out_dash_arr

    return Out_arr.astype('float32')/ADJUSTED_LENGTH_OF_SEQS
               
def pairwise_percent_identities(binned_seqs):
    
    matrices = {}
    

    #loop through bins
    for pi, seqs in binned_seqs.items():
        
        if len(seqs) > 1:
            matrices[pi] = GPU_PairwiseIdentity(seqs)

            fig, ax = plt.subplots(1, figsize = (8,5), sharex = True)
            ax.get_xaxis().set_visible("off")
            ax.get_yaxis().set_visible("off")
            ax.tick_params(axis='both', which='major', labelsize=13)
            ax.set_title("{0}% Bin, {1} Homologs".format(pi, len(matrices[pi].get())) , weight = 'bold')
            im  = ax.imshow(matrices[pi].get(), cmap = "Blues", vmin = 0, vmax = 1)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cb = plt.colorbar(im, cax=cax, orientation='vertical')
            cb.set_label("Pairwise Identity", labelpad=10, y=0.45, rotation = 90, fontsize = 15, weight = 'bold')
            plt.savefig(outpath + "/{0}%-PIM.jpg".format(pi))


        else:
            pass                
    
    return matrices


def linear_redundancy_reduction(matrix, f_min):
    
    A = matrix
    N = A.shape
    A = np.where(A <= f_min, 1, (1-A)/(1-f_min))   
    A = sum(sum(np.tril(A, -1)))        
    A = A * 2/N[0]
    C_non_redundant = 1 + A
 
    return C_non_redundant

def find_lowest_average_index(matrix, indices):
    n = matrix.shape[0]
    
    # Check if we have iterated through all the columns
    if len(indices) >= n:
        return indices
    
    # Calculate the average of the columns based on indices

    average_values = np.mean(matrix[:, indices], axis=1)
    
    # Find the index of the lowest average value

    lowest_index = np.argmin(average_values)
    
    # Append the lowest index to the indices list
    indices.append(lowest_index)
    new_row_values = [100] * n
    for i in indices:
        matrix[i] = new_row_values
    
    return find_lowest_average_index(matrix, indices)

def redundancy_reduction_extreme_insanity(dictionary, f_min, n_remote_homologs):
    sites_to_substitute = [i for i,x in enumerate(reference_sequence) if x == old_aa]
    #remove empty key:value pairs
    dictionary = {i:j for i,j in dist.items() if j.size > 0}
    percent_identities = [str(x) for x in dictionary.keys()]
    ranking_df = pd.DataFrame(columns = percent_identities, index = sites_to_substitute)
    linear_redundancy_reduction_df = pd.DataFrame(columns = percent_identities, index = sites_to_substitute)
    total_df = pd.DataFrame(columns = percent_identities, index = sites_to_substitute)
    
    for pi, mats in dictionary.items():
        for site in sites_to_substitute:
            matrix = []
            for seq in mats:
                if seq[site] == new_aa:
                    matrix.append(seq)
                else:
                    continue
            if len(matrix) > 0:    
                matrix = np.vstack(matrix)
            else:
                continue

            N = matrix.shape
            total_df.loc[site, str(pi)] = N[0]

            if N[0] > 1:
                A = GPU_PairwiseIdentity(matrix) 
                A = A.get()
                #np.savetxt("{0}-{1}.csv".format(pi,site),  A, delimiter=",")
                linear_redundancy_reduction_df.loc[site, str(pi)] = linear_redundancy_reduction(A, 0.3)
                lowest_pair = np.where(A == A.min())[0]
                indices = [lowest_pair[0], lowest_pair[1]]  # Initial list with the first column index
                A_copy = np.copy(A)

                while len(indices) < A.shape[0]:
                    indices = find_lowest_average_index(A_copy, indices)

                pairwise_combinations = list(combinations(indices[0:n_remote_homologs], 2))
                A_resampled = np.where(A <= f_min, 1, (1-A)/(1-f_min))
                A_not_resampled = 1 - A_resampled
                resampled = []
                notresampled = []

                for i, j in pairwise_combinations:
                    resampled.append(np.log2(A_resampled)[i,j])
                    notresampled.append(np.log2(A_not_resampled)[i,j])

                resampled = np.array(resampled)
                notresampled = np.array(notresampled)

                notresampled =  np.where(np.isinf(notresampled), -100, notresampled)
                resampled =  np.where(np.isinf(resampled), -100, resampled)

                counteOfresampled = np.array([0, 1]) 
                sumLog2probabilities = np.array([notresampled[0], resampled[0]])

                for g in range(2, len(resampled)):
                    counteOfresampled_plusone = counteOfresampled + 1
                    counteOfresampled = np.append(counteOfresampled, counteOfresampled_plusone)
                    sumLog2probabilities_pluslog2presampling = sumLog2probabilities + resampled[g-1]
                    sumLog2probabilities = sumLog2probabilities + notresampled[g-1]
                    sumLog2probabilities = np.append(sumLog2probabilities, sumLog2probabilities_pluslog2presampling)

                product_probabilities = 2**sumLog2probabilities
                count_weighted_product_probabilities = counteOfresampled * product_probabilities
                sum_terms = np.sum(count_weighted_product_probabilities)

                ranking_df.loc[site, str(pi)] = 1 + sum_terms
            
            elif  N[0] == 1:
                ranking_df.loc[site, str(pi)] = 1
                linear_redundancy_reduction_df.loc[site, str(pi)] = 1
            
            else:
                ranking_df.loc[site, str(pi)] = 0
                linear_redundancy_reduction_df.loc[site, str(pi)] = 0 

    ranking_df = ranking_df.fillna(0)
    total_df = total_df.fillna(0) 
    linear_redundancy_reduction_df = linear_redundancy_reduction_df.fillna(0)
    linear_redundancy_reduction_df["sites"] = linear_redundancy_reduction_df.index + 1 
    ranking_df.to_csv(outpath + "/Expectation_probabilities.csv")
    total_df.to_csv(outpath + "/total_homologs.csv")
    linear_redundancy_reduction_df.to_csv(outpath + "/linear_redundancy_correction.csv")

    return ranking_df, total_df, linear_redundancy_reduction_df 

def preferred_ranking(ranking_df):
    
    order = []
    
    for col in ranking_df.columns:
        if col != "99":
            icol = ranking_df[col].sort_values(ascending = False)
            #icol = icol[icol != 0]
            icol = icol[icol > 1]
            for idx in icol.index:
                if idx in order:
                    pass
                else:
                    order.append(idx)
    
    return [int(x) + 1 for x in order]

def preferred_ranking_2(ranking_df, total_df, rrc, threshold):
    
    order = {}
    
    for col in ranking_df.columns:
        if col != "99" and col != "sites":
            icol = ranking_df[col].sort_values(ascending = False)
            icol = icol[icol >= threshold]
            for idx in icol.index:
                if idx + 1 in order.keys():
                    continue
                else:
                    idx_actual = int(idx) + 1
                    order[idx_actual] = (col, total_df.loc[idx, col], rrc.loc[idx, col], ranking_df.loc[idx, col] )
                    
    
    return order

def valhalla_plots(dictionary_all, rrc, total_df, expectation):
    
    dictionary = dictionary_all[0]
    x = len(dictionary)
    key = list(dictionary.keys())[0]
    df = dictionary[key]
    y = len(df[df["REFSEQ"] == "{}".format(old_aa)])
    keys = [int(x) for x in list(dictionary.keys())]
    keys = sorted(keys, reverse=True)
    array = np.zeros(shape = (y,x))
    array2 = np.zeros(shape = (y,x))
    array3 = np.zeros(shape = (y,x))
                  

    rrc = rrc.sort_index(ascending=False)
    dims = rrc.T.shape
    x = dims[0] 
    y = dims[1] 
    ar = np.zeros(shape = (x,y))
    for i, col in enumerate(list(rrc.columns)):
 
        ar[i , :] = rrc[col]

        
    for i,k in enumerate(keys):
        temp = dictionary[str(k)][dictionary[str(k)]["REFSEQ"] == "{}".format(old_aa)]
        array[:,i] = temp["SHANNON_ENTROPY"]
    for i,k in enumerate(keys):
        temp = dictionary[str(k)][dictionary[str(k)]["REFSEQ"] == "{}".format(old_aa)]
        array2[:,i] = 1 - temp["FRACTION_{}".format(old_aa)]
    for i,k in enumerate(keys):
        temp = dictionary[str(k)][dictionary[str(k)]["REFSEQ"] == "{}".format(old_aa)]
        t = list(temp["FRACTION_{}".format(new_aa)] / temp["FRACTION_{}".format(old_aa)])
        for k,h in enumerate(t):
            if h > 2:
                t[k] = 2
            else:
                pass
                
      
        array3[:,i] = t

    keys = [str(x) +"%" for x in keys]
    #fig, (ax1,ax2,ax3, ax4) = plt.subplots(4, figsize = (10,10), sharex = True)
    fig, (ax1,ax2,ax3) = plt.subplots(3, figsize = (15,10))
    #fig.supylabel("Percent ID Bins (Exclusive)", x = 0)
    im  = ax1.imshow(array.T, cmap = "Blues")
    im2  = ax2.imshow(array2.T, cmap = "Blues")
    im3  = ax3.imshow(array3.T, cmap = "Blues", vmin = 0, vmax = 2)
    #im4  = ax4.imshow(ar, cmap = "Blues")


    divider = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    divider3 = make_axes_locatable(ax3)
    #divider4 = make_axes_locatable(ax4)

    cax = divider.append_axes('right', size='5%', pad=0.05)

    ax1.set_yticks(np.arange(len(keys)))
    ax1.set_yticklabels(keys)
    ax1.set_xticks(np.arange(len(temp["REFPOS"])))
    ax1.set_xticklabels(list(temp["REFPOS"]))
    ax1.set_xticklabels(list(temp["REFPOS"]), rotation = 43, size = 10)
    for key, spine in ax1.spines.items():
        spine.set_visible(False)
        cb = plt.colorbar(im, cax=cax, orientation='vertical')

    cb.set_label("Shannon Entropy \n $-\sum p(AA) ln(p(AA))$", labelpad=50, y=0.45, fontsize = 18, rotation = 270, weight = 'bold')


    ax2.set_yticks(np.arange(len(keys)))
    ax2.set_yticklabels(keys)
    ax2.set_xticks(np.arange(len(temp["REFPOS"])))
    ax2.set_xticklabels(list(temp["REFPOS"]))
    ax2.set_xticklabels(list(temp["REFPOS"]), rotation = 43, size = 10)
    for key, spine in ax2.spines.items():
        spine.set_visible(False)
    cax = divider2.append_axes('right', size='5%', pad=0.05)
        
    cb = fig.colorbar(im2, cax=cax, orientation='vertical')
    cb.set_label("Probability of \n Not Being {}".format(old_aa), labelpad=50, y=0.45, fontsize = 18, rotation = 270, weight = 'bold')

    ax3.set_yticks(np.arange(len(keys)))
    ax3.set_yticklabels(keys)
    ax3.set_xticks(np.arange(len(temp["REFPOS"])))
    ax3.set_xticklabels(list(temp["REFPOS"]), rotation = 43, size = 10)
    for key, spine in ax3.spines.items():
        spine.set_visible(False)
    cax = divider3.append_axes('right', size='5%', pad=0.05)

    cb = fig.colorbar(im3, cax=cax, orientation='vertical', ticks=[0,0.5, 1, 1.5, 2])
    cb.ax.set_yticklabels(['0', '0.5', '1', '1.5', '>= 2'])
    cb.set_label("p({0})/p({1})".format(new_aa,old_aa), labelpad=42, y=0.45, fontsize = 18, rotation = 270, weight = 'bold')

 
    #ax4.set_yticks(np.arange(len(keys)))
    #ax4.set_yticklabels(keys)
    #ax4.set_xticks(np.arange(len(temp["REFPOS"])), rotation = 45)
    #ax4.set_xticklabels(list(temp["REFPOS"]), rotation = 45)
    #for key, spine in ax4.spines.items():
    #    spine.set_visible(False)
    #cax = divider4.append_axes('right', size='5%', pad=0.05)

    #cb = fig.colorbar(im4, cax=cax, orientation='vertical')
    #cb.set_label("Redundancy \n Reduced \n {} \n Counts".format(new_aa), labelpad=50, y=0.45, rotation = 360)


    #ax4.set_xlabel("Reference {} Sites".format(old_aa))
   # ax3.set_xlabel("{} Sites in Target Protein".format(old_aa), size = 16, weight = 'bold')

    ax3.set_xlabel("({}) positions \n in target protein".format(old_aa), size = 16, weight = 'bold')


    plt.savefig(outpath + "/maps.jpg")
    
    rrc["sites"] = rrc.index.astype("int")  + 1
    rrc = rrc.sort_values(by=['sites'])
    rrc = rrc.loc[:, (rrc != 0).any(axis=0)]
    x = range(len(rrc))
    fig, (ax2,ax1, ax3) = plt.subplots(nrows =3, ncols = 1, figsize=(15,20))
    fig.subplots_adjust(hspace=.7)
    ax1.set_yscale('log', base=2)
    for i in rrc.columns:
        if i != "sites":
            y = rrc[i]
            y[y == 0] = 0.125
            ax1.plot(x, y, "-o", markersize = 10, label = i,  linewidth = '5')
    ax1.set_ylim([0.5, 10000])
    ax1.set_xticks(np.arange(len(rrc["sites"])))
    ax1.set_xticklabels(list(rrc["sites"]), fontsize = 20)
    ax1.tick_params(axis='y', labelsize=20)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize = 18)
    ax1.set_ylabel("Redundancy-Reduced {0} Count".format(new_aa, old_aa), fontsize = 18, weight = 'bold')

    dictionary_2 = total_df
    dictionary_2["sites"] = dictionary_2.index.astype("int")  + 1
    dictionary_2 = dictionary_2.sort_values(by=['sites'])
    dictionary_2 = dictionary_2.loc[:, (dictionary_2 != 0).any(axis=0)]
    x = range(len(dictionary_2))
    ax2.set_yscale('log', base=2)
    for i in dictionary_2.columns:
        if i != "sites":
            y = dictionary_2[i]
            y[y == 0] = 0.125
            ax2.plot(x, y, "-o", markersize = 10, label = i, linewidth = '5')
    ax2.set_ylim([0.5, 10000])
    ax2.set_xticks(np.arange(len(dictionary_2["sites"])))
    ax2.set_xticklabels(list(dictionary_2["sites"]), fontsize = 20)
    ax2.tick_params(axis='y', labelsize=20)
    plt.legend(bbox_to_anchor=(0.5, 3.0), loc='upper center', ncol = 10, borderaxespad=0, fontsize=18)
    ax2.set_ylabel("Total {0} Count".format(new_aa, old_aa), fontsize = 18, weight = 'bold')

    dictionary_3 = expectation
    dictionary_3["sites"] = dictionary_3.index.astype("int")  + 1
    dictionary_3 = dictionary_3.sort_values(by=['sites'])
    dictionary_3 = dictionary_3.loc[:, (dictionary_3 != 0).any(axis=0)]
    x = range(len(dictionary_3))
    #ax2.set_yscale('log', base=2)
    for i in dictionary_3.columns:
        if i != "sites":
            y = dictionary_3[i]
            y[y == 0] = 0.125
            ax3.plot(x, y, "-o", markersize = 10,label = i, linewidth = '5')
    ax3.set_xticks(np.arange(len(dictionary_3["sites"])))
    ax3.set_xticklabels(list(dictionary_3["sites"]), fontsize = 20)
    ax3.tick_params(axis='y', labelsize=20)
    ax3.set_yscale('log', base=2)
    #plt.legend(bbox_to_anchor=(0.5, 5.0), title = "Minimum % ID in Exclusive Bins", fontsize = 18, loc='upper center', ncol = 10, borderaxespad=0)
    legend = ax3.legend(bbox_to_anchor=(0.5, 5.0), title = "Minimum % ID in Exclusive Bins", fontsize = 20, loc='upper center', ncol = 10, borderaxespad=0)
    plt.setp(legend.get_title(),fontsize=25)

    ax3.set_ylabel("Independent Observations \n in 7 Most Remote Homologs".format(new_aa, old_aa), fontsize = 18, weight = 'bold')





    #plt.xlabel("{} Sites in Target Protein".format(old_aa), fontsize = 20, weight = 'bold')
    plt.xlabel("({}) positions \n in target protein".format(old_aa), fontsize = 30, weight = 'bold')
    plt.savefig(outpath + "/lines.jpg")

    return 

def bonding_partners(refseq, old_aa):

    
    site_indeces = np.array(find_index(refseq, old_aa))
    df = pd.DataFrame({"{}-SITES".format(old_aa): site_indeces + 1, 
                        "i-4": "","i-3": "","i-2": "","i-1": "","i": "{}".format(old_aa.lower()),"i+1": "",
                        "i+2": "","i+3": "","i+4": ""})
    
    E_bonding_partners = ["K", "R", "H", "N", "Q"]
    Q_bonding_partners = ["K", "R", "H", "N", "Q", "D", "E"]

    for j,i in enumerate(site_indeces):
        if i-4 >= 0:
            if old_aa == "D":
                if refseq[i - 4] in E_bonding_partners:
                    df.loc[j,"i-4"] = refseq[i - 4]
                else:
                    df.loc[j,"i-4"] = refseq[i - 4].lower()
            elif old_aa == "N":
                if refseq[i - 4] in Q_bonding_partners:
                    df.loc[j,"i-4"] = refseq[i - 4]
                else:
                    df.loc[j,"i-4"] = refseq[i - 4].lower()
        if i-3 >= 0:
            if old_aa == "D":
                if refseq[i - 3] in E_bonding_partners:
                    df.loc[j,"i-3"] = refseq[i - 3]
                else:
                    df.loc[j,"i-3"] = refseq[i - 3].lower()
            elif old_aa == "N":
                if refseq[i - 3] in Q_bonding_partners:
                    df.loc[j,"i-3"] = refseq[i - 3]
                else:
                    df.loc[j,"i-3"] = refseq[i - 3].lower()
        if i-2 >= 0:
            if old_aa == "D":
                if refseq[i - 2] in E_bonding_partners:
                    df.loc[j,"i-2"] = refseq[i - 2]
                else:
                    df.loc[j,"i-2"] = refseq[i - 2].lower()
            elif old_aa == "N":
                if refseq[i - 2] in Q_bonding_partners:
                    df.loc[j,"i-2"] = refseq[i - 2]
                else:
                    df.loc[j,"i-2"] = refseq[i - 2].lower()
        if i-1 >= 0:
            df.loc[j,"i-1"] = refseq[i - 1].lower()
        if i+1 < len(refseq):
            df.loc[j,"i+1"] = refseq[i + 1].lower()
        if i+2 < len(refseq):
            if old_aa == "D":
                if refseq[i + 2] in E_bonding_partners:
                    df.loc[j,"i+2"] = refseq[i + 2]
                else:
                    df.loc[j,"i+2"] = refseq[i + 2].lower()
            elif old_aa == "N":
                if refseq[i + 2] in Q_bonding_partners:
                    df.loc[j,"i+2"] = refseq[i + 2]
                else:
                    df.loc[j,"i+2"] = refseq[i + 2].lower()    
        if i+3 < len(refseq):    
            if old_aa == "D":
                if refseq[i + 3] in E_bonding_partners:
                    df.loc[j,"i+3"] = refseq[i + 3]
                else:
                    df.loc[j,"i+3"] = refseq[i + 3].lower()
            elif old_aa == "N":
                if refseq[i + 3] in Q_bonding_partners:
                    df.loc[j,"i+3"] = refseq[i + 3]
                else:
                    df.loc[j,"i+3"] = refseq[i + 3].lower()
        if i+4 < len(refseq):    
            if old_aa == "D":
                if refseq[i + 4] in E_bonding_partners:
                    df.loc[j,"i+4"] = refseq[i + 4]
                else:
                    df.loc[j,"i+4"] = refseq[i + 4].lower()
            elif old_aa == "N":
                if refseq[i + 4] in Q_bonding_partners:
                    df.loc[j,"i+4"] = refseq[i + 4]
                else:
                    df.loc[j,"i+4"] = refseq[i + 4].lower()

        df = df.replace(r'^\s*$', np.nan, regex=True)
    
    
    df.to_csv(outpath + "/BondingPartners.csv")
    return df
    


if __name__ == "__main__":
    t0 = time.time()
    file_path = "/mnt/c/Users/Nooriel/Desktop/"
    file_name = "pdi.aln"
    alignment = AlignIO.read(file_path + file_name, "clustal")
    protein_id = "pdi"
    old_aa = "K"
    new_aa = "R"
    outpath = file_path + "/" + protein_id
    if os.path.exists(outpath) == False:
        os.mkdir(outpath)
    else:
        pass
    reference_index = get_refseq_index(alignment, protein_id)
    print("Loding your alignment file: {}".format(alignment))
    column_idx = refseq_columns(alignment[reference_index])
    print("Processing")
    alignment_nogaps = clean_alignment(alignment, column_idx)
    reference_sequence = alignment_nogaps[reference_index].seq
    print("Cleaning")
    alignment_final = drop_sequences(alignment_nogaps, 0.25)
    print("Churning")
    sequence_summary = sequence_percent_identities(alignment_final, reference_index, outpath)
    print("Spooling")
    dist = bin_percent_identities(alignment_final)
    print("Mending")
    aa_counts_by_pi = AA_counting(dist, reference_index, old_aa, new_aa, reference_sequence, outpath)
    print("Yarning")
    matrices = pairwise_percent_identities(dist)
    print("Singing")
    f_min = 0.3
    extreme_expectation, total_df, linear_rrc = redundancy_reduction_extreme_insanity(aa_counts_by_pi[0], f_min, 7)
    #ar = valhalla_plots(aa_counts_by_pi[0], resampling)
    ar = valhalla_plots(aa_counts_by_pi, linear_rrc, total_df, extreme_expectation)    
    #sites_to_sub = preferred_ranking(extreme_expectation)
    sites_to_sub = preferred_ranking_2(extreme_expectation,  total_df, linear_rrc, threshold = 1.1)
    with open(outpath + '/preferred_sites_to_substitute.txt', 'w') as f:
    # loop through the list and write each item to the file
        f.write("Recommended order of sites for introduction of {0}-to-{1} mutations.\n".format(old_aa,new_aa))
        f.write("Position" + "      " + "%ID Bin" + "       " + "Total" + "     "+ "Redundancy Reduced" + "     "+ "7 most remote homologs" + "\n")
        for site, counts in sites_to_sub.items():
            f.write(str(site) + "   " + str(counts[0]) + "   " + str(counts[1]) + "   " + str(counts[2]) + "   " + str(counts[3]) + "\n")
    if (old_aa == "D") and (new_aa == "E"):
        bonding_partners(reference_sequence, old_aa)
    elif (old_aa == "N") and (new_aa == "Q"):
        bonding_partners(reference_sequence, old_aa)
    print("Here you go!")
    t1 = time.time()
    total = t1- t0

sys.exit()