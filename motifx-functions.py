#!/usr/bin/python

#######################################
# motif-x implementation in python    #
# Created on : 05-10-18               #
# Author     : Ziqi DENG              #
# Contact    : dengziqi1234@gmail.com #
#######################################

###########
# IMPORTS #
###########

import pandas as pd
#import numpy as np
import math
import scipy.stats
import sys
import re

#############
# FUNCTIONS #
#############

def GetSeqs(filename):
    """
    :param filename: file of pre-align sequences
    :return: list of sequences
    """

    seqs_set = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()
            seqs_set.append(line)
    return seqs_set


def GetAA():
    """
    :return: All amino acids array
    """
    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return AA


def BuildPWM(seqs, relative_freq=True):
    """
    :param seqs: list of sequences
    :param relative_freq: bool
    :return: A Dataframe of Position Weight Matrix
    """

    # Ensure same length characters
    seq_len = len(seqs[0])
    num_pos = seq_len
    for seq in seqs:
        seq_len_test = len(seq)
        if seq_len_test != num_pos:
            print('Unequal length of sequences')
            break

    # List of valid amino acids, sorted
    namespace = GetAA()

    # Make matrix of letters
    seq_matrix = [list(seq) for seq in seqs]
    df = pd.DataFrame(seq_matrix)

    # Get frequencies of aa
    pieces = []
    for column in df.columns:
        pos_series = df[column].value_counts()
        pos_series.name = column
        pieces.append(pos_series)
    pwm_matrix_df = pd.concat(pieces, axis=1,sort=True)

    # Reindex as Amino acids to remove non-standard letters
    pwm_matrix_df = pwm_matrix_df.reindex(namespace)
    #pwm_matrix_df = pwm_matrix_df.loc[namespace, :]

    pwm_matrix_df.fillna(0, inplace=True)

    # Do relative frequencies
    if relative_freq:
        pwm_matrix_df = pwm_matrix_df.loc[:, :] / pwm_matrix_df.sum()

    # Name to columns as positions
    pwm_matrix_df.columns = [i + 1 for i in range(num_pos)]

    return pwm_matrix_df


def MotifRegex(motifs, cent_regex, n=11):
    """
    :param motifs: significant residue and its position at the flank
    :param cent_regex: Central amino acid
    :param n: Length of peptide
    :return: regular expression of consensus motifs, e.g '..S..P.....'
    """

    # Central residue position
    cent_ind = math.ceil(n / 2)

    # Must have at least one row
    if motifs.shape[0] > 0:
        z = list('.' * n)
        for index, row in motifs.iterrows():
            aa = str(row['AA'])
            pos = row['pos'] - 1
            z[pos] = aa
            z[cent_ind - 1] = cent_regex
            z_str = ''.join(z)
        regex = ''.join(z)

    return regex


def FindMotif(fg_seqs, bg_seqs, min_seqs, pval_cutoff, cent_regex, verbose=True, perl_impl=True):

    AA = GetAA()

    fg_size = len(fg_seqs)
    bg_size = len(bg_seqs)

    kmer_length = len(fg_seqs[0])
    cent_index = math.ceil(kmer_length / 2)

    if verbose:
        print('Step 1: recursive motif building')

    motif = []
    pvals = []

    iter = 1
    while True:
        # Break if we dont have any sequences remaining
        if len(bg_seqs) == 0 or len(fg_seqs) == 0:
            break


        # Construct PFM and PWM
        ppwm = BuildPWM(fg_seqs, relative_freq=False)
        npwm = BuildPWM(bg_seqs)

        # Compute binomial matrix
        bin = 1 - scipy.stats.binom.cdf(ppwm - 1, len(fg_seqs), npwm)

        binomial_df = pd.DataFrame(bin, columns=[i + 1 for i in range(kmer_length)], index=AA)

        # Lowest p-value computed by pbinom is 10e-16, set any zero values to that
        if perl_impl is True:
            if verbose:
                print('Using perl implementation.')
            binomial_df[binomial_df == 0] = 1e-16  # this is what weblogo uses in perl, different for python
        else:

            binomial_df[binomial_df == 0] = sys.float_info.min  # Usually ~= 2.2e-308

        # Don't enrich on central residue or previous motifs
        binomial_df.loc[:, cent_index] = 1

        if len(motif) > 0:

            for mo in motif:
                binomial_df.loc[mo[0], mo[1] - 1] = 1


        # Give anything with a frequency less than the threshold a p-value of 1
        binomial_df[ppwm < min_seqs] = 1

        # Find the row/column of the min binomal p-value
        min_bin = binomial_df.min().min()


        mbin = []
        for index in AA:
            for column in range(kmer_length):
                cell_value = binomial_df.loc[index,column+1]
                if ( cell_value == min_bin) & (cell_value < pval_cutoff):
                    mbin.append([index,column+1]) #1-base index

        # No more significant, breaking
        if len(mbin) == 0:
            break

        # Case where there are more than one matches (likely pvalue=0)
        # Find match with greatest value in PFM
        if len(mbin) > 1:
            max_freq = 0
            for item in mbin:
                freq = ppwm.loc[item[0], item[1]]
                if freq > max_freq:
                    max_freq = freq
                    mbin_new = [item[0], item[1]]
                    mbin = [mbin_new]

        # format mbin = [('S',5)]

        aa = mbin[0][0]
        index = mbin[0][0]
        col = mbin[0][1]

        if verbose:
            print('\t Iteration {}, aa={} row={} col={}, pval={}, fg={}, bg={}'.format(iter, aa, index, col, min_bin,
                                                                                       len(fg_seqs), len(bg_seqs)))

        # Extract sequences
        fg_seqs = [seq for seq in fg_seqs if seq[col - 1] == aa]
        bg_seqs = [seq for seq in bg_seqs if seq[col - 1] == aa]
        motif.append([index, col])


        pvals.append(min_bin)

        iter += 1

    # Motif data: data frame with amino acid and positions for the motif
    motif_data = pd.DataFrame(motif, columns=['AA', 'pos'])

    # If we don't have any data, return NULL so that the parent function can deal with it
    if len(motif_data['AA']) == 0:
        return None

    # Find the regex of the motif given the row/col and aa values, e.g "....SP....."
    motif_regex = MotifRegex(motif_data, cent_regex, kmer_length)
    #print(motif_regex)

    # Compute the score of the motif
    pvals_score = []
    if len(pvals) == 0:
        motif_score = None
    else:
        pvals_score = [-math.log(pval, 10) for pval in pvals]
        motif_score = sum(pvals_score)

    # Compute fg/bg matches
    fg_matches = [line for line in fg_seqs if re.match(motif_regex, line)]
    bg_matches = [line for line in bg_seqs if re.match(motif_regex, line)]

    fg_matches_count = len(fg_matches)
    bg_matches_count = len(bg_matches)

    return fg_seqs, motif_data, motif_score, motif_regex, fg_matches_count, fg_size, bg_matches_count, bg_size


###########
#  motifx #
###########

"""
# ' Find overrepresented sequence motifs
# ' 
# ' Parameters:
# '
# ' fg.seqs Foreground k-mer sequences in a pre-aligned format. All k-mers must have same lengths.
# ' bg.seqs Background k-mer sequences in a pre-aligned format. All k-mers must have same lengths.
# ' central.res Central amino acid of the k-mer. Sequences without this amino acid in the centre position are filtered out. This can be one or more letter. For example, 'S', 'ST', 'Y', or 'STY'.
# ' min.seqs This threshold refers to the minimum number of times you wish each of your extracted motifs to occur in the data set. An occurrence threshold of 20 usually is appropriate, although this parameter may be adjusted to yield more specific or less specific motifs.
# ' pval.cutoff The p-value threshold for the binomial probability. This is used for the selection of significant residue/position pairs in the motif. A threshold of 0.000001 is suggested to maintain a low false positive rate in standard protein motif analyses.
# ' verbose If true, motifx will show textual details of the steps while running.
# ' perl.impl The original implementation of motifx in perl, P-values below 1e-16 cannot be computed and are thus set to zero. Motifx therefore sets any P-values of zero to the minimal P-value of 1e-16. In R, the minimal P-value is much lower (depending on the machine). If this option is set to TRUE, P-values with a value of zero are set to 1e-16, as in perl. Otherwise, the R P-value minimum will be used. For results identical to that of the webserver implementation, set to TRUE.
# ' 
# ' @return Data frame with seven columns containing overrepresented motifs. Motifs are listed in the order in which they are extracted by the algorithm, not with regard to statistical significance. Thus it should not be assumed that a motif found at a higher position in the list is more statistically significant than a motif found at a lower position in the list. The columns are as follows:
# ' \describe{
# '   \item{motif}{The overrepresented motif}
# '   \item{score}{The motif score, which is calculated by taking the sum of the negative log probabilities used to fix each position of the motif. Higher motif scores typically correspond to motifs that are more statistically significant as well as more specific }
# '   \item{fg.matches}{Frequency of sequences matching this motif in the foreground set}
# '   \item{fg.size}{Total number of foreground sequences}
# '   \item{bg.matches}{Frequency of sequences matching this motif in the background set}
# '   \item{bg.size}{Total number of background sequences}
# '   \item{fold.increase}{An indicator of the enrichment level of the extracted motifs. Specifically, it is calculated as (foreground matches/foreground size)/(background matches/background size).}
# ' }
"""

def motifx(fg_seqs, bg_seqs, central_res, min_seqs=20, pval_cutoff=1e-6, verbose=True, perl_impl=True):
    AA = GetAA()

    def CheckEmptySeqs():
        if len(fg_seqs) == 0:
            print('Could not find any foreground sequences!')
            return
        if len(bg_seqs) == 0:
            print('Could not find any background sequences!')
            return
        return

    CheckEmptySeqs()

    if len(central_res) > 1:
        cent_regex = '[' + central_res + ']'
    else:
        cent_regex = central_res

    c_res = central_res
    c_res = list(set(c_res).intersection(set(AA)))

    if len(c_res) == 0:
        print('Central residue must contain at least one amino acid character')
        return

    # Check sequence widths
    width = len(fg_seqs[0])
    if width < 3 or width > 35:
        print('Sequence width must be between 3 and 35!')
        return
    if width % 2 == 0:
        print('Sequence width must be an odd number!')
        return
    if width != len(bg_seqs[0]):
        print('Widths for foreground and background data must be equal!')
        return

    # Ensure k-mers have same lengths
    nc_pos = [len(seq) for seq in fg_seqs]
    if any(i != nc_pos[0] for i in nc_pos):
        print('Foreground k-mers must be same lentth.')
        return
    nc_bg = [len(seq) for seq in bg_seqs]
    if any(i != nc_pos[0] for i in nc_bg):
        print('Background k-mers must be same lentth.')
        return

    # Get central index
    ci = math.ceil(width / 2)

    # Filter the central residue for allowed residues
    fg_seqs = [seq for seq in fg_seqs if seq[ci - 1] in c_res]
    bg_seqs = [seq for seq in bg_seqs if seq[ci - 1] in c_res]

    if True:
        # Only set to true for exact match to motifx webserver
        # Remove non-amino acid residues
        fg_seqs = [seq for seq in fg_seqs if re.match('\-|\*|[^BJOUXZ]', seq)]
        bg_seqs = [seq for seq in bg_seqs if re.match('\-|\*|[^BJOUXZ]', seq)]
        bg_seqs = list(set(bg_seqs))

    CheckEmptySeqs()
    seqs = fg_seqs

    data = []
    while True:
        # Find the motif
        mt = FindMotif(fg_seqs=seqs, bg_seqs=bg_seqs, min_seqs=min_seqs,
                       pval_cutoff=pval_cutoff, cent_regex=cent_regex, verbose=verbose,
                       perl_impl=perl_impl)

        # mt output list(fg_seqs, motif_data, motif_score,motif_regex,fg_matches_count,fg_size,bg_matches_count, bg_size)
        if not mt:
            #print('No more motif to found')
            break

        # Append to list of data
        data.append(list(mt))

        # Remove stuff already apart of a motif
        if verbose:
            print('Step 2: positive and negative set reduction')
        motif_regex = mt[3]

        seqs = [seq for seq in seqs if not re.match(motif_regex, seq)]

        bg_seqs = [seq for seq in bg_seqs if not re.match(motif_regex, seq)]

        # No more sequences left to process, break
        if len(seqs) < min_seqs:
            print('No more sequences left to process')
            break

    if verbose:
        print('Converged, no more enrichments!')

    df_result = pd.DataFrame(data, columns=['fg_seqs', 'motif_data', 'motif_score', 'motif_regex',
                                            'fg_matches', 'fg_size', 'bg_matches', 'bg_size'])

    del df_result['fg_seqs']
    del df_result['motif_data']

    df_result['fold_increase'] = (df_result['fg_matches'] / df_result['fg_size']) / (
                df_result['bg_matches'] / df_result['bg_size'])

    return df_result

########
# MAIN #
########

if __name__ == '__main__':
    fg_seqs = GetSeqs(sys.argv[1])
    bg_seqs = GetSeqs(sys.argv[2])
    central_res = sys.argv[3]
    min_seqs = int(sys.argv[4])
    pval_cutoff = float(sys.argv[5])
    verbose = sys.argv[6]
    perl_impl = sys.argv[7]
    mot = motifx(fg_seqs, bg_seqs, central_res, min_seqs, pval_cutoff, verbose, perl_impl)
    print(mot)
