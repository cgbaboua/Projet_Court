import os
import sys 
import numpy as np

# Read all the sequences from a FASTA file and store them in a dictionary
def read_sequences(fasta_file) :
    # Dictionary to hold the sequences
    sequences = {}  
    with open(fasta_file,'r') as fasta :
        for line in fasta :
            # Identifier first line 
            if line.startswith('>') : 
                query_seq = []
                # Check the FASTA format (uniprot or ncbi)
                end_accession = line[4:].find("|")
                if end_accession == -1 :
                    end_accession = line.find(' ')
                    seq_accession = line[1:end_accession]
                else :
                    seq_accession = line[4:end_accession+4]
            else :
                # Sequence line
                query_seq.extend(list(line.replace(' ','').strip()))
                sequences[seq_accession]=query_seq
    return sequences

# Create the BLOSUM62 matrix stored in a dictionary
def matrix(matrix_file):
    # Dictionary to hold the BLOSUM62 values
	blosum62 = {}
    # List to hold the amino acids
	aa = []
    # List to hold the amino acid combinations
	aa_combi = []
	blosum_data = []
	with open (matrix_file,'r') as file :
		aa = file.readline().strip().split()
		for i in aa :
			for j in aa :
				aa_combi.append((i,j))
		for line in file :
			blosum_data+= line[1:].strip().split()
		for index in range(len(aa_combi)) :
				blosum62[aa_combi[index]] = int(blosum_data[index])
	return blosum62

# Needleman-Wunsch algorithm for sequence alignment
def needleman_wunsch(seq1,seq2,matrix,gap_penalty=-5,return_score=True):
    length_seq1 = len(seq1)
    length_seq2 = len(seq2)
    # Initialize a numpy array with zeros as a temporary placeholder
    score_matrix = np.zeros((length_seq1+1,length_seq2+1),dtype=int)
    # Initialize the score matrix
    for i in range(length_seq1+1) :
        score_matrix[i,0] = i*gap_penalty
    for j in range(length_seq2+1) :
        score_matrix[0, j] = j * gap_penalty
    # Fill in the score matrix
    for i in range(1,length_seq1+1) :
        for j in range(1,length_seq2+1) :
            match = score_matrix[i-1][j-1] + blosum62[(seq1[i-1], seq2[j-1])]
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
    # Traceback to find the optimal alignment        
    alignement1 = ""
    alignement2 = ""
    i = length_seq1
    j = length_seq2

    while i > 0 or j > 0 :
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + blosum62[(seq1[i-1], seq2[j-1])]:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else :
            alignment1 = '-' + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1
    # Return the score or the alignment
    if display_score == False :
        return alignment1, alignment2
    else :
        return score_matrix[length_seq1][length_seq2]

# Initialize a dictionary to store all the alignment scores 
def pairwise_alignment_scores(seq1_id,seq2_id,score,dict_alignment_scores) :
    alignment_scores[(seq1_id,seq2_id)] = score
    return alignment_scores

# Convert alignment scores to distance matrix 
def convert_score_to_distance(alignment_scores_dict,list_keys_sequences_dict) :
    # Initialize a numpy array with zeros as a temporary placeholder
    nb_seq = len(list_keys_sequences_dict)
    distance_matrix = np.zeros((nb_seq,nb_seq))
    # Find the maximum alignment score
    max_score = max(alignment_scores_dict.values())
    # Convert alignment scores to distances
    for i in range(nb_seq):
        for j in range(i+1, nb_seq):
            score = alignment_scores_dict((list_keys_sequences_dict[i], list_keys_sequences_dict[j]))
            distance = max_score - score
            distance_matrix[i, j] = distance
            # Matrix is symmetric
            distance_matrix[j, i] = distance  
    return distance_matrix



# Main program
sequences = read_sequences('SeqFastaCoreG.fasta')
blosum62 = matrix('blosum62.txt')
keys_list = list(sequences.keys())
alignment_scores = {}
for i in range(len(keys_list)) :
    for j in range(i+1,len(keys_list)) :
        score = needleman_wunsch(sequences[keys_list[i]],sequences[keys_list[j]],blosum62,-5,True)
        alignment_scores = pairwise_alignment_scores(keys_list[i],keys_list[j],score,alignment_scores)
        print(f"In progress {i}/{len(sequences)}", end="\r")
distance_matrix = convert_score_to_distance(alignment_scores,keys_list)
