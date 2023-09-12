import os
import sys 
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Read all the sequences from a FASTA file and store them in a dictionary
def read_sequences(fasta_file) :
    """
    Read all the sequences from a FASTA file and store them in a dictionary.
    
    Parameters:
        fasta_file (str): The path to the FASTA file.
        
    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences.
    """
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
    """
    Create the BLOSUM62 matrix stored in a dictionary.
    
    Parameters:
        matrix_file (str): The path to the BLOSUM62 matrix file in .txt format.
        
    Returns:
        dict: A dictionary representing the BLOSUM62 matrix.
    """
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


def needleman_wunsch(seq1,seq2,matrix,gap_penalty=-5,return_score=True):
    """
    Implement the Needleman-Wunsch algorithm for sequence alignment.
    
    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        matrix (dict): The BLOSUM62 matrix.
        gap_penalty (int, optional): The gap penalty. Defaults to -5.
        return_score (bool, optional): Whether to return the score or the alignment. Defaults to True.
        
    Returns:
        int or tuple: Either the alignment score or the aligned sequences, based on `return_score`.
    """
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
    alignment1 = ""
    alignment2 = ""
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
    if return_score == False :
        return alignment1, alignment2
    else :
        return score_matrix[length_seq1][length_seq2]


def pairwise_alignment_scores(seq1_id,seq2_id,score,dict_alignment_scores) :
    """
    Initialize a dictionary to store all the alignment scores.
    
    Parameters:
        seq1_id (str): The ID of the first sequence.
        seq2_id (str): The ID of the second sequence.
        score (int): The alignment score between the two sequences.
        dict_alignment_scores (dict): The dictionary to store the alignment scores.
        
    Returns:
        dict: The updated dictionary with all the alignment scores.
    """
    alignment_scores[(seq1_id,seq2_id)] = score
    return alignment_scores


def convert_score_to_distance(alignment_scores_dict,list_keys_sequences_dict) :
    """
    Convert alignment scores to a distance matrix.
    
    Parameters:
        alignment_scores_dict (dict): A dictionary containing pairwise alignment scores.
        list_keys_sequences_dict (list): A list of sequence IDs.
        
    Returns:
        np.array: A numpy array representing the distance matrix.
    """
    # Initialize a numpy array with zeros as a temporary placeholder
    nb_seq = len(list_keys_sequences_dict)
    distance_matrix = np.zeros((nb_seq,nb_seq))
    # Find the maximum and minimum alignment score
    max_score = max(alignment_scores_dict.values())
    # Convert alignment scores to distances
    for i in range(nb_seq):
        for j in range(i+1, nb_seq):
            score = alignment_scores_dict[list_keys_sequences_dict[i], list_keys_sequences_dict[j]]
            distance = max_score - score 
            distance_matrix[i, j] = distance
            # Matrix is symmetric
            distance_matrix[j, i] = distance  
    return distance_matrix


def upgma(distance_matrix, list_keys_sequences_dict, display_dendrogram=False):
    """
    Perform UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering.
    
    Parameters:
    - distance_matrix (numpy array): The matrix of pairwise distances between sequences.
    - list_keys_sequences_dict (list): List of sequence IDs.
    - display_dendrogram (bool,optiona): Whether to display the dendrogram. Defaults to False.
    
    Returns:
    - A string representation of the resulting tree.
    - The Dendogram if the 3rd argument is True
    """
    # Initialize the number of clusters
    n = len(distance_matrix)
    # Initialize cluster dictionary with sequence IDs
    clusters = {i: [list_keys_sequences_dict[i]] for i in range(n)}
    # Initialize the tree structure representation
    tree_structure = {i: list_keys_sequences_dict[i] for i in range(n)}
    # List to store the linkage matrix, which is needed for dendrogram
    Z = []
    # Main loop: Continue until there's more than one cluster
    while len(clusters) > 1:
        # Initialize minimum distance and the pair of clusters to be merged
        min_distance = float("inf")
        min_indices = (0, 0)
        # Find the pair of clusters with minimum distance
        for i in range(n):
            for j in range(i + 1, n):
                if i in clusters and j in clusters:
                    if distance_matrix[i][j] < min_distance:
                        min_distance = distance_matrix[i][j]
                        min_indices = (i, j)
        # Unpack the cluster indices
        x, y = min_indices
        # Get the actual clusters
        cluster_x, cluster_y = clusters[x], clusters[y]
        # Merge the two clusters
        new_cluster = cluster_x + cluster_y
        new_idx = n
        clusters[new_idx] = new_cluster 
        # Update the distance matrix
        updated_distances = np.zeros((n + 1, n + 1))
        updated_distances[:n, :n] = distance_matrix
        for i in range(n):
            if i in clusters:
                # Compute the new distance
                new_distance = (len(cluster_x) * distance_matrix[x, i] + len(cluster_y) * distance_matrix[y, i]) / (len(cluster_x) + len(cluster_y))
                updated_distances[new_idx, i] = updated_distances[i, new_idx] = new_distance
        distance_matrix = updated_distances
        # Increment cluster index and remove old clusters
        n += 1
        del clusters[x]
        del clusters[y]
        # Append the new merging information to Z
        Z.append([x, y, min_distance, len(new_cluster)])
        # Update the textual/tuple tree representation
        tree_structure[new_idx] = '(%s,%s)' % (tree_structure[x], tree_structure[y])
        del tree_structure[x], tree_structure[y]
    # Display the dendrogram if required
    if display_dendrogram:
        plt.figure(figsize=(20, 7))
        dendrogram(
            np.array(Z),
            leaf_rotation=90,
            leaf_font_size=4,
            leaf_label_func=lambda v: str(list_keys_sequences_dict[v]) if v < len(list_keys_sequences_dict) else '',
            labels=list_keys_sequences_dict
        )
        plt.title('Dendrogram Tree')
        plt.xlabel('Sequence IDs')
        plt.ylabel('Distances')
        plt.show()
    return tree_structure[new_idx]

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
tree_structure = upgma(distance_matrix, keys_list,display_dendrogram=True)

