import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
import argparse


def read_sequences(fasta_file):
    """
    Read all the sequences from a FASTA file and store them in a dictionary.

    Parameters:
        fasta_file (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences.
    """
    # Dictionary to hold the sequences
    sequences = {}
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            # Identifier first line
            if line.startswith(">"):
                query_seq = []
                # Check the FASTA format (uniprot or ncbi)
                end_accession = line[4:].find("|")
                if end_accession == -1:
                    end_accession = line.find(" ")
                    seq_accession = line[1:end_accession]
                else:
                    seq_accession = line[4: end_accession + 4]
            else:
                # Sequence line
                query_seq.extend(list(line.replace(" ", "").strip()))
                sequences[seq_accession] = query_seq
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
    with open(matrix_file, "r") as file:
        aa = file.readline().strip().split()
        for i in aa:
            for j in aa:
                aa_combi.append((i, j))
            for line in file:
                blosum_data += line[1:].strip().split()
            for index in range(len(aa_combi)):
                blosum62[aa_combi[index]] = int(blosum_data[index])
    return blosum62


def needleman_wunsch(seq1, seq2, blosum62, gap_penalty=-5, return_score=True):
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
    score_matrix = np.zeros((length_seq1 + 1, length_seq2 + 1), dtype=int)
    # Initialize the score matrix
    for i in range(length_seq1 + 1):
        score_matrix[i, 0] = i * gap_penalty
    for j in range(length_seq2 + 1):
        score_matrix[0, j] = j * gap_penalty
    # Fill in the score matrix
    for i in range(1, length_seq1 + 1):
        for j in range(1, length_seq2 + 1):
            match = score_matrix[i - 1][j - 1] + \
                blosum62[(seq1[i - 1], seq2[j - 1])]
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
    # Traceback to find the optimal alignment
    alignment1 = ""
    alignment2 = ""
    i = length_seq1
    j = length_seq2
    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and score_matrix[i][j]
            == score_matrix[i - 1][j - 1] + blosum62[(seq1[i - 1], seq2[j - 1])]
        ):
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = seq2[j - 1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        else:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j - 1] + alignment2
            j -= 1
    # Return the score or the alignment
    if return_score:
        return score_matrix[length_seq1][length_seq2]
    else:
        return alignment1, alignment2


def pairwise_alignment_scores(seq1_id, seq2_id, score, dict_alignment_scores):
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
    dict_alignment_scores[(seq1_id, seq2_id)] = score
    return dict_alignment_scores


def convert_score_to_distance(alignment_scores_dict, list_keys_sequences_dict):
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
    distance_matrix = np.zeros((nb_seq, nb_seq))
    # Find the maximum and minimum alignment score
    max_score = max(alignment_scores_dict.values())
    # Convert alignment scores to distances
    for i in range(nb_seq):
        for j in range(i + 1, nb_seq):
            score = alignment_scores_dict[
                list_keys_sequences_dict[i], list_keys_sequences_dict[j]
            ]
            distance = max_score - score
            distance_matrix[i, j] = distance
            # Matrix is symmetric
            distance_matrix[j, i] = distance
    return distance_matrix


def upgma(
    distance_matrix,
    list_keys_sequences_dict,
    display_dendrogram=False,
    dendrogram_file_path=None,
):
    """
    Perform UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering.

    Parameters:
        distance_matrix (numpy array): The matrix of pairwise distances between sequences.
        list_keys_sequences_dict (list): List of sequence IDs.
        display_dendrogram (bool,optiona): Whether to display the dendrogram. Defaults to False.
        dendrogram_file_path (str, optional): Path to save the dendrogram image. Only relevant if display_dendrogram is True.

    Returns:
        tuple: A tuple-based tree structure representing the hierarchical clustering.
        dendrogram : A dendrogram tree if display_dendrogram is True
    """
    # Initialize the number of clusters
    n = len(distance_matrix)
    # Initialize cluster dictionary with sequence IDs
    clusters = {i: [list_keys_sequences_dict[i]] for i in range(n)}
    # Initialize the tree structure representation
    tree_structure = {i: (list_keys_sequences_dict[i],) for i in range(n)}
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
                new_distance = (
                    len(cluster_x) * distance_matrix[x, i]
                    + len(cluster_y) * distance_matrix[y, i]
                ) / (len(cluster_x) + len(cluster_y))
                updated_distances[new_idx, i] = updated_distances[
                    i, new_idx
                ] = new_distance
        distance_matrix = updated_distances
        # Increment cluster index and remove old clusters
        n += 1
        del clusters[x]
        del clusters[y]
        # Append the new merging information to Z
        Z.append([x, y, min_distance, len(new_cluster)])
        # Update the textual/tuple tree representation
        # Update the textual/tuple tree representation
        if type(tree_structure[x]) == tuple and len(tree_structure[x]) == 1:
            tree_structure[x] = tree_structure[x][0]
        if type(tree_structure[y]) == tuple and len(tree_structure[y]) == 1:
            tree_structure[y] = tree_structure[y][0]
        tree_structure[new_idx] = (tree_structure[x], tree_structure[y])
        del tree_structure[x], tree_structure[y]
    # Display the dendrogram if required
    if display_dendrogram:
        plt.figure(figsize=(20, 7))
        dendrogram(
            np.array(Z),
            leaf_rotation=90,
            leaf_font_size=4,
            leaf_label_func=lambda v: str(list_keys_sequences_dict[v])
            if v < len(list_keys_sequences_dict)
            else "",
            labels=list_keys_sequences_dict,
        )
        plt.title("Dendrogram Tree")
        plt.xlabel("Sequence IDs")
        plt.ylabel("Distances")
        plt.savefig(dendrogram_file_path)
        plt.show()
    return tree_structure[new_idx], Z


def multiple_alignment(sequences, blosum62, Z_matrix, list_keys_sequences_dict):
    """
    Perform multiple sequence alignment using hierarchical clustering.

    Parameters:
        sequences (dict): The original sequences to align.
        blosum62 (dict): The BLOSUM62 substitution matrix.
        Z_matrix (list): The linkage matrix from hierarchical clustering.
        list_keys_sequences_dict (list): List of keys referring to sequences.

    Returns:
        list: The list of aligned sequences.
    """
    # Dictionary to hold aligned sequences.
    aligned_sequences = {}
    # List to hold aligned sequences in the order they were processed.
    aligned_sequences_list = []
    # Initialize clusters with individual sequences.
    clusters = {i: [i] for i in range(len(list_keys_sequences_dict))}
    # Loop through linkage matrix to merge clusters and perform alignments.
    for i, (cluster_idx1, cluster_idx2, _, _) in enumerate(Z_matrix):
        # Extract the individual clusters based on indices.
        cluster1 = clusters[cluster_idx1]
        cluster2 = clusters[cluster_idx2]
        # Retrieve either the aligned or the original sequences for these clusters.
        seqs1 = [aligned_sequences[list_keys_sequences_dict[idx]] if list_keys_sequences_dict[idx]
                 in aligned_sequences else sequences[list_keys_sequences_dict[idx]] for idx in cluster1]
        seqs2 = [aligned_sequences[list_keys_sequences_dict[idx]] if list_keys_sequences_dict[idx]
                 in aligned_sequences else sequences[list_keys_sequences_dict[idx]] for idx in cluster2]
        # Case 1: Both clusters have only one sequence; use classic Needleman-Wunsch algorithm.
        if len(cluster1) == 1 and len(cluster2) == 1:
            new_alignment = needleman_wunsch(
                seqs1[0], seqs2[0], blosum62, -5, False)
        # Case 2: Both clusters have more than one sequence; align the clusters.
        elif len(cluster1) > 1 and len(cluster2) > 1:
            new_alignment = needleman_wunsch_ali_to_ali(seqs1, seqs2, blosum62)
        # Case 3: One cluster has only one sequence; align it with the multiple-sequence cluster.
        else:
            if len(cluster1) == 1:
                single_seq = seqs1[0]
                multi_seqs = seqs2
            else:
                single_seq = seqs2[0]
                multi_seqs = seqs1
            new_alignment = needleman_wunsch_ali_to_seq(
                multi_seqs, single_seq, blosum62)
        # Update 'aligned_sequences' with the new alignments.
        new_cluster = cluster1 + cluster2
        for idx, sequence in zip(new_cluster, new_alignment):
            aligned_sequences[list_keys_sequences_dict[idx]] = sequence
        # Update 'clusters' with the new merged cluster.
        new_cluster_idx = len(list_keys_sequences_dict) + i
        clusters[new_cluster_idx] = new_cluster
    # Convert aligned_sequences dictionary to list.
    for key in list_keys_sequences_dict:
        aligned_sequences_list.append(aligned_sequences[key])
    return aligned_sequences_list


def needleman_wunsch_ali_to_seq(
    current_alignment, new_sequence, blosum62, gap_penalty=-5
):
    """
    Perform sequence alignment of a new sequence with an existing multiple sequence alignment
    using a modified Needleman-Wunsch algorithm.

    Parameters:
        current_alignment (list of str): The existing alignment sequences.
        new_sequence (str): The new sequence to align.
        blosum62 (dict): The BLOSUM62 scoring matrix.
        gap_penalty (int, optional): The gap penalty. Defaults to -5.

    Returns:
        list: The updated alignment including the new sequence.
    """
    length_alignment = len(
        current_alignment[0]
    )  # Length of sequences in the existing alignment
    length_new_seq = len(new_sequence)  # Length of the new sequence
    # Initialize the score matrix
    score_matrix = np.zeros(
        (length_alignment + 1, length_new_seq + 1), dtype=int)
    for i in range(length_alignment + 1):
        score_matrix[i, 0] = i * gap_penalty
    for j in range(length_new_seq + 1):
        score_matrix[0, j] = j * gap_penalty
    # Populate the score matrix
    for i in range(1, length_alignment + 1):
        for j in range(1, length_new_seq + 1):
            column = [seq[i - 1] for seq in current_alignment]
            average_score = sum(
                blosum62.get((aa, new_sequence[j - 1]), gap_penalty)
                for aa in column
                if aa != "-"
            ) / len(column)
            match_score = score_matrix[i - 1, j - 1] + average_score
            insert_gap = score_matrix[i - 1, j] + gap_penalty
            align_gap = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(match_score, insert_gap, align_gap)
    # Traceback
    i, j = length_alignment, length_new_seq
    new_aligned_sequence = ""
    while i > 0 or j > 0:
        column = [seq[i - 1] for seq in current_alignment]
        average_score = sum(
            blosum62.get((aa, new_sequence[j - 1]), gap_penalty)
            for aa in column
            if aa != "-"
        ) / len(column)
        if (
            i > 0
            and j > 0
            and score_matrix[i, j] == score_matrix[i - 1, j - 1] + average_score
        ):
            new_aligned_sequence = new_sequence[j - 1] + new_aligned_sequence
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i, j] == score_matrix[i - 1, j] + gap_penalty:
            new_aligned_sequence = "-" + new_aligned_sequence
            i -= 1
        else:
            new_aligned_sequence = new_sequence[j - 1] + new_aligned_sequence
            j -= 1
    # Add the new aligned sequence to the existing alignment
    current_alignment.append(new_aligned_sequence)
    return current_alignment


def needleman_wunsch_ali_to_ali(
    current_alignment, new_alignment, blosum62, gap_penalty=-5
):
    """
    Perform alignment of a new alignment with an existing multiple sequence alignment
    using a modified Needleman-Wunsch algorithm.

    Parameters:
        current_alignment (list of str): Existing alignment sequences.
        new_alignment (list of str): New alignment to be aligned.
        blosum62 (dict): BLOSUM62 scoring matrix.
        gap_penalty (int, optional): Gap penalty, defaults to -5.

    Returns:
        list: Updated alignment including the new alignment.
    """
    # Get the lengths of the current and new alignments
    length_current_ali = len(current_alignment[0])
    length_new_ali = len(new_alignment[0])
    # Initialize the score matrix with zeros
    score_matrix = np.zeros(
        (length_current_ali + 1, length_new_ali + 1), dtype=int)
    # Initialize the first column and row of the matrix
    for i in range(length_current_ali + 1):
        score_matrix[i, 0] = i * gap_penalty
    for j in range(length_new_ali + 1):
        score_matrix[0, j] = j * gap_penalty
    # Populate the score matrix using dynamic programming
    for i in range(1, length_current_ali + 1):
        for j in range(1, length_new_ali + 1):
            current_column = [seq[i - 1] for seq in current_alignment]
            new_column = [seq[j - 1] for seq in new_alignment]
            # Compute average score for matching both columns
            average_current = sum(
                blosum62.get((aa, bb), gap_penalty)
                for aa in current_column
                if aa != "-"
                for bb in new_column
                if bb != "-"
            ) / (len(current_column) * len(new_column))
            match_score = score_matrix[i - 1, j - 1] + average_current
            delete = score_matrix[i - 1, j] + gap_penalty
            insert = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(match_score, delete, insert)
    # Initialize lists to store the final aligned sequences
    aligned_current_alignment = [[] for _ in range(len(current_alignment))]
    aligned_new_alignment = [[] for _ in range(len(new_alignment))]
    # Traceback to find the optimal alignment
    i, j = length_current_ali, length_new_ali
    while i > 0 or j > 0:
        current_column = [seq[i - 1] for seq in current_alignment]
        new_column = [seq[j - 1] for seq in new_alignment]
        # Calculate average score for the current cell
        average_current = sum(
            blosum62.get((aa, bb), gap_penalty)
            for aa in current_column
            if aa != "-"
            for bb in new_column
            if bb != "-"
        ) / (len(current_column) * len(new_column))
        # Check which move is optimal and update accordingly
        if (
            i > 0
            and j > 0
            and score_matrix[i, j] == score_matrix[i - 1, j - 1] + average_current
        ):
            for seq_idx in range(len(current_alignment)):
                aligned_current_alignment[seq_idx].append(
                    current_column[seq_idx])
            for seq_idx in range(len(new_alignment)):
                aligned_new_alignment[seq_idx].append(new_column[seq_idx])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i, j] == score_matrix[i - 1, j] + gap_penalty:
            for seq_idx in range(len(current_alignment)):
                aligned_current_alignment[seq_idx].append(
                    current_column[seq_idx])
            for seq_idx in range(len(new_alignment)):
                aligned_new_alignment[seq_idx].append("-")
            i -= 1
        else:
            for seq_idx in range(len(current_alignment)):
                aligned_current_alignment[seq_idx].append("-")
            for seq_idx in range(len(new_alignment)):
                aligned_new_alignment[seq_idx].append(new_column[seq_idx])
            j -= 1
    # Reverse and join each aligned sequence to get the final alignment
    for seq_idx in range(len(current_alignment)):
        aligned_current_alignment[seq_idx] = "".join(
            reversed(aligned_current_alignment[seq_idx])
        )
    for seq_idx in range(len(new_alignment)):
        aligned_new_alignment[seq_idx] = "".join(
            reversed(aligned_new_alignment[seq_idx])
        )
    return aligned_current_alignment + aligned_new_alignment


# Main program
def main(args):
    """
    Main function to perform multiple sequence alignment.

    Parameters:
        args: ArgumentParser object containing command-line arguments
    """
    # Verify the existence of the input files
    if not os.path.exists(args.seq_file) or not os.path.exists(args.blosum_file):
        print("Error: Specified file paths do not exist.")
        exit(1)
    # Create the output directory if it doesn't exist
    if not os.path.exists("Outputs"):
        os.makedirs("Outputs")
    # Generate a unique counter for output files
    counter = 1
    while os.path.exists(os.path.join("Outputs", f"MultiAlignmt_{counter}.txt")):
        counter += 1
    # Use this counter to name the output files
    multialignmt_file_path = os.path.join(
        "Outputs", f"MultiAlignmt_{counter}.txt")
    dendrogram_file_path = os.path.join("Outputs", f"Dendrogram_{counter}.png")
    # Read sequences and BLOSUM62 matrix from input files
    sequences = read_sequences(args.seq_file)
    blosum62 = matrix(args.blosum_file)
    # Get keys of sequences dictionary to prepare for alignment
    keys_list = list(sequences.keys())
    alignment_scores = {}
    # Compute alignment scores for each sequence pair
    for i in range(len(keys_list)):
        for j in range(i + 1, len(keys_list)):
            score = needleman_wunsch(
                sequences[keys_list[i]], sequences[keys_list[j]], blosum62
            )
            alignment_scores = pairwise_alignment_scores(
                keys_list[i], keys_list[j], score, alignment_scores
            )

    distance_matrix = convert_score_to_distance(alignment_scores, keys_list)
    tree_structure, Z_matrix = upgma(
        distance_matrix,
        keys_list,
        display_dendrogram=args.display_dendrogram,
        dendrogram_file_path=dendrogram_file_path,
    )
    result_alignmt = multiple_alignment(
        sequences, blosum62, Z_matrix, keys_list)

    # Saving the files
    with open(multialignmt_file_path, "w") as f:
        for seq in result_alignmt:
            f.write(seq + "\n")

    if args.display_dendrogram:
        print(
            f"Your multiple alignment file '{os.path.basename(multialignmt_file_path)}' and your dendrogram '{os.path.basename(dendrogram_file_path)}' are located here: {multialignmt_file_path}"
        )
    else:
        print(
            f"Your multiple alignment file '{os.path.basename(multialignmt_file_path)}' is located here: {multialignmt_file_path}"
        )


if __name__ == "__main__":
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="Perform multiple sequence alignment")
    # Define command-line arguments
    parser.add_argument(
        "--seq_file", required=True, help="Path to the FASTA file containing sequences"
    )
    parser.add_argument(
        "--blosum_file", required=True, help="Path to the BLOSUM62 matrix file"
    )
    parser.add_argument(
        "--display_dendrogram", action="store_true", help="Display the dendrogram"
    )
    # Parse command-line arguments
    args = parser.parse_args()
    # Call the main function
    main(args)
