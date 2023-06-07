#!usr/bin/env

'''
This script performs the following tasks:
    1. Counts the frequency of each nucleotide in each sequence
        in a provided FASTA file.

    2. Counts the frequency of each nucleotide across all sequences
        in a provided FASTA file.

    3. Calculates the frequency of transition and
        transversion mutations between pairwise comparison of
        DNA seqeunces in a provided FASTA file.

    4. Calculates the overall ratio of transition to transverion
        mutations for each pairwise comparison of sequences in a
        provided FASTA file.

    5. Records all nucleotide frequencies, mutation counts, and
        ratios in a single text file: "FrequencyData.txt".

Before beginning, it is important to set the working directory to
the folder containing your FASTA file. This can be done by clicking
Run -> Edit Configurations and manually selecting a new working directory
in PyCharm.

Users only need to edit one variable in the following script:
input_file (Line 37) should be updated to reflect the name
(including extension) of the target FASTA file.
'''


#Import relevant packages
from Bio import SeqIO

#Load input file and update to match your filename/path
input_file = open('lab4_termite.fasta', 'r')

#Create output file to store frequencies
output_file = open('FrequencyData.txt', 'w') #This file will be stored in the same directory as input_file
output_file.write('Accession\tA\tC\tG\tT\tN\tGaps\n') #Initializing headers for output file

#Initialize variables to hold counts of A,C,G,T, and N across all sequences in the FASTA file
total_A = 0
total_C = 0
total_G = 0
total_T = 0
total_N = 0
total_gap = 0

#For loop to calculate nucleotide frequencies
ids = [] #Initialize a variable to store accession numbers of each sequence for later identification
seqs = [] #Initialize a variable to store each DNA sequence without headers

for record in SeqIO.parse(input_file, "fasta"):
    acc_num = record.id #Store each accession number as an identifier

    #Count and store frequency of A,C,G,T, and N in each sequence in the FASTA file
    A_count = record.seq.count('a')
    C_count = record.seq.count('c')
    G_count = record.seq.count('g')
    T_count = record.seq.count('t')
    N_count = record.seq.count('n')
    gap_count = record.seq.count('-')
    #NOTE: Change capitalization of nucleotides to match the format of your FASTA sequences

    #Incrementally add each sequence's nucleotide count to the total counters
    total_A += A_count
    total_C += C_count
    total_G += G_count
    total_T += T_count
    total_N += N_count
    total_gap += gap_count

    #Add all IDs and sequences to ids and seqs respectively
    ids.append(record.id)
    seqs.append(record.seq)

    #Add each value as a new row to the output file, corresponding to the order of the headers in Line 41
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\n' % \
    (acc_num, A_count, C_count, G_count, T_count, N_count, gap_count)
    output_file.write(output_line)

#Write final output line to include total frequencies of each nucleotide.
final_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\n' % \
("     TOTAL", total_A, total_C, total_G, total_T, total_N, total_gap)
output_file.write(final_line)
#NOTE: The above code is outside of the for loop to allow each sequence frequency to be added to the total

#Write a line to describe the strings of transition/transverion frequency and ratio data
description_lines = "\nEach line below shows the ratio of transition to transversion mutations for " \
                   "every pairwise comparison of DNA sequences within the provided FASTA file. " \
                   "Raw frequencies of each mutation are given in parentheses.\n"
output_file.write(description_lines)

#Write a function that counts transitions and transversions for each pairwise combination of DNA sequences
def count_mutations(DNAseqs, IDlist): #The function accepts a list of sequences and a list of IDs for matching
    #Define the sets of purines and pyrimidines
    purines = {"a", "g"}
    pyrimidines = {"c", "t"}

    #Initialize holder lists
    mutations=[] #List to hold counts of transition and transversion for each pairwise combination
    pairs = [] #List to hold accession numbers of each sequence used in the corresponding pairwise combination

    #Create a for loop to generate each pairwise sequence comparison in seqs
    for i in range(len(DNAseqs)):
        for j in range(i+1, len(DNAseqs)):
            seq1 = DNAseqs[i]
            seq2 = DNAseqs[j]

            pairs.append((IDlist[i], IDlist[j])) #Add the pairwise-corresponding IDs to the "pairs" list

            #Initialize transition and transversion frequencies
            transition_freq = 0
            transversion_freq = 0

            #For loop to calculate frequencies of each type of mutation at each position in the sequences
            for k in range(len(seq1)):
                if seq1[k] != seq2[k]:
                    if seq1[k] in purines and seq2[k] in purines:
                        transition_freq +=1
                    elif seq1[k] in pyrimidines and seq2[k] in pyrimidines:
                        transition_freq +=1
                    elif seq1[k] in purines and seq2[k] in pyrimidines:
                        transversion_freq +=1
                    elif seq1[k] in pyrimidines and seq2[k] in purines:
                        transversion_freq += 1

            mutations.append((transition_freq, transversion_freq)) #Add frequencies as pairs to "mutations"

    #Set the output of the function as the set of transition/transversion frequencies and corresponding IDs
    return mutations, pairs

#Now, we run this function using the previously populated "seqs" (Line 77) and "ids" (Line 76) variables
results = (count_mutations(seqs, ids))

mut_counts = results[0] #Store just Transition/Transverion Pairs in mut_counts
corr_id = results[1] #Store just the corresponding IDs in corr_id

#For loop to combine this information and calculate ratios
for i in range(len(mut_counts)):
    count_set = mut_counts[i]
    id_set = corr_id[i]
    count_transition = count_set[0]
    count_transversion = count_set[1]

    if count_transversion != 0: #Ensuring non-zero denominator for ratio
        raw_ratio = count_transition/count_transversion #Calculate ratio
        rounded_ratio = ('%.2f' % raw_ratio)  # Round to two decimal places
    else:
        rounded_ratio = "No transversions!" #Prevents errors when trying to divide by zero


    #Write the transition and transverion frequencies and ratios as a single line
    ratio_line = "\nThe ratio of transitions (%i) to transversions (%i) for %s vs %s is: %s\n" \
                 % (count_transition, count_transversion,
                    id_set[0], id_set[1],
                    str(rounded_ratio))

    #Write each line to the existing output file
    output_file.write(ratio_line)

#Close all files
output_file.close()
input_file.close()