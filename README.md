# working-with-fasta-files
The Python script in this repo `Subramaniam.py` performs the following tasks:

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
        ratios in a single text file: `FrequencyData.txt`.

Before using this Python script, it is important to set the working directory to
the folder containing your FASTA file. This can be done by clicking
Run -> Edit Configurations and manually selecting a new working directory
in PyCharm.

Users only need to edit one variable in `Subramaniam.py`:
`input_file` in Line 37 should be updated to reflect the name
(including extension) of the target FASTA file.
