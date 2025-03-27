# Data Collection (Oct. 2024 - Nov. 2024)
### 1. Literature Review & NCBI
Based on the literature review, marine bacteria were listed up on excel file with the corresponded enzyme sequences from NCBI.  
By using the InterPro (InterPro 104.0), predict the presence and location of signal peptides, which is a short sequence of AA at the beginning of protein. If it is located in the extracellular region, it indicates that the enzyme will be secreted from the cell to function in the outside, which means the extracellular localization allows the enzymes to act on substrates in the environment around the cell.
  1. Copy paste of the protein sequence to the box 'Search by sequence'
  2. Press search
  3. From the result, it is possible to finds the conserved domains and signal peptides in the protein sequence.

From NCBI, checked the gene sequences to see if the listed bacteria contain any additional enzymes besides the one mentioned in each paper during literature review.

### 2. BLASTp
Identify protein sequences with the high identity of sequence to those in the listed bacteria.  
![Screen Shot 2025-03-27 at 6 37 33 PM](https://github.com/user-attachments/assets/1205a98b-5857-46aa-a8e9-c6b369024217)  
The result shows how similar a given protein sequence (query) is to other protein sequences in a database. It describe which proteins in a database are most similar to your query. 
For each row of sequence from bacteria, you can find E-value (for significance), Identity (for similarity), and the alignment (to understand how they match).  
![Screen Shot 2025-03-27 at 7 22 35 PM](https://github.com/user-attachments/assets/49c74cea-8be5-48f5-bf25-b75ccaf67dba)  
The high matched results were recorded on excel file together with collected data from NCBI previously.

### 3. BacDive
It confirm if the bacteria with available type strains also contain the proteins identified in the BLAST results.  
It also provide properties of specific bacteria (habitat, growth T, pathogenicity).  

# Phylogentic trees (Nov. 2024 - Dec. 2024)
### 1. Convert CSV file (.csv) to FASTA file (.fa)
All of the collected data was recorded on excel file (.csv). To produce the phylogenetic trees, first, the python is used to generate fasta files of the sequences.  
Version 
```python
# Python code starts here
import pdb
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
file = open("Bacteria_genetic_comparison.csv", 'r', encoding='UTF8')
header = file.readline()
count = 1
seq = ""
seq_start = 0
all_records = []
for line in file:
    line = line.rstrip()
    if line:
        line = line.split(',')
        if seq_start == 1:
            subseq = line[0].replace(' ', '').upper().rstrip()
            if subseq[-1] == '"':
                seq_start = 2
                subseq = subseq[:-1]
            seq += subseq

        if line[0] == str(count):
            seq_start = 1
            if line[1]:
                name = line[1]
            seq_id = line[2] + line[3].replace('"',' ').rstrip()
            seq = "" #protein name
            print(line)

        if seq_start == 2:
            record = SeqRecord(
                Seq(seq),
                id=seq_id,
                name=name,
                description=name,
            )
            all_records.append(record)
            seq_start = 0
            count += 1

with open('protein.fa', 'w') as handle:
    SeqIO.write(all_records, handle, 'fasta')

### 2. FastTree 
It was used to construct phylogenetic trees from multiple sequence alignments







