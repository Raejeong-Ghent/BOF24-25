# Data collection (Oct. 2024 - Nov. 2024)
### 1. Literature Review & NCBI
Based on the literature review, marine bacteria were listed up manually on excel file with the corresponded enzyme (Chitinase, Deacetylase, LPMO) sequences from NCBI.  +++ HOW
By using the InterPro (InterPro 104.0), predict the presence and location of signal peptides, which is a short sequence of AA at the beginning of protein. If it is located in the extracellular region, it indicates that the enzyme will be secreted from the cell to function in the outside, which means the extracellular localization allows the enzymes to act on substrates in the environment around the cell.
  1. Copy paste of the protein sequence to the box 'Search by sequence'
  2. Press search
  3. From the result, it is possible to finds the conserved domains and signal peptides in the protein sequence.

From NCBI, checked the gene sequences to see if the listed bacteria contain any additional enzymes besides the one mentioned in each paper during literature review.

### 2. BLASTp Searches
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
_Version 3.12.7_
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
```

### 2. FastTree 
It was used to construct phylogenetic trees from multiple sequence alignments. 
Terminal-based workflow was used to complete the tasks. 
The following terminal commands were executed using these versions:

- gawk: _GNU Awk 5.3.1_
- MAFFT: _v7.526_
- FastTree: _v2.1.11_

```shell
gawk '/^>/ {flag=0} /bacteria_name/ {flag=1} flag {print}' protein-longnames.fa > bacteria_name_proteins.fa
mafft bacteria_name_proteins.fa > bacteria_name_proteins-aln.fa
FastTree bacteria_name_proteins-aln.fa > bacteria_name_proteins.tree
```
### 3. FigTree
_Version 1.4.4_
The produced file (.tree) could be opened in FigTree and it is possible to edit them such as labeling, coloring and scaling, etc.  

<ins>Results</ins>  
In total, 5 different phylogenetic trees were produced, each corresponding to specific enzyme. Upcoming object is to select marine bacteria with the greatest diversity. Thus, prioritized bacterial species that appear in more than one phylogenetic trees for the different enzymes. In other words, bacteria that produce more than one specific enzyme were selectied.  

However, constructing phylogenetic tree uing manually selected dataset may not be entirely reliable, as we might have missed unannotated chitin-realted proteins or proteins annotated with unexpected names. This limitation led to explore the use of whole genome sequences.  

# Select marine bacteria based on whole genome sequence (Dec. 2024 - Jan. 2025)
### 1. Whole genome sequence collection
The whole genome sequences of bacteria from previously recorded list are searched easily from NCBI. (Example followed)  

<img width="744" alt="Screen Shot 2025-03-29 at 2 50 38 PM" src="https://github.com/user-attachments/assets/60187cd6-100e-4987-b3c4-ddaa9f71a751" />

++++ HOW

Whole genome sequence are too large to use online BLAST to identify protein sequences with the high identity. Thus, by suing local BLAST, compared whole genome sequences with previously constructed dataset.

### 2. Local BALST
_Version blast 2.16.0_  
Filtering rules are applied:  
  1. Removing sequences shorter or equal to 20 baseparis (bp)
  2. Ignoring sequences with 100% identity, since they are identical to known ones which are not our focus obviously.
  3. Filtering out results with a high number of mismatches.
     * Initially, arbitrarily set the mismatch limit to 50%. If mismatch is higher than that, filtering them out.
     * For proper filtration, use histogram to examine the distrubution of mismatches data and adjust the threshold accordingly. +++

```shell
awk -F'\t' '
BEGIN { OFS="\t" }
NR == 1 { print; next } # Print the header
$4 >= 10 && $3 < 100 && $5 <= ($4 / 2) { print }
' /mnt/data/Aquimarina_amphilecti_DSM25232.txt > filtered_output.txt
```


























