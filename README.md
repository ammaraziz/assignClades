## assignClades.py

Assign clades to HA influenza sequences. Heavily reliant on nextstrain seasonal influenza tools and scripts. 
    
### Usage 

assignClades.py -s seq.fasta -l h3n2 -b outputName
	-s, --sequence          Path to input sequence [YourSequences.fasta]
	-l, --lineage           Lineage of input strains [h1n1, h3n2, vic, yam]
	-b, --batchName         The batch name

### Example usage

assign_clades.py --sequences Batch999_01Jan20.fasta --lineage h1n1 --batchName Batch999_results.txt
assign_clades.py -s Batch999_01Jan20.fasta -l h1n1 -b Batch999_results.txt

### Output

batchName_clades.txt        complete clade provenance
batchName_results.txt       current clade and vaccine result 

Clade defintions are stored in /config/{lineage}.tsv.

Special thanks to the nextstrain crew for their amazing work on influenza.