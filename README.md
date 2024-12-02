# Introduction
The purpose of this repository is to describe the steps involved in analyzing the RuvB-like gene. This process was done throughout the semester but will be broken down and summarized into the key steps done throughout various labs.

#  Steps
## Lab 3: Homolog Identification Using BLAST
##### Objective:
To identify homologs of the RuvB-like gene family across selected species by performing BLAST searches against a proteome database. This lab prepares the foundational dataset for downstream analysis by filtering high-confidence matches based on e-values.

##### Description:
The lab involves creating a local BLAST database from proteomes, querying the RuvB-like gene to find homologous sequences, and filtering results to retain significant hits. These homologs will serve as input for multiple sequence alignment and phylogenetic analysis in later labs.

##### Lab clone was created with and directory navigated to by these commands:
```
git clone https://github.com/Bio312/lab03-$MYGIT

cd lab03-$MYGIT
```

##### Creating Blast Database
I first began by unzipping the downloaded proteomes from the copied lab filed and combinining them while in the home directory for lab 3.
```
gunzip proteomes/*.gz

cat proteomes/*.faa > allprotein.fas
```
Next, the blast database was created:
```
makeblastdb -in allprotein.fas -dbtype prot
```
##### RuvB-like gene Analysis
To begin the analysis of my gene, we start by creating a separate directory for my gene. I started this in lab 3, so the code I put into the VS Code terminal was:
```
mkdir ~/lab03-$MYGIT/ruvB-like
```
Now that the separate directory was created, I had to navigate to this directory in order organize all the data on my gene into this directory. I did this by using the cd command which changes your current directory. The input in terminal was:
```
cd ~/lab03-$MYGIT/ruvB-like
```
Now that we are in the correct directory, the first step is to download the query protein that will be worked on. To do this I used the input:
```
ncbi-acc-download -F fasta -m protein "NP_003698.1"
```
NP value represents the RefSeq protein accession number. These are unique identifiers assigned to protein sequences within the RefSeq (Reference Sequence) database.

This value was assignmed to me at the beginning of the semester, the code that I just ran downloaded the query sequence of this protein and the objective currently is to identify what family this gene belonged to. Prior information was found up to this point by searching for the NP on gen bank.

A blast search was then performed for my gene in order to identify homologs of the query protein: with the input:
``` 
blastp -db ../allprotein.fas -query NP_000549.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out globins.blastp.detail.out
```
The blast results were next filtered to only include significant macthes with an e-value of 1e-30
```
awk '{if ($6<1e-30) print $1 }' ruvB-like.blastp.detail.out > ruvB-like.blastp.detail.filtered.out
```
To count the filtered hits:
```
wc -l ruvB-like.blastp.detail.filtered.out
```
To find the number of homologs in each species:
```
grep -o -E "^[A-Z]\.[a-z]+" ruvB-like.blastp.detail.filtered.out | sort | uniq -c
```
All of the data files for these lab 3 inputs can be found under "Lab 3 homolog data"
## Lab 4: Sequence Alignment and Statistical Analysis
##### Objective:
To align the homologous sequences of the RuvB-like gene family and compute alignment statistics, including conservation and sequence identity, using tools like MUSCLE and AlignBuddy.

##### Description:
This lab focuses on extracting homolog sequences from the filtered BLAST results and performing multiple sequence alignment. Alignment quality and sequence similarity metrics are calculated to evaluate the evolutionary conservation of the gene family.
##### Start by cloning the lab and creating a directory for the RuvB-like gene family
```
git clone https://github.com/Bio312/lab04-$MYGIT

mkdir ~/lab04-$MYGIT/RuvB-like
```
##### Move into the directory
```
cd ~/lab04-$MYGIT/RuvB-like
```
##### Extract RuvB-like homolog sequences from the BLAST Output file
```
seqkit grep --pattern-file ~/lab03-$MYGIT/ruvB-like/ruvB-like.blastp.filtered.out ~/lab03-$MYGIT/allprotein.fas > ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.fas
```
##### Perform multiple sequence alignment using muscle
```
muscle -align ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.fas -output ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas
```
##### View alignment in alv
```
alv -kli ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas
```
##### Generate PDF alignment using R
```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas
```
##### Calculate alignment statistics using alignbuddy
```
alignbuddy -pi ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas
```
##### Calculate the length of the alignment after removing any column with gaps:
```
alignbuddy -trm all  ~/lab04-$MYGIT/RuvB-like/RuvB-like.homolog
```
##### Calculate the length of the alignment after removing invariant (completely conserved) positions:

```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs
```
##### Calculate the average percent identity among all sequences in the alignment using t_coffee:
```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas -output sim
```
##### Repeat calculating the average percent identity using alignbuddy.
```
 alignbuddy -pi ~/lab04-$MYGIT/globins/globins.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
## Lab 5: Phylogenetic Tree Construction
##### Objective:
To construct a phylogenetic tree for the RuvB-like gene family using IQ-TREE and visualize the tree with midpoint rooting and other graphical representations.

##### Description:
The lab involves refining aligned sequences from Lab 4, running IQ-TREE to generate a maximum likelihood phylogenetic tree, and applying midpoint rooting for interpretation. The resulting tree provides insights into the evolutionary relationships of the gene family.

##### Clone the lab and create the directory
```
git clone https://github.com/Bio312/lab05-$MYGIT

mkdir -p ~/lab05-$MYGIT/RuvB-like
```
##### Navigate to the gene family directory
```
cd ~/lab05-$MYGIT/RuvB-like
```
##### Clean and copy the alignment file from lab04
```
sed 's/ /_/g' ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.al.fas | \
seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas
```
##### Run IQ-TREE with ultrafast bootstrap approximation
```
iqtree -s ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas -bb 1000 -nt 2
```
##### Display the tree in ASCII format
```
nw_display ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas.treefile
```
##### Generate an unrooted tree using R script
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R \
~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas.treefile \
~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas.treefile.pdf 0.4 15
```
##### Midpoint root the tree
```
gotree reroot midpoint -i ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.fas.treefile \
-o ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile
```
##### Order and display the midpoint-rooted tree
```
nw_order -c n ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile | nw_display -
```
##### Create graphical representations of the tree
```
nw_order -c n ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile | \
nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.svg
```
##### Convert SVG to PDF
```
convert ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.svg \
~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.pdf
```
##### Generate a cladogram for comparison
```
nw_order -c n ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile | nw_topology - | \
nw_display -s -w 1000 > ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.midCl.treefile.svg
```
##### Convert the cladogram SVG to PDF
```
convert ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.midCl.treefile.svg \
~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.midCl.treefile.pdf
```
## Lab 6: Gene Tree Reconciliation
##### Objective:
To reconcile the RuvB-like gene tree with a species tree to identify gene duplication and loss events using Notung and visualize the reconciled tree.

##### Description:
The reconciled tree highlights the evolutionary history of the RuvB-like gene family, including duplication and loss events across species. This lab integrates phylogenetic and species tree data to infer functional and evolutionary significance.
##### Clone the lab
```
git clone https://github.com/Bio312/lab06-$MYGIT
```
##### These commands create a new environment (my_python27_env) with Python 2.7 and installed the necessary software like ete3 to handle phylogenetic trees.
```
mamba create -n my_python27_env python=2.7
conda activate my_python27_env
mamba install -y ete3
```
##### Create a folder for the RuvB-like gene family and copy the relevant gene tree from Lab 5 to Lab 6.
```
mkdir ~/lab06-$MYGIT/RuvB-like
cp ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile ~/lab06-$MYGIT/RuvB-like
```
##### Reconcile the RuvB-like gene tree with the species tree using Notung and save the resulting reconciled tree image and event summary in the RuvB-like directory.
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/RuvB-like/
```
##### Convert the Notung output into a RecPhyloXML format for visualizing gene trees reconciled within species trees.
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.rec.ntg --include.species
```
##### Generate an SVG graphic showing the reconciled gene tree within the species tree using thirdkind.
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.rec.svg
```
##### Convert the thirdkind SVG output into a PDF for easier viewing.
```
convert -density 150 ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile.rec.pdf
```
## Lab 8: Domain Analysis and Integration
##### Objective:
To identify and visualize conserved protein domains in homologous sequences of the RuvB-like gene family using RPS-BLAST and integrate domain information with the phylogenetic tree.

##### Description:
This lab identifies Pfam domains within the RuvB-like sequences and maps them onto the phylogenetic tree. The analysis reveals conserved functional elements within the gene family, providing insights into domain evolution and specialization.

##### Clone the lab and create a directory for the gene
```
git clone https://github.com/Bio312/lab08-$MYGIT

mkdir ~/lab08-$MYGIT/RuvB-like
```
##### Switch to the appropriate directory
```
cd ~/lab08-$MYGIT/RuvB-like
```
##### This command uses sed to find and remove all instances of the stop codon (*) in the file RuvB-like.homologs.fas. The cleaned file is saved to the RuvB-like directory.
```
sed 's/*//' ~/lab04-$MYGIT/RuvB-like/RuvB-like.homologs.fas > ~/lab08-$MYGIT/RuvB-like/RuvB-like.homologs.fas
```
##### This command runs RPS-BLAST on the cleaned sequences in RuvB-like.homologs.fas, searching the Pfam database. The output is saved in a tab-delimited format with columns for the query ID, sequence length, domain positions, e-value, and domain title. An e-value cutoff of 1e-10 is used for more precise results.
```
rpsblast -query ~/lab08-$MYGIT/RuvB-like/RuvB-like.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/RuvB-like/RuvB-like.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
##### This command copies the midpoint-rooted tree file from the Lab 5 directory into your working directory for Lab 8. This tree will be used for plotting domain annotations.
```
cp ~/lab05-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile ~/lab08-$MYGIT/RuvB-like
```
##### This command runs an R script that generates a visualization of the phylogenetic tree alongside the predicted Pfam domains. The --vanilla flag ensures R runs without saving or restoring workspaces. The script takes three inputs: the phylogenetic tree, the RPS-BLAST output, and the output PDF file.
```
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/RuvB-like/RuvB-like.homologsf.al.mid.treefile ~/lab08-$MYGIT/RuvB-like/RuvB-like.rps-blast.out ~/lab08-$MYGIT/RuvB-like/RuvB-like.tree.rps.pdf
```
##### This command extracts the sixth column (Pfam domain titles) from the RPS-BLAST output using cut. It then sorts the results and counts the occurrences of each unique domain using uniq -c.
```
cut -f 6 ~/lab08-$MYGIT/RuvB-like/RuvB-like.rps-blast.out | sort | uniq -c
```
##### This command calculates the length of each domain by subtracting the start position from the end position. It prints the gene ID and the domain length, then sorts the list by domain length in descending order.
```
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/RuvB-like/RuvB-like.rps-blast.out | sort -k2nr
```
##### Explanation: This command extracts the first column (query sequence ID) and fifth column (e-value) from the RPS-BLAST output file. It helps compare e-values across different sequences to identify which gene has the best (smallest) e-value.
```
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/RuvB-like/RuvB-like.rps-blast.out
```
## Conclusion
This repository documents the full workflow for analyzing the RuvB-like gene family, providing all necessary commands, outputs, and explanations to enable reproducibility.

All necessary data is included in its respective lab folder in this repository, including all allignments, reconciliations and gene trees.