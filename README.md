# Hg-MATE-Db: Hg-cycling Microorganisms in Aquatic and Terrestrial Ecosystems Database
Microorganisms play a significant role in regulating the form and fate of mercury (Hg) in aquatic and terrestrial ecosystems. Microbes with the hgcAB gene pair can produce a more toxic, and bioaccumulative form of Hg, methylmercury (MeHg). Microbes that possess the mer operon can demethylate and/or reduce Hg species as part of a detoxification mechanism. Improved techniques for capturing hgcAB and mer presence and diversity are necessary for identifying the major microbial players in environmental Hg cycling. The primary goal of the database Hg-MATE is to provide an up-to-date collated resource of Hg-cycling genes from pure culture and environmental microbial genomes and meta-omic datasets. The current <b>version 1</b> contains an hgcAB dataset with resources for identifying key microbial producers of the toxin MeHg. Future versions  will include a mer gene dataset, which will contain resources for identifying genes of the mer-operon that encode for demethylation of organomercurials (merB), reduction of inorganic Hg(II) (merA), as well as operon regulation (merR), and Hg transport across the cell (merTPC).

## Resources
The version 1 of the Hg-MATE database (v1.01142021) contains:
* <b> A catalogue of 1053 HgcAB amino acid sequences </b>. There, HgcAB amino acid sequences are categorized into four types depending on whether they were encoded in pure culture/environmental microbial isolates (<i>ISO</i>), single-cell genome sequences (<i>CEL</i>), metagenome-assembled genomes (<i>MAGs</i>) and environmental meta-omic contig (<i>CON</i>). Included in the database are amino acid sequences of HgcA, HgcB, and concatenated HgcA and HgcB. If hgcB is not co-localized with hgcA in the genome and/or cannot be identified, then ‘na’ will be listed in the ‘HgcB’ sequence column. We collated the HgcAB databases from <a href="https://doi.org/10.3389/fmicb.2020.541554" target="_blank"><b>Gionfriddo et al. 2020</b></a> and <a href="https://doi.org/10.1128/mSystems.00299-20" target="_blank"><b>McDaniel et al. 2020</b></a> and added HgcAB amino acid sequences pulled from three public data repositories: <a href="https://www.ncbi.nlm.nih.gov/genbank/" target="_blank"><b>NCBI GenBank</b></a>, <a href="https://gold.jgi.doe.gov/" target="_blank"><b>JGI-IMG GOLD</b></a> and <a href="https://gtdb.ecogenomic.org/" target="_blank"><b>GTDB release 89</b></a> obtained on 23 October 2020. HgcAB amino acid sequences were identified in these databases by hmmsearch with HgcA and HgcB HMM profiles from <a href="https://doi.org/10.3389/fmicb.2020.541554" target="_blank"><b>Gionfriddo et al. 2020</b></a>. Other resources generated from pure culture/environmental isolates, single-cell genome sequences, and metagenome-assembled genomes (‘ISOCELMAG’) include 
* <b>FASTA files </b> containing amino acid sequences of HgcA (‘_HgcA.fas’), HgcB (‘_HgcB.fas’), and concatenated HgcA-HgcB sequences (‘_Hgc.fas’). FASTA files with either unaligned and aligned (msa) amino acid sequences are provided.
* <b>Hidden Markov models</b> containing amino acid sequences of HgcA (‘_HgcA.hmm’), HgcB (‘_HgcB.hmm’), and concatenated HgcA-HgcB sequences (‘_Hgc.hmm’). 
* <b>Reference packages</b> that can be used to identify and classify: 1) the cap-helix encoding region of HgcA (‘_HgcA_CH.refpkg‘), for example in Desulfovibrio desulfuricans ND132, this encompasses the CdhD-like encoding region, sites ~37-156 of HgcA (https://www.uniprot.org/uniprot/F0JBF0); 2) full HgcA (‘_HgcA_Full.refpkg‘); and 3) concatenated HgcA and HgcB (‘_HgcA-HgcB.refpkg‘). Each reference package contains sequence alignments, HMM model, phylogenetic tree, and NCBI taxonomy.


## Recommended Usage
We recommend using the resources provided by Hg-MATE (HMM profiles and reference packages) to detect and taxonomically identify hgcAB genes in metagenomes, metatranscriptomes and MAGs. The consensus protocol provided in <a href="https://www.biorxiv.org/content/10.1101/2022.03.14.484253v1.abstract" target="_blank"><b>Capo et al. 2022</b></a> is recommended for cross-comparaison between metagenomics studies. Intermediate files generated from raw data (fastq files) that are required to detect, count and identify the genes are : proteins.faa (obtained by prodigal or prokka) and counts.tsv (obtained by featureCounts). See <a href="https://www.biorxiv.org/content/10.1101/2022.03.14.484253v1.abstract" target="_blank"><b>Capo et al. 2022</b></a> for detailed description. The reference packages can be used for phylogenetic analysis of HgcA(B) sequences from amplicon sequencing or meta-omic datasets. 

An example workflow for how to use the Hg-MATE-Db reference packages for identifying and classifying HgcA(B) sequences using <a href="http://hmmer.org/" target="_blank"><b>HMMER</b></a> and <a href="https://doi.org/10.1186/1471-2105-11-538" target="_blank"><b>pplacer </b></a>.

* Copy your files (sample_1.fastq and sample_2.fastq) in the folder marky-coco. Use 'gunzip' if your fastq are in fastq.gz format.
```
cp /remote/folder/sample_1.fastq .
cp /remote/folder/sample_2.fastq .
```
* Activate the conda environment
```
conda activate coco
```
* Run the marky script
```
bash marky.sh sample
```

Test files are here https://figshare.com/articles/online_resource/Test_files_for_marky-coco_pipeline/19221213
```
wget https://figshare.com/ndownloader/articles/19221213/versions/1
unzip 1
bash marky.sh MG01
```

## SLURM USAGE
* Run the marky_to_slurm file. You can modify the requested amount of time in the file.
```
sbatch marky_to_slurm.sh sample
```


## ADVANCED USAGE
Standards input files are paired-end fastq files (sample_1.fastq & sample_2.fastq). Alternately, intermediate files can be used with marky.sh because the pipeline in based on a snakemake structure. To work, you need to put your files in the sample_tmp folder as following:
* sample_tmp/sample_P1.fastq & sample_tmp/sample_P2.fastq # cleaned fastq files
* sample_tmp/sample_megahit/final.contigs.fa # megahit outputs
* sample_tmp/sample.bam # bowtie2 outputs
* sample_tmp/sample_proteins.faa # prodigal (or prokka) outputs
* sample_tmp/sample_counts.tsv # featureCounts outputs


## OUTPUTS
This software will produce a folder {sample}_outputs including:
* a file hgcA_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and the amino acid sequences. See <b>IMPORTANT NOTES</b> for data intepretation.
* a file hgcB_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp) and the amino acid sequences. 
* a file merA_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and the amino acid sequences.
* a file merB_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp) and the amino acid sequences. 
* a file rpoBb_final.txt (bacterial rpoB genes) and a file rpoBa_final.txt (archaeal rpoB genes) including coverage values (nb of read/bp).
* a file fastp.html and and a file fastp.html that are outputs for the fastp step.
* a file bowtie2.log that is output for the bowtie2 step.


## IMPORTANT NOTES
* You can modify the parameters for each step of the software in the file workflow/Snakefile and/or workflow/genesearch.sh.
* You can use all intermediate files produced in the sample_tmp folder to perform other type of analysis. See the method section below to see what they are.
* <b>True hgcA genes</b> include one of the amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* <b>True hgcB genes</b> include one of the amino acids motifs: CMECGA and CIEGCA
* To <b>find hgcB genes side-by-side with hgcA genes</b> in the same contig (so co-located in a microbial genome), look the 3nd number in their gene_id (corresponding to contigs id). Note that the 3rd number of the gene_id corresponds to the number of the genes on contigs. Ex k141_6000_1 and k141_6000_2 would be co-located genes.
* <b>merA and merB gene´s</b> homologs are detected in this pipeline using HMM profiles creating from the amazing database from Christakis et al. (2021). This allow only to screen roughly your metagenome and maybe detect and count merA and merB homologs. Manual inspection (following Christakis, Barkay and Boyd 2021 paper) is required to fully identify "true" merA and "true" merB genes but not described (yet) here. DO NOT USE these output table saying you detect merA and merB genes because it will be clearly wrong.
* To <b>normalize hgc coverage values</b>, sum the coverage values obtained from bacterial and archaeal rpoB genes.
* To <b>assign NCBI txid to the corresponding taxonomy</b>, you can use the R script below  or do a manual assignment with db/db_txid_2202220  if you have only few  hgcA gene homologs. If the database do not include the txid you found in your sample, check the identity here https://www.ncbi.nlm.nih.gov/taxonomy/
```
R
> db <- read.table("db/db_txid_220220.txt",h=T)
> a <- read.table("{sample}_outputs/{sample}_hgcA_final.txt", h=T)
> b <- merge(db, a, by="txid", all.x=F, all.y=T)
> write.table(b, file="{sample}_outputs/{sample}_hgcA_final2.txt", sep="\t", row.names=F)
> quit()
```


## METHODS
<div class="intro">
<p> The detection, counting and taxonomic identification of hgcAB genes was done with <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>marky-coco</b></a>. The metagenomes were trimmed and cleaned using <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>fastp</b></a> (Chen et al. 2018) with following parameters: -q 30 -l 25 --detect_adapter_for_pe --trim_poly_g --trim_poly_x. A de novo single assembly approach was applied using the assembler <a href="https://github.com/voutcn/megahit" target="_blank"><b>megahit</b></a> 1.1.2 (Li et al 2016) with default settings. The annotation of the contigs for prokaryotic protein-coding gene prediction was done with the software <a href="https://github.com/hyattpd/Prodigal" target="_blank"><b>prodigal</b></a> 2.6.3 (Hyatt et al 2010) (Hyatt et al., 2010). The DNA reads were mapped against the contigs with <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" target="_blank"><b>bowtie2</b></a> (Langdmead and Salzberg 2012), and the resulting .sam files were converted to .bam files using <a href="http://www.htslib.org/" target="_blank"><b>samtools</b></a> 1.9 (Li et al 2009). The .bam files and the prodigal output .gff file were used to estimate read counts by using <a href="https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html" target="_blank"><b>featureCounts</b></a>  (Liao et al 2014). In order to detect hgc homologs, HMM profiles derived from the <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE.db.v1</b></a> were applied to the amino acid FASTA file generated from each assembly with the function hmmsearch from <a href="http://hmmer.org/" target="_blank"><b>hmmer</b></a> 3.2.1 (Finn et al 2011). The reference package ‘hgcA’ from Hg-MATE.db.v1 was used for phylogenetic analysis of the HgcA amino acid sequences. Briefly, amino acid sequences from gene identified as hgcA gene homolog were (i) compiled in a FASTA file, (ii) aligned to Stockholm formatted alignment of hgcA sequences from the reference package with the function hmmalign from hmmer 3.2.1 (iii) placed onto the HgcA reference tree with the function pplacer and (iv) classified using the functions rppr and guppy_classify from the program <a href="https://matsen.fhcrc.org/pplacer/" target="_blank"><b>pplacer</b></a> (Masten et al. 2010).
