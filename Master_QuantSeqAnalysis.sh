#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/
#9/23/2019
#######################################################################################
#######################################################################################

#make poly A tail file to filter out 
#printf ">polyA\nAAAAAAAAAAAAA\n>polyT\nTTTTTTTTTTTTT\n" | gzip - >  polyA.fa.gz

#make sample files.txt for each sample
#from within fastq file folder: 
ls -R *fastq.gz > fastq.txt
#trim with:
#
cat fastq.txt | awk -F '_' '{print $1}' > samples.txt
#remove undetermined manually

#######################################################################################
#QC on individual sequencing runs before concatenating fastq files
#######################################################################################
mkdir fastqc_pretrim_run4

for sample in `cat samples.txt`
do
R1=${sample}
echo ${R1} "pre qc"

#fastqc on pre-trimmed file
fastqc ${R1}*.fastq.gz --outdir fastqc_pretrim_run4

echo ${R1} "post qc"

done

multiqc fastqc_pretrim_run4 --filename PreTrim_multiqc_report_run4.html --ignore-samples Undetermined* --interactive

#all fastq files looked good. concatenate biological replicates for hippocampus:
mkdir catfastq

cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA136H47_UMI_S47_R1_001.fastq.gz PA136H47_UMI_S52_R1_001.fastq.gz > catfastq/PA136H47_UMI_cat_R1_001.fastq.gz

cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA165H59_UMI_S59_R1_001.fastq.gz PA165H59_UMI_S60_R1_001.fastq.gz > catfastq/PA165H59_UMI_cat_R1_001.fastq.gz

cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA138H48_UMI_S48_R1_001.fastq.gz PA138H48_UMI_S53_R1_001.fastq.gz > catfastq/PA138H48_UMI_cat_R1_001.fastq.gz


cat cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/run1fastqs/PA145H51_UMI_S51_R1_001.fastq.gz /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/ook25hx4/Unaligned/Project_PARM_Tag107P_Moreno/PA145H51_UMI_S54_R1_001.fastq.gz > PA145H51_UMI_cat_R1_001.fastq.gz


cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA149H52_UMI_S52_R1_001.fastq.gz PA149H52_UMI_S55_R1_001.fastq.gz > catfastq/PA149H52_UMI_cat_R1_001.fastq.gz


cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA151H53_UMI_S53_R1_001.fastq.gz PA151H53_UMI_S56_R1_001.fastq.gz > catfastq/PA151H53_UMI_cat_R1_001.fastq.gz


cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA158H55_UMI_S55_R1_001.fastq.gz PA158H55_UMI_S57_R1_001.fastq.gz > catfastq/PA158H55_UMI_cat_R1_001.fastq.gz


cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA160H56_UMI_S56_R1_001.fastq.gz PA160H56_UMI_S58_R1_001.fastq.gz > catfastq/PA160H56_UMI_cat_R1_001.fastq.gz

cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/PA161H57_UMI_S57_R1_001.fastq.gz PA161H57_UMI_S59_R1_001.fastq.gz > catfastq/PA161H57_UMI_cat_R1_001.fastq.gz

#recheck fastqc pretim on cat files > all good
#manually replace run 1 HC files with new cat files 

#######################################################################################
#Trim and UMI Index
#######################################################################################
#intall packages if don't have them:
#install fastqc with brew or conda
#install multiqc with: pip install multiqc
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.
#from bbmap:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/ download bbduk.sh
#######################################################################################
#make folders > one set for each raw seq folder > 3
#path to raw fastq files: /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno

#/Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/brcnrfiug7/Unaligned/Project_PARM_Tag0107P_Moreno

#/Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/grkm4x3k32/Unaligned/Project_PARM_Tag0107P_Moreno

#cat all sample files togeter into one:

cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/0v99codvud/Unaligned/Project_PARM_Tag0107P_Moreno/samples.txt /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/brcnrfiug7/Unaligned/Project_PARM_Tag0107P_Moreno/samples.txt /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/grkm4x3k32/Unaligned/Project_PARM_Tag0107P_Moreno/samples.txt > samples.txt

#path to analysis: /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/
#path to output file: /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/
#scripts: /Users/aciernia/Desktop/MIATcellRNAseq/scripts/
#run code from in /Users/aciernia/Desktop/MIATcellRNAseq/

mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_pretrim
mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_posttrim
mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/trimmed_sequences
#run trim.sh srcipt:
#./trimHC.sh
#./trimFC.sh
#./trimCB.sh

#!/bin/bash
#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/grkm4x3k32/Unaligned/Project_PARM_Tag0107P_Moreno/samples.txt`
do
R1=${sample}
echo ${R1} "pre qc"

#fastqc on pre-trimmed file
fastqc /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/grkm4x3k32/Unaligned/Project_PARM_Tag0107P_Moreno/${R1}*.fastq.gz --outdir /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_pretrim

#index UMIs > add UMI sequence for each read to the read identifier
#The program umi2index pre-processes fastq files by extracting the UMI from a read and storing it as a read index. This step is performed right after demultiplexing to ensure that artificial UMI sequences do not influence subsequent analysis of endogenous read sequences. The program is invoked as follows:
#umi2index <fastq_gz_file_in> <fastq_gz_file_out>
echo ${R1} "umi2index"
./umi2index /Users/aciernia/Desktop/MIATcellRNAseq/rawsequencing/grkm4x3k32/Unaligned/Project_PARM_Tag0107P_Moreno/${R1}*.fastq.gz /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/${R1}.UMI.fastq.gz


### remove the adapter contamination, polyA read through, and low quality tails
##https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming
## quality-trim to Q10 using the Phred algorithm,
#set Java's memory usage to 24g with -Xmx24g
echo ${R1} "trimming"

/Users/aciernia/Desktop/programs/bbmap/bbduk.sh -Xmx24g in=/Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/${R1}.UMI.fastq.gz out=/Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/${R1}.trim.fastq.gz ref=/Users/aciernia/Desktop/scripts/polyA.fa.gz,/Users/aciernia/Desktop/scripts/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 threads=12


#fastqc on post-trimmed file
echo ${R1} "post qc"
fastqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/${R1}.trim.fastq.gz --outdir /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_posttrim


done

#######################################################################################
#collect all qc together with:

multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_pretrim/ --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/PreTrim_multiqc_report.html --ignore-samples Undetermined* --interactive

multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/fastqc_posttrim/ --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/PostTrim_multiqc_report.html --ignore-samples Undetermined* --interactive


#move to trimmed sequences folder
mv /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/*.trim.fastq.gz /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/trimmed_sequences
#######################################################################################
##align to mm10 with STAR
######################################################################################
#install STAR
#https://github.com/alexdobin/STAR
#git clone https://github.com/alexdobin/STAR.git
#brew install gcc
#Build STAR:
# note that the path to c++ executable has to be adjusted to its current version > 9
#cd source
#make STARforMacStatic CXX=/usr/local/Cellar/gcc/9.1.0/bin/g++-9
#run with ./STAR -h
#or add to path: export PATH=$PATH:/the_dir_with_STAR
#export PATH=$PATH:/Users/aciernia/Desktop/programs/STAR/source
#echo $PATH
#echo 'export PATH=/Users/aciernia/Desktop/programs/STAR/source:$PATH' >>~/.bash_profile
#source .bash_profile
# build indicies with mm10STARbuild.sh
#######################################################################################

#run trim.sh srcipt:
#./STAR_alignmm10.sh
mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out

#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIATcellRNAseq/samples.txt`
do
R1=${sample}

echo ${R1} "unzipping"

gunzip /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq.gz

echo ${R1} "mapping started"

#allows 10000 files to be open at once:
ulimit -n 10000

STAR --runThreadN 12 --genomeDir /Users/aciernia/Desktop/programs/STAR_libs/mm10/star_indices/ --readFilesIn  /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}

#rezip fastq
gzip /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq

echo ${R1} "mapping completed"
done


#######################################################################################
#mapping qc
#takes Log.final.out files and compiles into HTML

multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/*Log.final.out --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/STARAlignmentLogs.html --interactive


#######################################################################################
#UMI Filter after alignment
#######################################################################################
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.

#######################################################################################
#make folders

#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIATcellRNAseq/samples.txt`
do
R1=${sample}
echo ${R1} "index bam"

#Indexed bam files are necessary for many visualization and downstream analysis tools

samtools index /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bai

echo ${R1} "flagstat on unfiltered bam"
#flagstat on reads
samtools flagstat /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam > /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.preUMIcollapse.flagstat.txt

#collapse UMIs
#Program collapse_UMI_bam takes a bam file as input which is generated by aligning a fastq file to a reference genome. This fastq file is derived from <fastq_gz_file_out>, usually by adapter and quality trimming. The output bam file is used for downstream analysis. This version of collapse_UMI_bam has been tested with bam files generated by the star aligner
#collapse_UMI_bam <bam_file_in> <bam_file_out>

echo ${R1} "umi2index"
./collapse_UMI_bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.bam


echo ${R1} "coordinate sort and index"
samtools sort  /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.bam -o /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam

samtools index /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bai

#remove unmapped reads
#samtools view -h -F 4 -b ${R1}.collapseUMI.sort.bam > ${R1}.collapseUMI.mapped.bam

echo ${R1} "flagstat on filtered bam"
#flagstat on reads
samtools flagstat /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam > /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.flagstat.txt


done

#######################################################################################
#collect all qc flagstats together with:

multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/*.flagstat.txt --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/STARAlignment.html --interactive

#######################################################################################
#additional alignment QC with Qualimap
#######################################################################################
#intall packages if don't have them:
#install qualimap
#add the location of the Qualimap tool to our PATH variable
#echo 'export PATH=/Users/aciernia/Desktop/programs/qualimap_v2.2.1/:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#qualimap rnaseq --help
#######################################################################################
#make folders
mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_preUMI
mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_postUMI

#By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the unset DISPLAY command
unset DISPLAY



#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIATcellRNAseq/samples.txt`
do
R1=${sample}
echo ${R1} "qualimapQC RNAseq"

##-outdir: output directory for html report
#-a: Counting algorithm - uniquely-mapped-reads(default) or proportional (each multi-mapped read is weighted according to the number of mapped locations)
#-bam: path/to/bam/file(s)
#-p: Sequencing library protocol - strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
#-gtf: path/to/gtf/file - needs to match the genome build and GTF used in alignment
#--java-mem-size=: set Java memory

#Note that Qualimap must be run with the -outdir option as well as -outformat HTML (which is on by default). MultiQC uses files found within the raw_data_qualimapReport folder (as well as genome_results.txt).

qualimap rnaseq -bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_preUMI/${R1}.preUMI --java-mem-size=8G -p strand-specific-forward

qualimap rnaseq -bam /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_postUMI/${R1}.postUMI --java-mem-size=8G -p strand-specific-forward


done

#######################################################################################
#collect all qc together with:

multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_preUMI/*.preUMI/ --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/PreUMI_multiQualimapRNAseq.html --interactive
multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/qualimapQC_postUMI/*.postUMI/ --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/PostUMI_multiQualimapRNAseq.html --interactive

#######################################################################################
#count with subread feature counts
#######################################################################################
#install if don't have:
#http://bioinf.wehi.edu.au/featureCounts/
#download binary and add to path
#echo 'export PATH=/Users/aciernia/Desktop/programs/subread-1.6.5-MacOSX-x86_64/bin:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#######################################################################################

#counting

mkdir -p /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/counts


#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIATcellRNAseq/samples.txt`
do
R1=${sample}


#quant seq kit is FWD stranded

echo ${R1} "count started"
#0 (unstranded), 1 (stranded) and 2 (reversely stranded)
# Number of CPU threads
#set for gene_id > ensembl id mm10
#exon level counts

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/counts/${R1}.preUMI.counts.txt /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/counts/${R1}.collapsedUMI.counts.txt /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam

echo ${R1} "count completed"
done

#######################################################################################
#collect all qc together with:
multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/counts/*.preUMI.counts.txt.summary --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/Counts.preUMI --interactive
multiqc /Users/aciernia/Desktop/MIATcellRNAseq/analysis_files/counts/*.collapsedUMI.counts.txt.summary --filename /Users/aciernia/Desktop/MIATcellRNAseq/analysis_output/Counts.collapsedUMI --interactive
#######################################################################################
#######################################################################################
