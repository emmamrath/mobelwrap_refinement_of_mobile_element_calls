Mobelwrap refinement of mobile element calls
============================================

## Overview of Mobster+Mobelwrap pipeline to identify mobile element insertions

[Mobster](https://github.com/jyhehir/mobster) is a software that identifies mobile element insertions (MEIs) in a human genome BAM file of aligned sequencing reads. 
[Mobelwrap](https://github.com/emmamrath/mobelwrap_refinement_of_mobile_element_calls) is a software that carries out further processing of the Mobster-called MEIs, 
to genotype, filter, and further refines positions of MEIs, and derive an elongated sequence consisting of the MEI breakpoint. 
This processing can use multiple cores to speed up processing time. 
Subsequent processing then merges MEIs calls for multiple samples into a VCF file.  
  
There are 2 versions of Mobster. Please download both.  
* 1) Mobster-0.1.6 		from [https://sourceforge.net/projects/mobster/](https://sourceforge.net/projects/mobster/)  
* 2) MobileInsertions-0.2.2	from [https://github.com/jyhehir/mobster](https://github.com/jyhehir/mobster)  
  
Mobelwrap software is contained in this repository. 
Mobelwrap takes as input the Mobster results file and the bam file used by Mobster, to confirm true positives and filter false positives. 
Mobelwrap also derives the sequence spanning the mobile element insertion point.
This sequence consists of reference sequence and appended mobile element sequence, as observed in the bam, 
and is a concatenation of multiple reads in the bam.
Mobelwrap records this extended sequence in the Mobelwrap output. 
Mobelwrap output is the same tab-delimited format as Mobster output, with extra columns.
A subsequent program converts the Mobelwrap output to VCF file.  
  
This Mobelwrap pipeline uses the results from Mobster version MobileInsertions-0.2.2. 
As of May 2017, data files from the Mobster-0.1.6 package were needed in order to run MobileInsertions-0.2.2,
as the data files in the MobileInsertions-0.2.2 package were corrupted.  
  
In addition to reading this overview, please read the Mobster README.md file in the Mobster package 
and the Mobster FAQ at [https://jyhehir.github.io/mobster/faq.html](https://jyhehir.github.io/mobster/faq.html).  
  
The citation for Mobster is:  
Genome Biol. 2014;15(10):488. Mobster: accurate detection of mobile element insertions in next generation sequencing data.  
Thung DT, de Ligt J, Vissers LE, Steehouwer M, Kroon M, de Vries P, Slagboom EP, Ye K, Veltman JA, Hehir-Kwa JY.  
[https://www.ncbi.nlm.nih.gov/pubmed/25348035](https://www.ncbi.nlm.nih.gov/pubmed/25348035)  

Mobelwrap results were used in the following paper which can be considered the citation for the Mobelwrap software:  
[The Medical Genome Reference Bank: Whole genomes and phenotype of 2,570 healthy elderly](https://www.biorxiv.org/content/10.1101/473348v1)  
  
Mobster resources as of May 2017:  
* The Mobster FAQ is at:			[https://jyhehir.github.io/mobster/faq.html](https://jyhehir.github.io/mobster/faq.html)  
* The Mobster source code is at:		[https://github.com/jyhehir/mobster](https://github.com/jyhehir/mobster)  
* The Mobster contact is:			[Jayne.Hehir-Kwa@radboudumc.nl](mailto:Jayne.Hehir-Kwa@radboudumc.nl)  
  
Mobster is Java programs in a jar file. 
Mobster reads a BAM file on input, and produces an output tab-delimited text file of mobile element insertions (MEI) it has identified. 
The Mobster jar file also uses the Picard and MOSAIK programs that are supplied in the Mobster packages.  
  
Terminology: Mobster and Mobelwrap report the border5 and border3 genomic position of the MEI.  
* Border5 is the position in the reference genome whose nucleotide is the same as the ref.seq. and to the left of it is the MEI sequence.  
* Border3 is the position in the reference genome whose nucleotide is the same as the ref.seq. and to the right of it is the MEI sequence.  
  
## Input parameters to the Mobelwrap program

* -i <in_Mobster> (required) Input file containing Mobster predictions  
* -b <in_bam_reads> (required) Input BAM file containing mapped reads. Needs to be indexed. BAM reads in the region of each Mobster mobile-element-insertion-prediction will be retrieved from this BAM file to determine the exact insertions points of the Mobster mobile-element-insertion-prediction  
	parser.add_argument('-l', action="store", dest="read_length_of_in_bam_reads", required=True, help='The general read length of the reads in the input BAM file')
	parser.add_argument('-r', action="store", dest="in_reference_genome", required=True, help='Reference genome. A fasta file indexed by samtools faidx. BAM reads will be compared to this reference to determine the exact insertion points for Mo
bster mobile-element-insertion-predictions')
	parser.add_argument('-e', action="store", dest="existing_mobile_element_regions", required=True, help='Mobster\'s repmask file \'hg19_alul1svaerv.txt\' consisting of 1 line header and 1 line per mobile element region already existing in the
 reference genome.')
	parser.add_argument('-o', action="store", dest="out_Mobster_file", required=True, help='Output file will contain refined Mobster mobile-element-insertion-predictions, and where possible, extra information appended to the record to precisely
 identify the insertion points. (And if this refinement program does not find extra information for Mobster MEIs in existing ME region, then the MEI is dropped.)')
	parser.add_argument('-s', action="store", dest="restart_Mobster_file", required=False, help='If this file is present, then keep the MEI calls from this file and continue processing input Mobster calls from the next input Mobster MEI after t
he last MEI call seen in this file.')
	parser.add_argument('-q', action="store", dest="quality_filter", required=False, help='If present, then filter MEI calls, keeping calls having this or higher quality score.')
	parser.add_argument('-d', action="store", dest="discard_Mobster_file", required=False, help='If this output file is present, then record all discarded Mobster MEI calls to this file.')
	parser.add_argument('-c', action="store", dest="num_cores", required=False, help='Number of concurrent processors for multiprocessor processing')
	parser.add_argument('-k', action="store", dest="keep_chr_prefix", required=False, help='Mobster puts chr in front of chromosome and this program removes it. If this flag is present then this program will not remove the chr prefix.')
	parser.add_argument('-p', action="store", dest="python_path_2", required=False, help='The path where the python program called by this program is found. Alternatively, set PYTHONPATH2 environmental variable to this value.')
	args = parser.parse_args()


## How to prepare a Mobster+Mobelwrap pipeline
  
```
cd /your/software/git_directory  
git clone https://emmrat@kccg.garvan.org.au/bitbucket/scm/gcr/mobster_for_mobile_elements.git  
```
  
The code from git repo and the run directory need to be different to each other
because the install requires running script install.sh 
which will not have permission to create directories it needs in the cloned directory.  
```
cp -r /your/software/git_directory /your/software/run_directory  
```

Install into your directory for which you have write permission.  
```
cd /your/software/run_directory  
cd mobster_for_mobile_elements  
unzip Mobster-0.1.6.zip  
ln -s Mobster-0.1.6 Mobster  
tar -xvf MOSAIK-2.2.3-Linux-x64.tar  
ln -s MOSAIK-2.2.3-Linux-x64 MOSAIK  
tar -xvf MOSAIK-2.2.3-source.tar  
cp MOSAIK-2.2.3-source/networkFile/2.1.26.pe.100.0065.ann Mobster/MOSAIK  
cp MOSAIK-2.2.3-source/networkFile/2.1.26.se.100.005.ann Mobster/MOSAIK  
export PATH=/your/mobster/installation/path/mobster_for_mobile_elements/MOSAIK:$PATH  
```

The provided precompiled MobileInsertions-0.2.2.jar was compiled on a Linux-x64 machine. 
If the precompiled MobileInsertions-0.2.2.jar does not work on your system,
then carry out the following additional steps to compile it.

```
tar -xvf mobster-master.tar &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; or &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; unzip mobster-master.zip  
cd mobster-master  
./install.sh  
cp target/MobileInsertions-0.2.2.jar ..  
```

Copy and edit jars/Mobster.properties to specify input file and output file prefix for the Mobster jobs.  
Also make sure that paths to other programs in the Mobster.properties file are specified correctly,  
for picard-1.73/CollectInsertSizeMetrics.jar, MosaikBuild, MosaikAligner, files in their command lines, REPEATMASK_FILE, and uncomment TEMP_DIR.  
The lines are already set correctly for MAPPING_TOOL=bwa and MOBIOME_MAPPING_TOOL=mosaik  

```
cd /your/software/directory  
cp mobster_for_mobile_elements/Mobster/jars/Mobster.properties /your/work/directory/your_sample_Mobster.properties  
nano your_sample_Mobster.properties  
	SAMPLENAME=test_sample									==> 		SAMPLENAME=Your_Sample_Name  
	REPEATMASK_FILE=../repmask/hg19_alul1svaerv.txt						==>		REPEATMASK_FILE=/your/software/directory/mobster_for_mobile_elements/extra_data_for_Mobster/hg19_alul1svaerv_empty.txt  
	IN_FILE=../test_files/chr12_wgs_3000ins_100excl_10x_bwadef.bam				==>		IN_FILE=/your/bam/directory/your_sample_noH.bam		(your bam file with hard clips removed as explained below)  
	OUT_FILE=../test_output/chr12_wgs_3000ins_100excl_usingproperties			==>		OUT_FILE=/your/work/directory/your_sample_MobsterOUTPUT  
	#TEMP_DIR=/home/djie/data/public_trio/tmp/						==>		TEMP_DIR=/your/tmp  
	MAPPING_TOOL=bwa									==>		MAPPING_TOOL=bwa		(make sure you have the correct value according to which tool was previously used to map reads in your bam file)  
	PICARD_COLLECT_INSERT_SIZE_METRICS_JAR=../picard-1.73/CollectInsertSizeMetrics.jar	==>		PICARD_COLLECT_INSERT_SIZE_METRICS_JAR=/your/software/directory/mobster_for_mobile_elements/Mobster/picard-1.73/CollectInsertSizeMetrics.jar  
	MOBIOME_MAPPING_CMD=MosaikBuild -q (FASTQ) -st illumina -out (DAT_FILE) -quiet && MosaikAligner -in (DAT_FILE) -out (OUT_FILE) -ia ../mobiome/54_mobiles_inclHERVK.dat -hs 9 -mmp 0.1 -act 20 -j ../mobiome/54_mobiles_inclHERVK_hs9 -p 2 -annpe ../MOSAIK/2.1.26.pe.100.0065.ann -annse ../MOSAIK/2.1.26.se.100.005.ann -quiet  
	==>  
	MOBIOME_MAPPING_CMD=/your/software/directory/mobster_for_mobile_elements/MOSAIK/MosaikBuild -q (FASTQ) -st illumina -out (DAT_FILE) -quiet && /your/software/directory/mobster_for_mobile_elements/MOSAIK/MosaikAligner -in (DAT_FILE) -out (OUT_FILE) -ia /your/software/directory/mobster_for_mobile_elements/Mobster/mobiome/54_mobiles_inclHERVK.dat -hs 9 -mmp 0.1 -act 20 -j /your/software/directory/mobster_for_mobile_elements/Mobster/mobiome/54_mobiles_inclHERVK_hs9 -p 2 -annpe /your/software/directory/mobster_for_mobile_elements/Mobster/MOSAIK/2.1.26.pe.100.0065.ann -annse /your/software/directory/mobster_for_mobile_elements/Mobster/MOSAIK/2.1.26.se.100.005.ann -quiet
```

The following 3 jobs in the jobs_for_running_Mobster_Refined_pipeline will run the pipeline for multiple samples.  
You would need to put the list of samples (bam_file and sample_id) in list_of_sample_files.txt 
and modify the 5 scripts PBS_run_Mobster2_for_multiple_samples.sh, PBS_run_Mobster2_for_one_sample.sh, PBS_Mobelwrap_run_for_multiple_samples.sh, PBS_Mobelwrap_run_for_one_sample.sh, merge_VCF_files_then_annotate.sh  
```
./PBS_run_Mobster2_for_multiple_samples.sh  
./PBS_Mobelwrap_run_for_multiple_samples.sh  
./merge_VCF_files_then_annotate.sh  
```

## How to run a Mobster-Mobelwrap pipeline, from calling mobile-element-insertions (MEIs) to final VCF file containing MEIs for all samples

Remove reads from the BAM file that will make Mobster crash.   
Remove hard clips because they are the second read to match its read pair, it's soft clip read is the first read to match their read pair, and Mobster will crash thinking that the hard clip has no matching read pair because the read pair will have already been matched to the soft clip.  
Remove supplementary reads because Mobster AnchorClusterer will crash with them present.  
For a 75G BAM file, this produces a 352G BAM file because a lot of it is not actually recompressed.  
```
samtools view -h -F 2048 /your/work/directory/your_sample.bam | awk -F'\t' '( ($6 !~ /H/) || ($1 ~ /^@/) )' | samtools view -Shu - > /your/work/directory/your_sample_noH.bam
```

Another way to remove reads from the BAM file that will make Mobster crash is the following.  
Write out uncompressed sam and then compress to bam, because samtools compression is slower than awk output buffering and thus piping awk output to samtools compression will result in most of file not being compressed.  
For a 75G BAM file, it first produces a 477G SAM file, then produces a 77G BAM file.  
```
samtools view -h -F 2048 /your/work/directory/your_sample.bam | awk -F'\t' '( ($6 !~ /H/) || ($1 ~ /^@/) )' > /your/work/directory/your_sample_noH.sam  
samtools view -bS /your/work/directory/your_sample_noH.sam > /your/work/directory/your_sample_noH.bam  
```

Then run Mobster Mobster.jar which uses Mobster.properties as input.  
This runs Picard's CollectInsertSizeMetrics.jar, Mobster's PotentialMEIReadFinder, Mosaiks's MosaikBuild, Mosaik's MosaikAligner, Mobster's RefAndMEPairFinder, and Mobster's AnchorClusterer  
```
java -Xmx64g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/your/tmp -jar /your/software/directory/mobster_for_mobile_elements/Mobster/jars/Mobster.jar -properties /your/work/directory/your_sample_Mobster.properties
```

Alternatively, run the following Mobster jobs, which does not first run Picard CollectInsertSizeMetrics.jar  
```
java -Xmx8g -Djava.io.tmpdir=/your/tmp -jar /your/software/directory/mobster_for_mobile_elements/Mobster/jars/PotentialMEIReadFinder.jar -in  /your/work/directory/your_sample_noH.bam -maxclip 7 -min_avg_qual 10 -minclip 35 -out  /your/work/directory/your_sample_noH_Mobster_OUT -split -tool bwa  
/your/software/directory/mobster_for_mobile_elements/MOSAIK/MosaikBuild -q /your/work/directory/your_sample_noH_Mobster_OUT_potential.fq -st illumina -out /your/work/directory/your_sample_noH_MOSAIK_DAT_FILE -quiet  
/your/software/directory/mobster_for_mobile_elements/MOSAIK/MosaikAligner -in /your/work/directory/your_sample_noH_MOSAIK_DAT_FILE -out /your/work/directory/your_sample_noH_MosaikAligner_OUT -ia /your/software/directory/mobster_for_mobile_elements/Mobster/mobiome/54_mobiles_inclHERVK.dat -hs 9 -mmp 0.1 -act 20 -j /your/software/directory/mobster_for_mobile_elements/Mobster/mobiome/54_mobiles_inclHERVK_hs9 -p 2 -annpe /your/software/directory/mobster_for_mobile_elements/MOSAIK/2.1.26.pe.100.0065.ann -annse /your/software/directory/mobster_for_mobile_elements/MOSAIK/2.1.26.se.100.005.ann -quiet  
java -Xmx8g -Djava.io.tmpdir=/your/tmp -jar /your/software/directory/mobster_for_mobile_elements/Mobster/jars/RefAndMEPairFinder.jar -max_memory 10000000 -multiple /your/work/directory/your_sample_noH_MosaikAligner_OUT.bam -out /your/work/directory/your_sample_noH_Mobster_OUT_RefAndMEPairFinder_OUT_ANCHORS -p -potential /your/work/directory/your_sample_noH_Mobster_OUT_potential.bam -samplename mySample -single /your/work/directory/your_sample_noH_MosaikAligner_OUT.bam -tool mosaik  
java -Xmx8g -Djava.io.tmpdir=/your/tmp -jar /your/software/directory/mobster_for_mobile_elements/Mobster/jars/AnchorClusterer.jar -in /your/work/directory/your_sample_noH_Mobster_OUT_RefAndMEPairFinder_OUT_ANCHORS_anchors.bam -insplit /your/work/directory/your_sample_noH_Mobster_OUT_RefAndMEPairFinder_OUT_ANCHORS_splitanchors.bam -max_memory 10000000 -maxclust 384 -maxdist 701 -mfl 444 -mintotal 5 -out /your/work/directory/your_sample_noH_Mobster_OUT_AnchorClusterer -overlap 50 -repmask /your/software/directory/mobster_for_mobile_elements/Mobster/repmask/hg19_alul1svaerv.txt -rpc 1 -sample mySample -sd 117 -search 200 -splithits 2  
```

### An example of the contents of a script to run Mobster using the one Mobster properties to run all steps including an initial Picard CollectInsertSizeMetrics step:

```
 #!/bin/bash  
 #PBS -l ncpus=1  
 #PBS -l jobfs=2GB  
 #PBS -l mem=64GB  
 #PBS -l walltime=38:00:00  
 #PBS -P wq2  
 #PBS -q normal  
 module load java/jdk1.7.0_25  
 module load samtools/1.3.1  
 module load R/3.2.2  
 export PATH=/your/software/directory/mobster_for_mobile_elements/MOSAIK:$PATH  

 samtools view -h -F 2048 /your/samples/your_sample.bam | awk -F'\t' '(($6 !~ /H/) || ($1 ~ /^@/))' | samtools view -Shu - > /your/work/directory/your_sample_noH.bam  

 java -Xmx64g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/your/tmp -jar /your/software/directory/mobster_for_mobile_elements/Mobster/jars/Mobster.jar -properties /your/work/directory/your_sample.properties  

 cp /your/work/directory/your_sample_MobsterOUTPUT_predictions.txt /your/work/directory/your_sample_predictions.txt  
 # rm /your/work/directory/your_sample_MobsterOUTPUT*  
 # rm /your/work/directory/your_sample_noH.bam  
 # rm /your/work/directory/your_sample.properties  
```

### An example of the contents of a script to run the alternative sequence of Mobster jobs, that removes intermediate files to save disk space:

```
 #!/bin/bash
 #PBS -l ncpus=1
 #PBS -l jobfs=2GB
 #PBS -l mem=20GB
 #PBS -P wq2
 #PBS -q normal
 module load java/jdk1.8.0_60
 module load samtools/1.3.1
 module load R/2.15.2-gcc
 export PATH=/short/wq2/emr913/mobile_elements_2016dec/Mobster_for_running/MOSAIK:$PATH

 samtools view -h -F 2048 /g/data3/wq2/MGRB/alignments/phase1_samples/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam | awk -F'\t' '(($6 !~ /H/) || ($1 ~ /^@/))' | samtools view -Shu - > /g/data3/wq2/MGRB_mobile_elements_2016dec/run_Mobster_4steps444_MGRBp1_2016dec/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.noHardClips.bam

 java -Xmx8g -Djava.io.tmpdir=/tmp -jar jars/PotentialMEIReadFinder.jar -in your_sample_noH.bam -maxclip 7 -min_avg_qual 10 -minclip 35 -out your_sample_noH_OUT -split -tool bwa  
MOSAIK/MosaikBuild -q your_sample_noH_OUT_potential.fq -st illumina -out your_sample_noH_OUT_potential_fq_MOSAIK_DAT_FILE -quiet  
 MOSAIK/MosaikAligner -in your_sample_noH_OUT_potential_fq_MOSAIK_DAT_FILE -out your_sample_noH_OUT_potential_fq_MosaikAligner_OUT_FILE -ia Mobster/mobiome/54_mobiles_inclHERVK.dat -hs 9 -mmp 0.1 -act 20 -j Mobster/mobiome/54_mobiles_inclHERVK_hs9 -p 2 -annpe Mobster/MOSAIK/2.1.26.pe.100.0065.ann -annse Mobster/MOSAIK/2.1.26.se.100.005.ann -quiet  
 java -Xmx8g -Djava.io.tmpdir=/tmp -jar Mobster/jars/RefAndMEPairFinder.jar -max_memory 10000000 -multiple your_sample_noH_OUT_potential_fq_MosaikAligner_OUT_FILE.bam -out your_sample_noH_OUT_RefAndMEPairFinder_OUT_ANCHORS -p -potential your_sample_noH_OUT_potential.bam -sample_name MGRB -single your_sample_noH_OUT_potential_fq_MosaikAligner_OUT_FILE.bam -tool mosaik
java -Xmx8g -Djava.io.tmpdir=/tmp -jar jars/AnchorClusterer.jar -in your_sample_noH_OUT_RefAndMEPairFinder_OUT_ANCHORS_anchors.bam -insplit your_sample_noH_OUT_RefAndMEPairFinder_OUT_ANCHORS_splitanchors.bam -max_memory 10000000 -maxclust 384 -maxdist 701 -mfl 444 -mintotal 5 -out your_sample_noH_OUT_AnchorClusterer -overlap 50 -repmask Mobster/repmask/hg19_alul1svaerv.txt -rpc 1 -sample MGRB -sd 117 -search 200 -splithits 2  
```
### How to run Mobelwrap on Mobster output

We now have a text file of mobile element insertions (MEIs) for each sample.  
The file may contain duplicate MEIs for the same sample.   
Any duplicates will be removed by GATK tools CombineVariants further down when it is used to merge multiple VCF files into one VCF file.  

For each sample, run a python program to further refine the Mobster MEI calls.  
The program is convert_Mobster_predictions_into_Refined_predictions and it calls identify_split_read_clusters.py on the PYTHONPATH2 path.  
The flag &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -l 150 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; is specifying the read-length of the reads in the bam file.  
The flag &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -e hg19_alul1svaerv.txt &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; is specifying the file containing existing ME regions in the reference genome. This is the 'repmask' file is from Mobster.  
In the previous step in this Refined Mobster pipeline, the repmask file used to run Mobster was empty (contained only a header) so that Mobster calls MEIs even if they are within 90 base-pairs of an existing ME.  
The refinement program will drop an MEI called by Mobster if no split-read evidence is found for its existance and it falls within 90 base-pairs of an existing ME.  
The input BAM file to this inhouse refinement program can be the original BAM file or it can be the BAM file having hard clipped reads removed (that was necesssary for the Mobster step).  
In theory, the inhouse refinement program should give the same result regardless of which of the 2 BAM files is used.  
In practise, the resulting read counts and precise pin-pointing of MEI border5 and border3 differ slightly and insignificantly depending on which BAM file is used.
The -l parameter is the length of the bam reads.  
The optional -q parameter will cause the program to drop MEI calls having a QUAL quality_score less than the parameter's value.  
If the option -d is present, dropped MEI calls will be written to that file.  
Provide a fastq of human reference genome (same one used for bwa mem), make sure it is indexed (it was already indexed for running bwa mem).  

```
export PYTHONPATH2=/path/to/python/programs  
python convert_Mobster_predictions_into_Refined_predictions.py \  
	-i Mobster_calls_for_sample1 \  
	-b sample1.bam \  
	-l 150 \  
	-r indexed_reference_genome_fasta \  
	-e hg19_alul1svaerv.txt \  
	-o refined_calls_for_sample1 \  
	-q 13  
python convert_Mobster_predictions_into_Refined_predictions.py \  
	-i Mobster_calls_for_sample2 \  
	-b sample2.bam \  
	-l 150 \  
	-r indexed_reference_genome_fasta \  
	-e hg19_alul1svaerv.txt \  
	-o refined_calls_for_sample2 \  
	-q 13  
python convert_Mobster_predictions_into_Refined_predictions.py \  
	-i Mobster_calls_for_sample3 \  
	-b sample3.bam \  
	-l 150 \  
	-r indexed_reference_genome_fasta \  
	-e hg19_alul1svaerv.txt \  
	-o refined_calls_for_sample3 \  
	-q 13  
```

If the refinement program convert_Mobster_predictions_into_Refined_predictions.py crashes partway through running (due to exceeding time limits),
then it can be restarted to continue where it left off:

```
export PYTHONPATH2=/path/to/python/programs  
python convert_Mobster_predictions_into_Refined_predictions.py \  
	-i Mobster_calls_for_sample1 \  
	-b sample1.bam \  
	-l 150 \  
	-r indexed_reference_genome_fasta \  
	-e hg19_alul1svaerv.txt \  
	-f refined_calls_for_sample1_output_from_crashed_run \  
	-o refined_calls_for_sample1_continuation_of_output \  
	-q 13  
```

Convert each sample's MEI file to a VCF file.  
In -s parameter, provide the sample_id to appear in the VCF header  
Provide a fastq of human reference genome (same one used for bwa mem), make sure it is indexed (it was already indexed for running bwa mem).  
The input MEI calls is either the output file from Mobster, or the output file from the Refined program.  

```
python create_VCF_file_from_Mobster_predictions_for_one_sample.py \  
	-i refined_calls_for_sample1 \  
	-r indexed_reference_genome_fasta \  
	-s sample1_ID \  
	-o MEI_calls_for_sample1.vcf  
python create_VCF_file_from_Mobster_predictions_for_one_sample.py \  
	-i refined_calls_for_sample2 \  
	-r indexed_reference_genome_fasta \  
	-s sample2_ID \  
	-o MEI_calls_for_sample2.vcf  
python create_VCF_file_from_Mobster_predictions_for_one_sample.py \  
	-i refined_calls_for_sample3 \  
	-r indexed_reference_genome_fasta \  
	-s sample3_ID \  
	-o MEI_calls_for_sample3.vcf  
```

The next step is to merge the multiple VCF files (one per sample) into 1 VCF file (contains all the samples)  
Run GATK tools CombineVariants to merge the multiple VCF files into one VCF file.  
It requires the reference sequence fasta file (indexed by samtools faidx?).  
It requires the VCF files to be in the same chromosome order as the reference sequence fasta file   
and for all REF values to be present in the VCF (so no POS set to zero or less than zero).  
The output from Mobster is not necessarily in the chromosome order of the ref.seq. fasta file  
and is even not guaranteed to be in ascending POS order within a chromosome/contig.  
So run jobs to do those sorting things first.  

```
python prepare_VCF_before_running_GATK.py -i MEI_calls_for_sample1.vcf -r /reference_genomes/hs37d5x/hs37d5x.fa -o MEI_calls_for_sample1_fixPosAndRef.vcf  
grep '^#' MEI_calls_for_sample1_fixPosAndRef.vcf > MEI_calls_for_sample1_fixPosAndRef_hdr.vcf  
grep -v '^#' MEI_calls_for_sample1_fixPosAndRef.vcf > MEI_calls_for_sample1_fixPosAndRef_nohdr.vcf  
sort -k1 -nk2 MEI_calls_for_sample1_fixPosAndRef_nohdr.vcf > MEI_calls_for_sample1_fixPosAndRef_nohdr_sorted.vcf  
cat MEI_calls_for_sample1_fixPosAndRef_hdr.vcf MEI_calls_for_sample1_fixPosAndRef_nohdr_sorted.vcf > MEI_calls_for_sample1_fixPosAndRef_sorted.vcf  
python sort_VCF_by_order_found_in_input_list.py -i MEI_calls_for_sample1_fixPosAndRef_sorted.vcf -l hs37d5x_sort_order.txt -o MEI_calls_for_sample1_fixPosAndRef_sortedByRefSeqOrder.vcf  
  
python prepare_VCF_before_running_GATK.py -i MEI_calls_for_sample2.vcf -r /reference_genomes/hs37d5x/hs37d5x.fa -o MEI_calls_for_sample2_fixPosAndRef.vcf  
grep '^#' MEI_calls_for_sample2_fixPosAndRef.vcf > MEI_calls_for_sample2_fixPosAndRef_hdr.vcf  
grep -v '^#' MEI_calls_for_sample2_fixPosAndRef.vcf > MEI_calls_for_sample2_fixPosAndRef_nohdr.vcf  
sort -k1 -nk2 MEI_calls_for_sample2_fixPosAndRef_nohdr.vcf > MEI_calls_for_sample2_fixPosAndRef_nohdr_sorted.vcf  
cat MEI_calls_for_sample2_fixPosAndRef_hdr.vcf MEI_calls_for_sample2_fixPosAndRef_nohdr_sorted.vcf > MEI_calls_for_sample2_fixPosAndRef_sorted.vcf  
python sort_VCF_by_order_found_in_input_list.py -i MEI_calls_for_sample2_fixPosAndRef_sorted.vcf -l hs37d5x_sort_order.txt -o MEI_calls_for_sample2_fixPosAndRef_sortedByRefSeqOrder.vcf  
  
python prepare_VCF_before_running_GATK.py -i MEI_calls_for_sample3.vcf -r /reference_genomes/hs37d5x/hs37d5x.fa -o MEI_calls_for_sample3_fixPosAndRef.vcf  
grep '^#' MEI_calls_for_sample3_fixPosAndRef.vcf > MEI_calls_for_sample3_fixPosAndRef_hdr.vcf  
grep -v '^#' MEI_calls_for_sample3_fixPosAndRef.vcf > MEI_calls_for_sample3_fixPosAndRef_nohdr.vcf  
sort -k1 -nk2 MEI_calls_for_sample3_fixPosAndRef_nohdr.vcf > MEI_calls_for_sample3_fixPosAndRef_nohdr_sorted.vcf  
cat MEI_calls_for_sample3_fixPosAndRef_hdr.vcf MEI_calls_for_sample3_fixPosAndRef_nohdr_sorted.vcf > MEI_calls_for_sample3_fixPosAndRef_sorted.vcf  
python sort_VCF_by_order_found_in_input_list.py -i MEI_calls_for_sample3_fixPosAndRef_sorted.vcf -l hs37d5x_sort_order.txt -o MEI_calls_for_sample3_fixPosAndRef_sortedByRefSeqOrder.vcf  
```

Merging the MEI calls of all samples into a big VCF file.
```
java -jar /GATK_tools/GenomeAnalysisTK.jar \  
	-T CombineVariants \  
	-R Reference_Sequence.fa \  
	--variant MEI_calls_for_sample1_fixPos_sortedByRefSeqOrder.vcf \  
	--variant MEI_calls_for_sample2_fixPos_sortedByRefSeqOrder.vcf \  
	--variant MEI_calls_for_sample3_fixPos_sortedByRefSeqOrder.vcf \  
	-o all_samples_MEI_calls.vcf  

python prepare_VCF_after_running_GATK.py -i all_samples_MEI_calls.vcf -o all_samples_MEI_calls_afterGATK.vcf
```

MEI calls that are close to each other (within 100 base-pairs) are likely the same MEI, so merge them.  
Variants from all samples are now in one VCF file.  
Where 2 samples had an exact match for a mobile element insertion call, they appear in the one VCF record.  
Where 2 samples had a similar match for a mobile element insertion call, they remain in 2 separate records, even though they are probably the same mobile element insertion.  
So now run a program to fuzzy match such variants.  
If 2 different mobile element insertion calls are within a certain number of base pairs of each other, then this program will merge them into the one VCF record.  
The number of base pairs is specfied on input to the program.  
After looking at examples, 100 seems to be a good window. Mobster mobile element insertion calls within 100 bp of each other seem to be the same mobile element insertion.  

```
cat all_samples_MEI_calls_afterGATK.vcf | python merge_a_multi_sample_VCF_for_similar_position_variants.py 100 indexed_reference_genome_fasta > all_samples_MEI_calls_merged100bp.vcf
```

Optionally, some annotations can be added to each VCF MEI call, specifying whether MEI falls in a gene exon, intron, regulatory region, etc.  
The python program to do this annotation is in the annotation repository, and is not in the Mobster repository.  
```
cat all_samples_MEI_calls_merged100bp.vcf \
	| python annotate_VCF_with_an_extra_database.py \
	-t ME-UCSC-Gene \
	-d hgtables_2017mar27_Genes_for_MEI_annotation.txt \
	> all_samples_MEI_calls_merged100bp_annotated.vcf
```

Here is the file containing the output; a VCF file of mobile element insertion calls for all samples.
```
less -S all_samples_MEI_calls_merged100bp_annotated.vcf
```

## What to do when Mobster takes too long to run or requires too much memory.

What to do when Mobster takes too long to run or requires too much memory.  
Split up BAM file by chromosome 
making sure that reads from this chromosome having discordant read pairs from other chromosomes are included in this chromosomes BAM file,
then run the Mobster pipeline, then merge predictions with merge_Mobster_predictions.py.  
When spliting by chromosome, make sure the chromosome's BAM file has discordant read pairs mapped elsewhere and having discordant read mapping to this chromosome of interest.
Here is an example of spliting BAM file into 2 groups by chromosome.  

```
sample_id='sample1000_AAAAB'
#bam_file='AAAAB.sorted.dupmarked.bam'
output_directory='/nvme/mobile_elements_2016dec/debug_MGRB_picard_header'
bam_file_noHardClips=${output_directory}/'MGRBp1_sample1000_AAAAB.sorted.dupmarked.bam_noHardClips.bam'
```
just the header
```
echo "EXECUTE just the header"
cat "${bam_file_noHardClips}" \  
	| awk -F'\t' '($1 ~ /^@/)' \  
	> "${output_directory}"/MGRBp1_"${sample_id}"_hdr.sam  
```
just reads mapped to chr 1:11
```
echo "EXECUTE just reads mapped to chr 1:11"  
samtools view "${bam_file_noHardClips}" \  
	| awk -F'\t' '(($3 == "1") || ($3 == "2") || ($3 == "3") || ($3 == "4") || ($3 == "5") || ($3 == "6") || ($3 == "7") || ($3 == "8") || ($3 == "9") || ($3 == "10") || ($3 == "11"))' \  
	> "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam  
```
just reads mapped to all chromosomes except chr 1:11
```
echo "EXECUTE just reads mapped to all chromosomes except chr 1:11"  
samtools view "${bam_file_noHardClips}" \  
	| awk -F'\t' '(($1 !~ /^@/) && ($3 != "1") && ($3 != "2") && ($3 != "3") && ($3 != "4") && ($3 != "5") && ($3 != "6") && ($3 != "7") && ($3 != "8") && ($3 != "9") && ($3 != "10") && ($3 != "11"))' \
	> "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam  
```
the IDs of reads mapped to chr 1:11
```
echo "EXECUTE the IDs of reads mapped to chr 1:11"  
cat "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam \  
	| cut -f1 > "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_IDs.txt  
```
the IDs of reads mapped to all chromosomes except chr 1:11
```
echo "EXECUTE the IDs of reads mapped to all chromosomes except chr 1:11"  
cat "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam \  
	| cut -f1 > "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_IDs.txt  
```
the non-chr1:11 read-pairs of reads mapped to chr 1:11
```
echo "EXECUTE the non-chr1:11 read-pairs of reads mapped to chr 1:11"
```
while read in_ID_line; do  
	grep $in_ID_line "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam  
done < "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_IDs.txt > "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_read_pairs_for_chr1o11.sam
```
(This grep crashes - too much data for grep)
```
cat "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam \  
	| grep -f "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_IDs.txt \  
	> "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_read_pairs_for_chr1o11.sam
```
the chr1:11 read-pairs of reads not mapped to chr 1:11
```
echo "EXECUTE the chr1:11 read-pairs of reads not mapped to chr 1:11"
while read in_ID_line; do
	grep $in_ID_line "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam  
done < "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_IDs.txt > "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_read_pairs_for_NOTchr1o11.sam
``
(This grep crashes - too much data for grep)
```
cat "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam \  
	| grep -f "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_IDs.txt > "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_read_pairs_for_NOTchr1o11.sam
```

create indexed sam file of reads mapped to chr1:11 including all their read-pairs
```
echo "EXECUTE create indexed sam file of reads mapped to chr1:11 including all their read-pairs"
cat "${output_directory}"/MGRBp1_"${sample_id}"_hdr.sam \
	"${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam \
	"${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_read_pairs_for_chr1o11.sam \
	> "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11.sam
samtools view -bS "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11.sam \
	> "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11.bam
```

create indexed sam file of reads not mapped to chr1:11 including all their read-pairs
```
echo "EXECUTE create indexed sam file of reads not mapped to chr1:11 including all their read-pairs"  
cat "${output_directory}"/MGRBp1_"${sample_id}"_hdr.sam \  
	"${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_read_pairs_for_NOTchr1o11.sam \  
	"${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam \  
	> "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11.sam  
samtools view -bS "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11.sam \  
	> "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11.bam  
```

Clean up
```
echo "EXECUTE end of processing, don't remove files"
rm "${output_directory}"/MGRBp1_"${sample_id}"_hdr.sam
rm "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_nohdr.sam
rm "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_nohdr.sam
rm "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_IDs.txt
rm "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_IDs.txt
rm "${output_directory}"/MGRBp1_"${sample_id}"_NOTchr1to11_read_pairs_for_chr1o11.sam
rm "${output_directory}"/MGRBp1_"${sample_id}"_chr1to11_read_pairs_for_NOTchr1o11.sam
```

Then run the Mobster pipeline above.
Then merge the predictions into one file. 
Merging predictions is necessary because MEIs from discordant read pairs appearing in more than one BAM file will be called for each BAM file,
so need to merge, choosing the one with the best support (most reads seen) so that this MEI appears only once in the final merged file.
Here is an example.
```
python merge_Mobster_predictions.py output_Mobster_predictions.txt input_Mobster_1_predictions.txt input_Mobster_2_predictions.txt input_Mobster_2_predictions.txt ...
```

