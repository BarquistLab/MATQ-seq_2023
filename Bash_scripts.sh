##############################
#                            #
#     MATQ-seq - 2023        #
#                            #
##############################


#Used to count the number of reads before and after trimming
parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *.gz


#Within primers file (matseq_primers.fa)
>GAT27_dt
GTGAGTGATGGTTGAGGATGTGTGGAGNNNNN
>GAT275N
GTGAGTGATGGTTGAGGATGTGTGGAGNNNNN
>GAT21_6N3G
GATGGTTGAGGATGTGTGGAGNNNNNNGGG
>GAT27_PCR
GTGAGTGATGGTTGAGGATGTGTGGAG

#nextera_and_primers.fa
#The sequences above added to the default bbduk adapters file


#---------
#
# BBDuk
#
#---------

dir=~/Single_cell/

mkdir -p $dir/BBDuk_L
mkdir -p $dir/BBDuk_L_R

for i in {1..384}
  do

  report="$(find $dir/sequencing_reads -name '*.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub(".fastq.gz", "", $0); print}')"
  r1="${report}.fastq.gz"


  #trim left
  bbduk.sh -Xmx10g t=20 in=$dir/sequencing_reads/$r1 out=$dir/BBDuk_L/$r1 ref=$dir/matseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 trimpolya=30

  #Trim right
  bbduk.sh -Xmx10g t=20 in=$dir/BBDuk_L/$r1 out=$dir/BBDuk_L_R/$r1 ref=$dir/nextera_and_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=r k=17 mink=11 hdist=1

done



#---------
#
# Bowtie2
#
#---------

dir=~/Single_cell/
mkdir -p $dir/bowtie2_aligned

#make index
bowtie2-build salmonella_sl1344.fa sl1344

for i in {1..384}
  do

  report="$(find $dir/BBDuk_L_R -name '*.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub(".fastq.gz", "", $0); print}')"
  read1="${report}.fastq.gz"
  outfile="${report}.bam"

  #Align
  bowtie2 -p 30 --local -x $dir/sl1344 -N 1 -U $dir/BBDuk_L_R/$read1 | samtools view -@ 20 -bS -h - > $dir/bowtie2_aligned/$outfile
done





#---------
#
# featureCounts
#
#---------

#specify version 2.0.1 when installing

dir=~/Single_cell/
mkdir -p $dir/featureCounts


for i in {1..384}
do

  report="$(find $dir/BBDuk_L_R/ -name '*.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub(".fastq.gz", "", $0); print}')"
  bam_file="${report}.bam"
  count_file="${report}.count"

  #Align
  featureCounts -T 10 -a $dir/salmonella_sl1344.gff \
  -t "CDS,sRNA,tRNA,rRNA,pseudogene,3UTR,5UTR" -g "Name" \
  -o $dir/featureCounts/$count_file $dir/bowtie2_aligned/$bam_file
done


