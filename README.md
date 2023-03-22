# kmerDedup
Reduce genome assembly redundancy using shared mapped k-mer method.
This pipeline was inspired by khaper (https://github.com/tangerzhang/khaper, https://github.com/lardo/khaper) but withoud k-mer size and genome size limitations, and I have never seen the source code of the core program "Compress" in khaper pipeline.

Use cases
1. assembling software (flye, canu, ipa, hifiasm, masurca, wtdbg2 or any others) generate larger genome size than expected and purge_dups (https://github.com/dfguan/purge_dups) can not work well (no distinct depth peaks when using purge_dups or purge_haplotigs);
2. pooled sample assembly (assembling software may give you 2-4X larger genome size than expected);
3. there are several assemblies generated from different softwares, and you may want to cat them together and choose longer non-redundancy contig/scaffold set.

Notes
1. polish your assembly first if it was generated from long noisy reads;
2. don't use long noisy reads in this pipeline. You can use pacbio-hifi reads and/or ngs reads (and/or ont duplex reads).

Limits
1. only work for diploid, I think;
2. can not deal with chimera contigs/scaffolds because kmerDedup do not break and don't know how to break chimera.

There are 4 steps to run this pipeline, 
1. use jellyfish2 to count and dump k-mer sequences, than filter the dump file using kmerFilter.pl script; 
2. map k-mer sequences to contigs/scaffolds using bowtie2 and use BamDeal to calculate coverage;
3. use kmerDedup.pl (core program in this pipeline) to select longest contig/scaffold set;
4. (optional) use purge_dups to purge kmerDedup.pl results if the kmerDedup.pl result is still larger than expected.

# Dependences
1. software

perl https://www.perl.org

parallel https://www.gnu.org/software/parallel

jellyfish2 https://github.com/gmarcais/Jellyfish

bowtie2 https://github.com/BenLangmead/bowtie2

BamDeal https://github.com/BGI-shenzhen/BamDeal

samtools http://www.htslib.org

pigz https://github.com/madler/pigz

optional software (if you want to run purge_dups after kmerDedup)

purge_dups https://github.com/dfguan/purge_dups

bwa-mem2 (for short reads mapping) https://github.com/bwa-mem2/bwa-mem2

minimap2 (for long reads mapping) https://github.com/lh3/minimap2

2. perl modules

use strict;

use warnings;

use Getopt::Long;

use File::Spec;

use File::Basename qw(basename dirname);

use FindBin qw($Bin $Script);

use Data::Dumper;

use List::Util qw{any first};

use Storable;

use Bio::SeqIO;

# Testing data
The writer of pseudohaploid (https://github.com/schatzlab/pseudohaploid) provide a good testing data. You can download canu results from this link (http://labshare.cshl.edu/shares/schatzlab/www-data/pseudohaploid/arabidopsis) and download corresponding short reads data from NCBI (SRR3703081, SRR3703082, SRR3703105).
```
wget https://labshare.cshl.edu/shares/schatzlab/www-data/pseudohaploid/arabidopsis/rawcanu.fa.gz && gunzip rawcanu.fa.gz
```
and you will get rawcanu.fa

Choose your best way to download SRR3703081, SRR3703082, SRR3703105 from NCBI, then

```
cat SRR3703081_1.fastq.gz SRR3703082_1.fastq.gz SRR3703105_1.fastq.gz > ath_R1.fq.gz
or
zcat SRR3703081_1.fastq.gz SRR3703082_1.fastq.gz SRR3703105_1.fastq.gz |pigz -p 6 -c > ath_R1.fq.gz

cat SRR3703081_2.fastq.gz SRR3703082_2.fastq.gz SRR3703105_2.fastq.gz > ath_R2.fq.gz
or
zcat SRR3703081_2.fastq.gz SRR3703082_2.fastq.gz SRR3703105_2.fastq.gz |pigz -p 6 -c > ath_R2.fq.gz
```

# step 1 count and filter k-mers
k-mer size, you can use 17 19 21 23 or others, and you can try difference sizes if you want, I don't known what size is the best so far, I think all sizes may work.
## 1.1 count k-mers
```
/Bio/software/anaconda3/envs/jellyfish/bin/jellyfish count -m 21 -s 10G -t 30 -c 8 -C /dev/fd/0 --timing ath.time -o ath.count <(/Bio/software/anaconda3/envs/pigz/bin/pigz -p 6 -d -c ath_R1.fq.gz) <(/Bio/software/anaconda3/envs/pigz/bin/pigz -p 6 -d -c ath_R2.fq.gz)
```
if you have high quality long reads, feel free to feed into jellyfish count
```
/Bio/software/anaconda3/envs/jellyfish/bin/jellyfish count -m 21 -s 10G -t 30 -c 8 -C /dev/fd/0 --timing ath.time -o ath.count <(/Bio/software/anaconda3/envs/pigz/bin/pigz -p 6 -d -c ath_R1.fq.gz) <(/Bio/software/anaconda3/envs/pigz/bin/pigz -p 6 -d -c ath_R2.fq.gz) <(cat ath_hifi.fa)
```

if jellyfish gives you more than one count file, you need merge them first
```
if [ -e ath.count_1 ]; then
  /Bio/software/anaconda3/envs/jellyfish/bin/jellyfish merge -v -o ath.count.jf ath.count_*
else
  mv ath.count ath.count.jf
fi
```

## 1.2 stat and histo (you can skip)
```
/Bio/software/anaconda3/envs/jellyfish/bin/jellyfish stats -o ath.stats ath.count.jf
/Bio/software/anaconda3/envs/jellyfish/bin/jellyfish histo -t 30 ath.count.jf|/Bio/bin/perl -lane 'my ($dpt, $cnt) = split(/\s+/, $_); my $nn = $dpt * $cnt;print "$dpt\t$cnt\t$nn"' > ath.histo
```
## 1.3 dump k-mers
-l 3 (no inculding, it means all k-mers with depth <= 3 will be deleted) is suitable for many case, -u can be very large, I think 10000-50000 may suitable for many case unless you are dealing with a highly repetitive genome.
You don't need to split ath.filt.fa if you don't want to run bowtie2 parallel.
```
/Bio/software/anaconda3/envs/jellyfish/bin/jellyfish dump -c -t -o ath.dump.all ath.count.jf
/Bio/bin/perl kmerDedup/kmerFilter.pl -d ath.dump.all -o ath.filt.fa -l 3 -u 10000
/Bio/bin/perl kmerDedup/splitFasta.pl -f ath.filt.fa -o split -k ath.kmer
```

# step 2 mapping k-mers
Use fa2fa.pl to format contig file (especially when you cat difference software results together), -l is use to limit minimum contig length, contig with length lower than this value will be discarded.
## 2.1 build reference index
```
/Bio/bin/perl kmerDedup/fa2fa.pl -f rawcanu.fa -o ath.format.fa -c F -n F -l 1000
source /Bio/software/anaconda3/envs/bowtie2/../../bin/activate bowtie2;/Bio/software/anaconda3/envs/bowtie2/bin/bowtie2-build --threads 30 ath.format.fa ath.format
```
## 2.2 mapping k-mers to reference
Use absolute path if you want to parallelize bowtie2 task submit. -k can be larger if you are dealing with a highly repetitive genome but larger -k may slow down the mapping step and the kmerDedup step.
-L can be equal to k-mer size or smaller than k-mer size.
Please don't change --very-sensitive --score-min L,-0.6,-0.2 --end-to-end or add other parameters to allow mismatch unless you really know what you are doing.
Don't forget to use samtools view -F 4 to filter out unmap k-mers.
```
ls split/|grep ".fa$"|sort -V |sed 's/.fa$//g'|/Bio/bin/perl -lane 'print "source /Bio/software/anaconda3/envs/bowtie2/../../bin/activate bowtie2;/Bio/User/anaconda3/envs/bowtie2/bin/bowtie2 --very-sensitive -k 1000 --score-min L,-0.6,-0.2 --end-to-end --reorder -L 21 --rg-id ath --rg SM:ath -p 6 -f /WORK/Bio/Project/pipe/assembly/dedup/pipe/02.kmer/split/$_.fa -x /WORK/Bio/Project/pipe/assembly/dedup/pipe/01.ref/ath.format |/Bio/bin/samtools-1.14 view -@ 4 -F 4 -bS - > /WORK/Bio/Project/pipe/assembly/dedup/pipe/03.mapping/$_.bam"'> bowtie2.sh
nohup /Bio/bin/parallel -j 10 < bowtie2.sh > bowtie2.sh.log &
or
submit bowtie2.sh task parallelize when you are using a cluster system.
```

## 2.3 merge mapping results and calculate contig depth
please don't forget -n parameter when using samtools merge, or you will need samtools sort -n to sort bam file by read names.
```
ls /WORK/Bio/Project/pipe/assembly/dedup/pipe/03.mapping/*.bam > /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.bamlist.xls
/Bio/bin/samtools-1.14 merge -@ 30 -n -f -b /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.bamlist.xls /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.sort.bam
/Bio/software/BamDeal-0.27/bin/BamDeal_Linux statistics Coverage -i /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.sort.bam -r /WORK/Bio/Project/pipe/assembly/dedup/pipe/01.ref/ath.format.fa -q 0 -o /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.cov
```

# step 3 kmerDedup

## 3.1 run kmerDedup.pl from bam file

```
/Bio/bin/perl kmerDedup/kmerDedup.pl
	Description: remove redundancy using k-mer shared percentage

	-f    <file>   fasta file                  [force]
	-k    <str>    output prefix               [force]
	-o    <dir>    output dir                  [force]
	-bam  <file>   bam file                    [force if no -d]
	-dum  <file>   path to dump hash           [force if no -b]
	-cov  <file>   coverage file from BamDeal  [optional]
	-wtl  <file>   whitelist to keep           [optional]
	-bll  <file>   blacklist to remove         [optional]
	-mpr  <float>  max duplication percentage  [0.3]
	-mcv  <int>    min k-mer coverage(%)       [30]
	-kmer <int>    k-mer size                  [17]
	-mode <1/2>    1:ratio only; 2:ratio * cov [2]
	-s    <bin>    samtools binary             [/Bio/bin/samtools-1.14]
	-t    <int>    threads for samtools view   [4]
	-h           Help document

```
You can try -mpr 0.3 first.
```
/Bio/bin/perl kmerDedup/kmerDedup.pl -k ath -mpr 0.3 -mcv 30 -kmer 21 -o 05.kmerdedup -f /WORK/Bio/Project/pipe/assembly/dedup/pipe/01.ref/ath.format.fa -bam /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.sort.bam -cov /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.cov.stat -s /Bio/bin/samtools-1.14
```
Running time dependant on bam line number and how many scaffod need to compare, for this example data (200M genome (1.7X to expected), about 206 M bam lines), it cost 15 minutes to get results. For a 0.62G genome and 2X redundancy, 854 M bam lines, it cost 2.5 hours to get results. For a 1G genome and 5.5X redundancy (cat difference software results together or pooled sample), ~4000 M bam lines, it may cost ~30 hours to get results, but worth waiting.

Results
```
├── ath.all.xls     ## all contig records
├── ath.bam.log     ## progressing log
├── ath.drop.fa     ## droped contig in fasta format
├── ath.dump.hash   ## mapped k-mer count and share mapped k-mer count hash, use this file if you want to adjust -mpr
├── ath.keep.fa     ## keep contig in fasta format
└── ath.keep.xls    ## same as ath.all.xls but only kept contigs
```

## 3.2 adjust -mpr
Use 05.kmerdedup/ath.dump.hash instead of bam to save running time.
Don't use too small -mpr, you can try difference values, but stop using smaller value when BUSCO C value decreasing dramatically, if kmerDedup result still larger than expected, use largest BUSCO C version to run purge_dups (for example, raw assembly BUSCO C is 98%, D is 58%; -mpr 0.25 BUSCO C is 92%, D is 5%; -mpr 0.3 BUSCO C is 94%, D is 6%; -mpr 0.4 BUSCO C is 97%, D is 10%, then I will choose -mpr 0.4 results to run purge_dups).
```
/Bio/bin/perl kmerDedup/kmerDedup.pl -k ath -mpr 0.5 -mcv 30 -kmer 21 -o 05.kmerdedup0.5 -f /WORK/Bio/Project/pipe/assembly/dedup/pipe/01.ref/ath.format.fa -dum 05.kmerdedup/ath.dump.hash -cov /WORK/Bio/Project/pipe/assembly/dedup/pipe/04.merge/ath.cov.stat -s /Bio/bin/samtools-1.14 -s /Bio/bin/samtools-1.14
```

# Example results

| Nx      | rawcanu.fa | pseudohap.fa | kmerDedup-mpr 0.3 -mode 2 | kmerDedup -mpr 0.5 -mode 2 |
|---------|------------|--------------|---------------------------|----------------------------|
| Total   | 214700023  | 143490505    | 123325563                 | 130255158                  |
| Count   | 2074       | 505          | 216                       | 311                        |
| Average | 103519.78  | 284139.61    | 570951.68                 | 418826.87                  |
| Median  | 41735      | 82149        | 296376                    | 138132                     |
| N00     | 4128881    | 4128881      | 4128881                   | 4128881                    |
| N10     | 2878449    | 3214456      | 3214456                   | 3214456                    |
| N20     | 1852867    | 2412251      | 2436491                   | 2436491                    |
| N30     | 1128124    | 1852867      | 2120621                   | 1949687                    |
| N40     | 715569     | 1204449      | 1599836                   | 1470490                    |
| N50     | 350182     | 950253       | 1164438                   | 1109900                    |
| N60     | 138132     | 715569       | 934741                    | 900377                     |
| N70     | 75340      | 424773       | 703219                    | 639541                     |
| N80     | 50143      | 259779       | 453978                    | 388804                     |
| N90     | 36026      | 107428       | 279847                    | 199957                     |
| N100    | 12153      | 13139        | 13139                     | 13139                      |


BUSCO V 5.3.1 embryophyta_odb10

| Method                              | rawcanu        | pseudohap      | kmerDedup -mpr 0.3 -mode 2 | kmerDedup -mpr 0.5 -mode 2 |
|-------------------------------------|----------------|----------------|----------------------------|----------------------------|
| Complete BUSCOs (C)                 | 1598 (99.01%)  | 1593 (98.70%)  | 1584 (98.14%)              | 1591 (98.57%)              |
| Complete and single-copy BUSCOs (S) | 812 (50.31%)   | 1456 (90.21%)  | 1537 (95.23%)              | 1531 (94.86%)              |
| Complete and duplicated BUSCOs (D)  | 786 (48.70%)   | 137 (8.49%)    | 47 (2.91%)                 | 60 (3.72%)                 |
| Fragmented BUSCOs (F)               | 8 (0.50%)      | 8 (0.50%)      | 11 (0.68%)                 | 10 (0.62%)                 |
| Missing BUSCOs (M)                  | 8 (0.50%)      | 13 (0.81%)     | 19 (1.18%)                 | 13 (0.81%)                 |
| Total BUSCO groups searched         | 1614 (100.00%) | 1614 (100.00%) | 1614 (100.00%)             | 1614 (100.00%)             |


# step 4 run purge_dups if really needed.

# cites
I am working in industry and my english is so so, so I am not so professional to write papers, and there is no related paper in the near further. You can cite this github repository directly, but in your own risk.

You may need to cite the corresponding papers of the software which you have used in this pipeline.
