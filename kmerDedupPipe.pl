#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Data::Dumper;
use List::Util qw{any first};
############# GetOptions #################
my ($outdir, $key, $config_file, $shdir);
my $run;
GetOptions(
	"h|?"=>\&help,
	"c:s"=>\$config_file,
	"k:s"=>\$key,
	"o:s"=>\$outdir,
	"s:s"=>\$shdir,
	"run:s"=>\$run,
) || &help;
&help unless ($outdir && $key && $config_file);

sub help
{
	print"
	Description: kmer-base redundancy removing pipeline

	-c    <file>  config file           [force]
	-k    <file>  output prefix         [force]
	-o    <dir>   output directory      [force]
	-s    <dir>   shell directory       [-o]
	-run  <T/F>   run or not            [F]
	-h            Help document
";
	exit;
}
=pod
Tag Meaning
NM     Edit distance
MD     Mismatching positions/bases
AS     Alignment score
BC     Barcode sequence
X0     Number of best hits
X1     Number of suboptimal hits found by BWA
XN     Number of ambiguous bases in the reference
XM     Number of mismatches in the alignment
XO     Number of gap opens
XG     Number of gap extentions
XT     Type: Unique/Repeat/N/Mate-sw
XA     Alternative hits; format: (chr,pos,CIGAR,NM;)*
XS     Suboptimal alignment score
XF     Support from forward/reverse alignment
XE     Number of supporting seeds
=cut
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$config_file = &abs_dir($config_file);

if(!-d $outdir){
	mkdir $outdir;
}
$outdir = &abs_dir($outdir);

$shdir ||= "$outdir/shell";
if(!-d $shdir){
	mkdir $shdir;
}
$shdir = &abs_dir($shdir);

my %config;
&read_config($config_file, \%config);

my %readfile;
my @gzfile;
my @uzfile;
open(IN, $config{readlist});
while (<IN>) {
	chomp;
	my ($type, $read, $file) = split(/\s+/, $_);
	$type = lc($type);
	$read = lc($read);
	if($file =~ /.gz$/){
		push(@gzfile, $file);
	}else{
		push(@uzfile, $file);
	}
	$readfile{$type}{$read} = $file;
}
close(IN);

&main();

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&Runtime($begin_time);
##########################################
sub main()
{
	&make_dir();
	my $shrun = "$shdir/step_by_step.sh";
	open(PIPE, ">$shrun");
	print PIPE "set -e\n";
	my $cmd = &step01_ref();
	print PIPE "$cmd\n";
	$cmd = &step02_kmer();
	print PIPE "$cmd\n";
	$cmd = &step03_mapping();
	print PIPE "$cmd\n";
	$cmd = &step04_merge();
	print PIPE "$cmd\n";
	$cmd = &step05_kmerdedup();
	print PIPE "$cmd\n";
	if($config{purgeDups} eq "T"){
		$cmd = &step06_purgeDups();
		print PIPE "$cmd\n";
	}
	close(PIPE);
}
sub make_dir()
{
	mkdir "$outdir/01.ref";
	mkdir "$outdir/02.kmer";
	mkdir "$outdir/03.mapping";
	mkdir "$outdir/04.merge";
	mkdir "$outdir/05.kmerdedup";
	mkdir "$outdir/06.purgeDups";

	mkdir "$shdir/step01";
	mkdir "$shdir/step02";
	mkdir "$shdir/step03";
	mkdir "$shdir/step04";
	mkdir "$shdir/step05";
	mkdir "$shdir/step06";
}

#
sub step01_ref()
{
	my $shfa2fa = "$shdir/step01/q01.fa2fa.sh";
	open(SH, ">$shfa2fa");
	my $cmd = "$config{perl} $Bin/kmerDedup/fa2fa.pl -f $config{genome} ";
	$cmd .= "-o $outdir/01.ref/$key.format.fa -c F -n F -l $config{mlength}";
	print SH "$cmd\n";
	close(SH);

	my $shidx = "$shdir/step01/q02.buildIndex.sh";
	open(SH, ">$shidx");
	$cmd = "source $config{bowtie2}/../../bin/activate $config{bowtie2};";
	$cmd .= "$config{bowtie2}/bin/bowtie2-build --threads $config{thread_large} $outdir/01.ref/$key.format.fa ";
	$cmd .= "$outdir/01.ref/$key.format";
	print SH "$cmd\n";

	my $shrun = "$shdir/step01.ref.sh";
	open(SH, ">$shrun");
	print SH "bash $shfa2fa\n";
	print SH "bash $shidx\n";
	close(SH);
	$cmd = "bash $shrun";
	return($cmd);
}

sub step02_kmer()
{
	my $kmerdir = "$outdir/02.kmer";
	my $shkmer = "$shdir/step02/q01.jfcount.sh";
	open(SH, ">$shkmer");
	my $cmd = "$config{jellyfish} count -m $config{kmer} -s 10G ";
	$cmd .= "-t $config{thread_large} -c 8 -C /dev/fd/0 ";
	$cmd .= "--timing $kmerdir/$key.time -o $kmerdir/$key.count ";
	if(@gzfile > 0){
		my $gfile = "";
		foreach my $zz (@gzfile) {
			$gfile .= "<($config{pigz} -p $config{thread_small} -d -c $zz) ";
		}
		$cmd .= $gfile;
	}
	if(@uzfile > 0){
		my $ufile = "";
		foreach my $uu (@uzfile) {
			$ufile .= "<(cat $uu) ";
		}
		$cmd .= $ufile;
	}
	print SH "$cmd\n";
	close(SH);

	my $shmerge = "$shdir/step02/q02.jfmerge.sh";
	open(SH, ">$shmerge");
	$cmd = "if [ -e $kmerdir/$key.count_1 ]; then\n";
	$cmd .= "  $config{jellyfish} merge -v -o $kmerdir/$key.count.jf $kmerdir/$key.count_*\n";
	$cmd .= "else\n";
	$cmd .= "  mv $kmerdir/$key.count $kmerdir/$key.count.jf\n";
	$cmd .= "fi";
	print SH "$cmd\n";
	close(SH);

	my $shstats = "$shdir/step02/q03.jfstats.sh";
	open(SH, ">$shstats");
	$cmd = "$config{jellyfish} stats -o $kmerdir/$key.stats $kmerdir/$key.count.jf";
	print SH "$cmd\n";
	close(SH);

	my $shhisto = "$shdir/step02/q04.jfhisto.sh";
	open(SH, ">$shhisto");
	$cmd = "$config{jellyfish} histo -t $config{thread_large} $kmerdir/$key.count.jf|";
	$cmd .= "$config{perl} -lane \'my (\$dpt, \$cnt) = split(/\\s+/, \$_); my \$nn = \$dpt * \$cnt;";
	$cmd .= "print \"\$dpt\\t\$cnt\\t\$nn\"\' > $kmerdir/$key.histo";
	print SH "$cmd\n";
	close(SH);

	my $shdump = "$shdir/step02/q05.jfdump.sh";
	open(SH, ">$shdump");
	$cmd = "$config{jellyfish} dump -c -t -o $kmerdir/$key.dump.all $kmerdir/$key.count.jf";
	print SH "$cmd\n";
	close(SH);

	my $shfilter = "$shdir/step02/q06.dumpFilter.sh";
	open(SH, ">$shfilter");
	$cmd = "$config{perl} $Bin/kmerDedup/kmerFilter.pl -d $kmerdir/$key.dump.all ";
	$cmd .= "-o $kmerdir/$key.filt.fa -l $config{lowkmer} -u $config{upkmer}";
	print SH "$cmd\n";
	close(SH);

	my $shsplit = "$shdir/step02/q07.split.sh";
	open(SH, ">$shsplit");
	$cmd = "$config{perl} $Bin/kmerDedup/splitFasta.pl -f $kmerdir/$key.filt.fa ";
	$cmd .= "-o $kmerdir/split -k $key.kmer";
	print SH "$cmd\n";
	close(SH);

	my $shrun = "$shdir/step02.kmer.sh";
	open(SH, ">$shrun");
	print SH "set -e\n";
	print SH "bash $shkmer\n";
	print SH "bash $shmerge\n";
	print SH "bash $shstats\n";
	print SH "bash $shhisto\n";
	print SH "bash $shdump\n";
	print SH "bash $shfilter\n";
	print SH "bash $shsplit\n";
	close(SH);
	$cmd = "bash $shrun";
	return($cmd);
}

sub step03_mapping()
{
	my $shcmd = "$shdir/step03/q01.cmd.sh";
	open(SH, ">$shcmd");
	my $cmd = "ls $outdir/02.kmer/split|grep \".fa\$\"|sed \'s/.fa\$//g\'|$config{perl} -lane \'print \"";
	$cmd .= "source $config{bowtie2}/../../bin/activate $config{bowtie2};";
	$cmd .= "$config{bowtie2}/bin/bowtie2 $config{lpara} -L $config{kmer} --rg-id $key --rg SM:$key ";
	$cmd .= "-p $config{thread_small} -f $outdir/02.kmer/split/\$_\.fa -x $outdir/01.ref/$key.format ";
	$cmd .= "|$config{samtools} view -\@ 4 -F 4 -bS - > $outdir/03.mapping/\$_\.bam";
	$cmd .= "\"\'> $shdir/step03/bowtie2.sh";
	print SH "$cmd\n";
	close(SH);

	my $shqsub = "$shdir/step03/q02.qsub.sh";
	open(SH, ">$shqsub");
	$cmd = "$config{perl} $config{qsub} --convert no --queue $config{queue} ";
	$cmd .= "--maxjob $config{maxjob} $shdir/step03/bowtie2.sh";
	print SH "$cmd\n";
	close(SH);
	my $shrun = "$shdir/step03.mapping.sh";
	open(SH, ">$shrun");
	print SH "set -e\n";
	print SH "bash $shcmd\n";
	print SH "bash $shqsub\n";
	close(SH);
	$cmd = "bash $shrun";
	return($cmd);
}

sub step04_merge()
{
	my $shlist = "$shdir/step04/q01.list.sh";
	open(SH, ">$shlist");
	my $cmd = "ls $outdir/03.mapping/\*.bam > $outdir/04.merge/$key.bamlist.xls";
	print SH "$cmd\n";
	close(SH);

	my $shmerge = "$shdir/step04/q02.merge.sh";
	open(SH, ">$shmerge");
	$cmd = "$config{samtools} merge -\@ $config{thread_large} -n ";
	$cmd .= "-f -b $outdir/04.merge/$key.bamlist.xls -o $outdir/04.merge/$key.sort.bam";
	print SH "$cmd";
	close(SH);

	my $shcov = "$shdir/step04/q03.bamdeal.sh";
	open(SH, ">$shcov");
	$cmd = "$config{bamdeal} statistics Coverage -i $outdir/04.merge/$key.sort.bam ";
	$cmd .= "-r $outdir/01.ref/$key.format.fa -q 0 -o $outdir/04.merge/$key.cov";
	print SH "$cmd\n";
	close(SH);
	my $shrun = "$shdir/step04.merge.sh";
	open(SH, ">$shrun");
	print SH "set -e\n";
	print SH "bash $shlist\n";
	print SH "bash $shmerge\n";
	print SH "bash $shcov\n";
	close(SH);
	$cmd = "bash $shrun";
	return($cmd);
}

sub step05_kmerdedup()
{
	my $shdedup = "$shdir/step05/q01.kmerDedup.sh";
	open(SH, ">$shdedup");
	my $cmd = "$config{perl} $Bin/kmerDedup/kmerDedup.pl ";
	$cmd .= "-k $key -mpr $config{mindup} -mcv $config{mkcov} -kmer $config{kmer} ";
	$cmd .= "-o $outdir/05.kmerdedup -f $outdir/01.ref/$key.format.fa ";
	$cmd .= "-bam $outdir/04.merge/$key.sort.bam -cov $outdir/04.merge/$key.cov.stat ";
	if(-e $config{whitelist}){
		$cmd .= "-wtl $config{whitelist} ";
	}
	if(-e $config{blacklist}){
		$cmd .= "-bll $config{blacklist} ";
	}
	$cmd .= "-s $config{samtools}";
	print SH "$cmd\n";
	close(SH);
	my $shrun = "$shdir/step05.kmerdedup.sh";
	open(SH, ">$shrun");
	print SH "bash $shdedup\n";
	close(SH);
	$cmd = "bash $shrun";
	return($cmd);
}

sub step06_purgeDups()
{
	my $shref = "$shdir/step06/q01.ref.sh";
	open(SH, ">$shref");
	my $cmd = "ln -sf $outdir/05.kmerdedup/$key.keep.fa $outdir/06.purgeDups/$key.kdedup.fa";
	print SH "$cmd\n";
	$cmd = "$config{purge_dups}/bin/split_fa $outdir/06.purgeDups/$key.kdedup.fa > $outdir/06.purgeDups/$key.split.fa";
	print SH "$cmd\n";
	if($config{purgeData} eq "ngs"){
		$cmd = "$config{bwamem2} index $outdir/06.purgeDups/$key.kdedup.fa";
		print SH "$cmd\n";
	}
	close(SH);

	my $shmsplit = "$shdir/step06/q02.mapsplit.sh";
	open(SH, ">$shmsplit");
	$cmd = "$config{minimap2} -t $config{thread_large} -I 64G -xasm5 -DP $outdir/06.purgeDups/$key.split.fa ";
	$cmd .= "$outdir/06.purgeDups/$key.split.fa|$config{pigz} -p $config{thread_small} -c - ";
	$cmd .= "> $outdir/06.purgeDups/$key.split.paf.gz";
	print SH "$cmd\n";
	close(SH);

	my $shmapreads = "$shdir/step06/q03.mapreads.sh";
	open(SH, ">$shmapreads");
	if($config{purgeData} eq "ngs"){
		$cmd = "$config{bwamem2} mem -t $config{thread_large} -M $outdir/06.purgeDups/$key.kdedup.fa ";
		$cmd .= "$readfile{ngs}{r1} $readfile{ngs}{r2} |";
		$cmd .= "$config{samtools} fixmate -m -\@ $config{thread_small} - ";
		$cmd .= "$outdir/06.purgeDups/$key.kdedup.fixmate.bam";
		print SH "$cmd\n";
	}else{
		my $preset = "map-hifi";
		if($config{tgstype} eq "clr"){
			$preset = "map-pb";
		}
		if($config{tgstype} eq "ont"){
			$preset = "map-ont";
		}
		$cmd = "$config{minimap2} -t $config{thread_large} -I 64G -x $preset ";
		$cmd .= "$outdir/06.purgeDups/$key.kdedup.fa $readfile{tgs}{r1} |";
		$cmd .= "$config{pigz} -p $config{thread_small} -c - > ";
		$cmd .= "$outdir/06.purgeDups/$key.tgs.paf.gz";
		print SH "$cmd\n";
	}
	close(SH);

	my $shcov = "$shdir/step06/q04.stat.sh";
	my $statfile;
	my $basecov;
	open(SH, ">$shcov");
	print SH "cd $outdir/06.purgeDups\n";
	if($config{purgeData} eq "ngs"){
		$cmd = "$config{purge_dups}/bin/ngscstat $outdir/06.purgeDups/$key.kdedup.fixmate.bam";
		print SH "$cmd\n";
		$statfile = "$outdir/06.purgeDups/TX.stat";
		$basecov = "$outdir/06.purgeDups/TX.base.cov";
	}else{
		$cmd = "$config{purge_dups}/bin/pbcstat $outdir/06.purgeDups/$key.tgs.paf.gz";
		print SH "$cmd\n";
		$statfile = "$outdir/06.purgeDups/PB.stat";
		$basecov = "$outdir/06.purgeDups/PB.base.cov";
	}
	close(SH);

	my $shcutoff = "$shdir/step06/q05.cutoffs.sh";
	open(SH, ">$shcutoff");
	$cmd = "$config{purge_dups}/bin/calcuts $statfile ";
	$cmd .= "> $outdir/06.purgeDups/cutoffs 2>$outdir/06.purgeDups/calcults.log";
	print SH "$cmd\n";
	$cmd = "$config{purge_dups}/bin/python $config{purge_dups}/bin/hist_plot.py ";
	$cmd .= "-c $outdir/06.purgeDups/cutoffs $statfile ";
	$cmd .= "$outdir/06.purgeDups/$key.kdedup.cov.pdf";
	print SH "$cmd\n";
	close(SH);

	my $shpurge = "$shdir/step06/q06.purge_dups.sh";
	open(SH, ">$shpurge");
	$cmd = "$config{purge_dups}/bin/purge_dups -2 -T $outdir/06.purgeDups/cutoffs ";
	$cmd .= "-c $basecov $outdir/06.purgeDups/$key.split.paf.gz ";
	$cmd .= " > $outdir/06.purgeDups/dups.bed 2>$outdir/06.purgeDups/purge_dups.log";
	print SH "$cmd\n";
	close(SH);

	my $shgetseq = "$shdir/step06/q07.get_seqs.sh";
	open(SH, ">$shgetseq");
	print SH "cd $outdir/06.purgeDups\n";
	$cmd = "$config{purge_dups}/bin/get_seqs -e $outdir/06.purgeDups/dups.bed ";
	$cmd .= "$outdir/06.purgeDups/$key.kdedup.fa";
	print SH "$cmd\n";
	$cmd = "cat |perl -lane \'if(/>/){\$_ =~ s/_1\$//g; print \$_}else{print \$_}' > $key.purged.fa";
	close(SH);
	my $shrun = "$shdir/step06.purge_dups.sh";
	open(SH, ">$shrun");
	print SH "set -e\n";
	print SH "bash $shref\n";
	print SH "bash $shmsplit\n";
	print SH "bash $shmapreads\n";
	print SH "bash $shcov\n";
	print SH "bash $shcutoff\n";
	print SH "bash $shpurge\n";
	print SH "bash $shgetseq\n";
	close(SH);
	$cmd = "bash $shrun\n";
	return($cmd);
}

#
sub read_config()
{
	my ($file, $config) = @_;
	my $marker_flag = 0;
	open(IN, $file);
	while (<IN>) {
		chomp;
		next if(/^$/);
		next if(/^\s+$/);
		next if(/^\#/);
		next if(/^\r/);
		my ($key, $value) = split(/:=/, $_);
		$key =~ s/^\s+|\s+$//g;
		$value =~ s/^\s+|\s+$//g;
		$$config{$key} = $value;
	}
	close(IN);
}

#
sub formatSeq {
	my ($seq_ref, $len) = @_;
	$len ||= 80;
	$$seq_ref =~ s/(.{$len})/$1\n/g;
	chomp $$seq_ref;
}
# sub date format
sub date_format()
{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec);
}

# sub Runtime
sub Runtime()
{
	my ($begin_time) = @_;
	my $now_time = time();
	my $total_time = $now_time - $begin_time;
	my $sec = 0; my $minu = 0; my $hour = 0;
	if($total_time >= 3600){
		$hour = int($total_time/3600);
			my $left = $total_time % 3600;
		if($left >= 60){
			$minu = int($left/60);
			$sec = $left % 60;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		} else {
			$minu = 0;
			$sec = $left;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		}
	} else {
		if($total_time >= 60){
			$minu = int($total_time/60);
			$sec = $total_time % 60;
			$total_time = $minu."m\-".$sec."s";
		} else {
			$sec = $total_time;
			$total_time = $sec."s";
		}
	}
	print "Total elapsed time [$total_time]\n\n";
}

# sub absolutely directory
sub abs_dir()
{
	my ($in) = @_;
	my $current_dir = `pwd`;
	chomp($current_dir);
	my $return_dir = "";
	if(-f $in){
		my $in_dir = dirname($in);
		my $in_file = basename($in);
		chdir $in_dir;
		$in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = "$in_dir/$in_file";
	} elsif(-d $in) {
		chdir $in;
		my $in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = $in_dir;
	} else {
		die("ERROR: there is no file or dir called [$in], please check!\n");
	}
	chdir $current_dir;
	return $return_dir;
}

# show log
sub show_log()
{
	my ($text) = @_;
	my $current_time = &date_format(localtime());
	print "$current_time: $text\n";
}

# run or die
sub run_or_die()
{
	my ($cmd) = @_;
	&show_log($cmd);
	my $flag = system($cmd);
	if($flag != 0){
		&show_log("ERROR: command $cmd");
		exit(1);
	}
	&show_log("done\n");
}

# qsub
sub qsub()
{
	my ($sh, $vf, $maxjob, $queue) = @_;
	$vf ||= "2G";
	$maxjob ||= 50;
	$queue ||= "all.q";
	my $cmd_qsub = "/Bio/bin/qsub-sge.pl --convert no --queue $queue --maxjob $maxjob --resource vf=$vf $sh";
	&run_or_die($cmd_qsub);
}
