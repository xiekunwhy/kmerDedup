#!/usr/bin/env perl
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
############# GetOptions #################
my ($fasta, $key, $outdir, $bam, $whitelist, $blacklist);
my ($mpercent, $mcov, $covfile, $dumpfile);
my ($samtools, $threads);
my ($kmer, $mode);

GetOptions(
	"h|?"=>\&help,
	"f:s"=>\$fasta,
	"k:s"=>\$key,
	"o:s"=>\$outdir,
	"bam:s"=>\$bam,
	"dum:s"=>\$dumpfile,
	"cov:s"=>\$covfile,
	"wtl:s"=>\$whitelist,
	"bll:s"=>\$blacklist,
	"mpr:s"=>\$mpercent,
	"mcv:s"=>\$mcov,
	"kmer:s"=>\$kmer,
	"mode:s"=>\$mode,
	"s:s"=>\$samtools,
	"t:s"=>\$threads,
) || &help;
&help unless ($fasta && $key && $outdir && ($bam || $dumpfile));

sub help
{
	print"
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
	-mcv  <int>    min k-mer coverage(\%)       [30]
	-kmer <int>    k-mer size                  [17]
	-mode <1/2>    1:ratio only; 2:ratio * cov [1]
	-s    <bin>    samtools binary             [/Bio/bin/samtools-1.14]
	-t    <int>    threads for samtools view   [4]
	-h           Help document
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$samtools ||= "/Bio/bin/samtools-1.14";
$mpercent //= 0.3;
$mcov //= 30;
$threads ||= 4;
$kmer ||= 17;
$mode ||= 1;

$fasta = &abs_dir($fasta);

if($bam){
	$bam = &abs_dir($bam);
}
if($dumpfile){
	$dumpfile = &abs_dir($dumpfile);
}

if($whitelist){
	$whitelist = &abs_dir($whitelist);
}

if($covfile){
	$covfile = &abs_dir($covfile);
}

if($blacklist){
	$blacklist = &abs_dir($blacklist);
}

if(!-d $outdir){
	mkdir $outdir;
}
$outdir = &abs_dir($outdir);

my $dumpF = "$outdir/$key\.dump.hash";

##### read white list
my %white;
if($whitelist){
	open(IN, $whitelist);
	while (<IN>) {
		chomp;
		my ($id, $other) = split(/\s+/, $_, 2);
		$white{$id} = 0;
	}
	close(IN);
}

##### read black list
my %black;
if($blacklist){
	open(IN, $blacklist);
	while (<IN>) {
		chomp;
		my ($id, $other) = split(/\s+/, $_, 2);
		$black{$id} = 0;
	}
	close(IN);
}

my %fcov;
if($covfile){
	open(IN, $covfile);
	while (<IN>) {
		chomp;
		next if(/^#/);
		my @line = split(/\s+/, $_);
		$fcov{$line[0]} = $line[4];
	}
	close(IN);
}

##### get fasta length
my %flen;
my $fin;
open($fin, $fasta);
my $in = Bio::SeqIO->new(
	-fh     =>  $fin,
	-format =>  'Fasta',
);
while (my $chunk = $in->next_seq) {
	my $seqid = $chunk->id;
	my $seq = $chunk->seq;
	next if(!defined $seq || length($seq) < 1);
	$flen{$seqid} = length($seq);
}
close($fin);

##### read bam
my $readN = 0;
my $fread = "";
my %share;
my %count;
my %chrname;
my @seqorder = sort{$flen{$b} <=> $flen{$a}} keys %flen;

foreach my $c (@seqorder) {
	if(!exists $fcov{$c}){
		$fcov{$c} = sprintf("%.03f", $count{$c} * $kmer/$flen{$c} * 100);
	}
}

if(!$dumpfile || !-e $dumpfile){
	open(my $lgfh, ">$outdir/$key\.bam.log");
	open(IN, "$samtools view -\@ $threads $bam|");
	while (<IN>) {
		my ($read, $idn, $contig, $other) = split(/\t/, $_, 4);
		if($readN == 0){
			$chrname{$contig} = 0;
			$fread = $read;
			$readN++;
			next;
		}
		if($read ne $fread){
			my @tmporder = sort{$flen{$b} <=> $flen{$a}} keys %chrname;
			foreach my $c1 (@tmporder) {
				$count{$c1}++;
				foreach my $c2 (@tmporder) {
					next if($c1 eq $c2);
					last if($flen{$c1} > $flen{$c2});
					$share{$c1}{$c2}++;
				}
			}
			undef %chrname;
		}
		$chrname{$contig} = 0;
		$fread = $read;
		$readN++;
		my $yu = $readN % 1000000;
		if($yu == 0){
			my $pn = $readN/1000000;
			$current_T = &date_format(localtime());
			print $lgfh "$current_T  $pn M bam line processed\n";
			$lgfh->flush();
			print "$current_T  $pn M bam line processed\n";
		}
	}
	close(IN);
	close($lgfh);
	## for the last kmer
	my @tmporder = sort{$flen{$b} <=> $flen{$a}} keys %chrname;
	foreach my $c1 (@tmporder) {
		$count{$c1}++;
		foreach my $c2 (@tmporder) {
			next if($c1 eq $c2);
			last if($flen{$c1} > $flen{$c2});
			$share{$c1}{$c2}++;
		}
	}
	undef %chrname;
	my %compact;
	%{$compact{share}} = %share;
	%{$compact{count}} = %count;
	store \%compact, $dumpF;
	undef %compact;
}else{
	my %compact = %{retrieve($dumpfile)};
	%share = %{$compact{share}};
	%count = %{$compact{count}};
}

##
my %delete;
foreach my $c1 (@seqorder) {
	next if(exists $delete{$c1});
	foreach my $c2 (@seqorder) {
		next if(exists $delete{$c2});
		next if($c1 eq $c2);
		next if($flen{$c2} > $flen{$c1});
		my $c2cnt = $count{$c2};
		my $shcnt = $share{$c2}{$c1};
		if(!defined $shcnt){
			$shcnt = 0;
		}
		if($c2cnt && $c2cnt > 0){
			my $frac = sprintf("%.6f", $shcnt/$c2cnt);
			my $cfrac = sprintf("%.6f", $frac * $fcov{$c2} / 100);
			if($frac > 1){
				print "$c1\t$c2\t$c2cnt\t$shcnt\n";
			}
			if($mode == 1){
				if($frac >= $mpercent){
					$delete{$c2}{to} = $c1;
					$delete{$c2}{sim} = $frac;
					$delete{$c2}{csim} = $cfrac;
				}
			}else{
				if($cfrac >= $mpercent){
					$delete{$c2}{to} = $c1;
					$delete{$c2}{sim} = $frac;
					$delete{$c2}{csim} = $cfrac;
				}
			}
		}else{
			$delete{$c2}{to} = $c1;
			$delete{$c2}{sim} = -1;
			$delete{$c2}{csim} = -1;
		}
	}
}

my %keep;
my $seqn = 0;
foreach my $c1 (@seqorder) {
	if(exists $keep{$c1}){
		next if($keep{$c1}{keep} eq "abandon");
	}
	if(exists $delete{$c1}){
		$keep{$c1}{to} = $delete{$c1}{to};
		$keep{$c1}{sim} = $delete{$c1}{sim};
		$keep{$c1}{csim} = $delete{$c1}{csim};
		$keep{$c1}{keep} = "abandon";
	}else{
		if($seqn == 0){
			$keep{$c1}{to} = "null";
			$keep{$c1}{sim} = 0;
			$keep{$c1}{csim} = 0;
			$keep{$c1}{keep} = "keep";
		}else{
			my %ftmp;
			my %ctmp;
			my $flag = 0;
			foreach my $c0 (@seqorder) {
				next if($c0 eq $c1);
				next if($flen{$c1} > $flen{$c0});
				if(exists $keep{$c0}){
					next if($keep{$c0}{keep} eq "abandon");
				}
				my $c1cnt = $count{$c1};
				my $shcnt = $share{$c1}{$c0};
				if(!defined $shcnt){
					$shcnt = 0;
				}
				if($c1cnt && $c1cnt > 0){
					my $frac = sprintf("%.6f", $shcnt/$c1cnt);
					$ftmp{$c0} = $frac;
					$ctmp{$c0} = sprintf("%.6f", $frac * $fcov{$c1}/100);
					$flag++;
				}
			}
			if($flag == 0){
				$keep{$c1}{to} = "null";
				$keep{$c1}{sim} = -1;
				$keep{$c1}{csim} = -1;
				$keep{$c1}{keep} = "abandon";
			}else{
				if($mode == 1){
					my @lseq = sort{$ftmp{$b} <=> $ftmp{$a}} keys %ftmp;
					$keep{$c1}{to} = $lseq[0];
					$keep{$c1}{sim} = $ftmp{$lseq[0]};
					$keep{$c1}{csim} = $ctmp{$lseq[0]};
					if($ftmp{$lseq[0]} >= $mpercent){
						$keep{$c1}{keep} = "abandon";
					}else{
						$keep{$c1}{keep} = "keep";
					}
				}else{
					my @lseq = sort{$ctmp{$b} <=> $ctmp{$a}} keys %ctmp;
					$keep{$c1}{to} = $lseq[0];
					$keep{$c1}{sim} = $ftmp{$lseq[0]};
					$keep{$c1}{csim} = $ctmp{$lseq[0]};
					if($ctmp{$lseq[0]} >= $mpercent){
						$keep{$c1}{keep} = "abandon";
					}else{
						$keep{$c1}{keep} = "keep";
					}
				}
			}
		}
	}
	$seqn++;
}

foreach my $c (@seqorder) {
	next if($keep{$c}{keep} eq "abandon" && $keep{$c}{sim} > 0);
	if($keep{$c}{sim} < 0 || $fcov{$c} < $mcov){
		$keep{$c}{keep} = "lowcov";
	}
}

foreach my $c (@seqorder) {
	if(exists $black{$c}){
		$keep{$c}{keep} = "black";
	}
}

foreach my $c (@seqorder) {
	if(exists $white{$c}){
		$keep{$c}{keep} = "keep";
	}
	if(exists $count{$c}){
		$keep{$c}{count} = $count{$c};
	}else{
		$keep{$c}{count} = 0;
	}
}

open(OUT, ">$outdir/$key.all.xls");
open(KP, ">$outdir/$key.keep.xls");
print OUT "#contig\tlength\tmapKmer\tcoverage\tcompareTo\tshareRatio\tshareCov\tfate\n";
print KP "#contig\tlength\tmapKmer\tcoverage\tcompareTo\tshareRatio\tshareCov\tfate\n";
foreach my $c (@seqorder) {
	print OUT "$c\t$flen{$c}\t$keep{$c}{count}\t$fcov{$c}\t$keep{$c}{to}\t$keep{$c}{sim}\t$keep{$c}{csim}\t$keep{$c}{keep}\n";
	if($keep{$c}{keep} eq "keep"){
		print KP "$c\t$flen{$c}\t$keep{$c}{count}\t$fcov{$c}\t$keep{$c}{to}\t$keep{$c}{sim}\t$keep{$c}{csim}\t$keep{$c}{keep}\n";
	}
}
close(OUT);
close(KP);

open(FA, ">$outdir/$key.keep.fa");
open(DROP, ">$outdir/$key.drop.fa");
open($fin, $fasta);
my $in2 = Bio::SeqIO->new(
	-fh     =>  $fin,
	-format =>  'Fasta',
);
while (my $chunk = $in2->next_seq) {
	my $seqid = $chunk->id;
	my $seq = $chunk->seq;
	next if(!defined $seq || length($seq) < 1);
	&formatSeq(\$seq);
	if($keep{$seqid}{keep} eq "keep"){
		print FA ">$seqid\n$seq\n";
	}else{
		print DROP ">$seqid\n$seq\n";
	}
}
close($fin);
close(FA);
close(DROP);

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&Runtime($begin_time);
##########################################
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
