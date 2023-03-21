#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Data::Dumper;
use List::Util qw{any first};
use Bio::SeqIO;
############# GetOptions #################
my ($fasta, $outfile, $clean, $rename, $keep);
my ($keepj, $minlen, $replace, $upcase);
GetOptions(
	"h|?"=>\&help,
	"f:s"=>\$fasta,
	"o:s"=>\$outfile,
	"c:s"=>\$clean,
	"n:s"=>\$rename,
	"u:s"=>\$upcase,
	"k:s"=>\$keep,
	"j:s"=>\$keepj,
	"l:s"=>\$minlen,
	"r:s"=>\$replace,
) || &help;
&help unless ($fasta && $outfile);

sub help
{
	print"
	Description: format fasta sequence

	-f  <file>  fasta sequence                       [force]
	-o  <file>  output fasta sequence                [force]
	-c  <T/F>   clean(remove sequence with * or .)   [T]
	-n  <T/F>   rename gene ids to avoid duplication [T]
	-u  <T/F>   up case all sequences                [T]
	-k  <T/F>   keep comment in head or not          [F]
	-j  <T/F>   keep J aa or not                     [F]
	-l  <int>   min sequence length                  [1]
	-r  <T/F>   replace | in sequence or not         [F]
	-h          Help document
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$fasta = &abs_dir($fasta);
$clean ||= "T";
$rename ||= "T";
$upcase ||= "T";
$keep ||= "F";
$keepj ||= "F";
$minlen ||= 1;
$replace ||= "F";
my $fin;
open($fin, $fasta);
my $in = Bio::SeqIO->new(
	-fh     =>  $fin,
	-format =>  'Fasta',
);
open(OUT, ">$outfile");
my $i = 1;
while (my $chunk = $in->next_seq) {
	my $seqid = $chunk->id;
	if($replace eq "T" && $seqid =~ /\|/){
		$seqid =~ s/\|/_/g;
	}
	if($rename eq "T"){
		$seqid .= "\_WXYZ\_" . sprintf("%08d", $i);
	}
	if($keep eq "T"){
		$seqid .= ' ' . $chunk->desc;
	}
	my $seq = $chunk->seq;
	next if(!defined $seq || length($seq) < 1);
	if(!defined $seq){
		print "ERROR: no seq : $seqid\n";
		die;
	}
	if($upcase eq "T"){
		$seq = uc($seq);
	}
	next if($seq =~ /J/ && $keepj ne "T");
	next if(length($seq) < $minlen);
	if($seq =~ /\./){
		$seq =~ s/\./\*/g;
	}
	if($seq =~ /\-/){
		$seq =~ s/\-/\*/g;
	}
	$seq =~ s/\*$//g;
	next if($seq =~ /\*/ && $clean eq "T");
	&formatSeq(\$seq);
	print OUT ">$seqid\n$seq\n";
	$i++;
}
close($fin);
close(OUT);
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
