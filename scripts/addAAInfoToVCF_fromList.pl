#!/usr/bin/perl

# addAAInfoToVCF_fromList.pl
# written by Linnéa Smeds                    8 Nov 2016
# =====================================================
# Takes a vcf file and a list with ancestral
# alleles, and adds this to the vcf INFO column for
# each position.
#
# The ancestral allele list should have the format
# CHR:POS ANC_ALLELE
# =====================================================
# usage: perl addAAInfoToVCF_fromList.pl file.vcf anc.list new.vcf

use strict;
use warnings;

# INPUT PARAMETERS
my $VCFFILE = $ARGV[0];	# The vcf file with all positions
my $LIST = $ARGV[1]; 	# List with ancestral alleles
my $OUT = $ARGV[2];		# The new vcf file

# Starting clock
my $time = time;

# OPEN FILE HANDLE TO OUTPUT FILE
open(OUT, ">$OUT");

# SAVE ANCESTRAL ALLELES IN HASH
my %hash=();
open(LIST, $LIST);
my $acnt=0;
while(<LIST>) {
	if($acnt==0) {
		unless($_ =~ m/\w+:\d+\s[ATCGN]/) {
			die "File with ancestral allele does not have the format CHR:POS ANC_ALLELE\n but looks like this: $_";
		}
	}
	my ($id,$anc)=split(/\s+/, $_);
	$hash{$id}=$anc;
	$acnt++;
}
print STDERR "$acnt ancestral alleles saved\n";

# Open FILE HANDLE TO INPUT FILE
my $cnt=0;
open(IN, $VCFFILE);
while(<IN>) {
	# Found the first INFO line
	if(/^##INFO/)	{
		print OUT $_;

		#NOTE! IF YOU ALREADY HAVE THE INFO AA LINE, COMMENT EVERYTHING FROM HERE TO */
		#Loop through all lines starting with ##INFO
		my $next = <IN>;
		while ($next =~ m/^##INFO/) {
			print OUT $next;
			$next=<IN>;
		}
		#Outside of the loop-> found a line that didn't start with INFO.
		#Print the extra info line (and the other line)
		print OUT "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
		print OUT $next;
		#*/
	}
	# Any other type of header line, just print
	elsif(/^#/) {
		print OUT $_;
	}
	# Anything that doesn't start with "#" -> a variant!
	else {

		my @tab=split(/\s+/, $_);
		my $scaf=$tab[0];
		my $pos=$tab[1];
		my $id=$scaf.":".$pos;
		#print "DEBUG: looking at $id\n";

		# Check if there is an ancestral allele
		if(exists $hash{$id}) {
			$tab[7]=$tab[7].";AA=".$hash{$id};
 			print OUT join("\t", @tab)."\n";
			$cnt++;
			#print "\t printed to output\n"; #DEBUG
		}
		# This part can be removed if other lines should be skipped, or
		# changed to print OUT $_; (if other lines should be printed)
		else {
			#die "Did not find ancestral allele for $id!!\n";
			print "$id: NO anc found!\n";

		}
	}
}
close(IN);
close(OUT);

print STDERR "Saved $cnt variants\n";
$time = time-$time;
if($time<60) {
	print STDERR "Total time elapsed: $time sec.\n";
}
else {
	$time = int($time/60 + 0.5);
	print STDERR "Total time elapsed: $time min.\n";
}
