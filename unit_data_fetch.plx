#this is for extracting the data of certain unit from the unit_combine file.
use strict;
use warnings;

print "enter unit locus_tag file:\t";
my $locus_file=<>;

print "enter target file:\t";
my $data_file=<>;

chomp ($data_file,$locus_file);

#outfile name
my @out_file=split/\./,$locus_file;
my $outf=$out_file[0].".da";

#readin locus tag
my @locus_tag;
my $i=-1;
open IN, "$locus_file" or die $!;
while(<IN>)
  {
  	chomp;
  	$locus_tag[++$i]=$1 if($_=~/(\w+)/);
  }
close IN;

#read the data file and fetch the infomation of the locus tags

open OUT, ">$outf" or die $!;
open IN, "$data_file" or die $!;
while(<IN>)
  {
  	chomp;
  	my $tmp=$_;
  	foreach $i (0..$i)
  	  {
  	    if($tmp=~/$locus_tag[$i]/)
  	      {
  	        print OUT "$tmp\n";
  	      }
  	  }
  }
close IN;
close OUT;
