#this is for counting the no inserted genes in all samples.

use strict;
use warnings;

print "enter the minimal sample number of the no_inserted gene:\n";
my $val=<>;

chomp $val;

my %id;#store the id of no inserted one

#read in all samples
opendir DIR,"." or die $!;
while(readdir DIR)
  {
    chomp;
    if($_=~/.no_insertion\.gene/)
      {
        open IN, "$_" or print $!;
        while(<IN>)
          {
          	chomp;
          	if($_=~/\w/)
          	  {
          	    my @tmp=split /\t/,$_;
          	    $id{$tmp[5]}++ if($id{$tmp[5]});
          	    $id{$tmp[5]}=1 if(!$id{$tmp[5]})	
          	  }
          }
        close IN;
       }
  }
close DIR;

#print out no_inserted gene
open OUT, ">no_inserted.name" or die $!;
foreach my $i(sort keys %id)
  {
    print OUT "$i\n" if($id{$i}>=$val);
  }
close OUT;