#This is for get the sum of read count number as well as counting the records number of the same genome site in the annotation file. 
use strict;
use warnings;

print "enter annotation.re file:\t";
my $file=<>;
chomp ($file);

my @filen=split /\_annotation.re/,$file;
my $outfile="$filen[0]"."_site_combine.re";

my %seqbs;#before sort
my $m=-1;
my $n=0;
open IN, "$file" or die $!;
while(<IN>)
  {
  	chomp;
  	#read_id...count...read_id_full	  genome	  identity	alignment length	mismatches	gap opens	 q. start	  q. end	 s. start	s. end	  evalue	bit score	relative location	    strain	          organelle	function_unit_type	function_unit_id	chromosome	function_unit site_easy	function_unit site_complex	function_unit_function	function_unit_note	
  	#1	       154252	   1|154252	   CP017555	    100	       94	                 0	        0	        14	      107	     47989	 47896  	4.47E-44	174	           out	    Yarrowia lipolytica		               protein	         YALI1_C00464g       	1C	           46482..47651		                                hypothetical protein	Compare to YALI0C00407g, Isopropyl Malate Dehydrogenase (LEU2),uniprot|P18120 Yarrowia lipolytica beta-isopropyl- malate dehydrogenase, similar to Saccharomyces cerevisiae LEU2 (YCL018W); ancestral locus
  	#1...........2.........3.............4.........5..............6..............7..........8..........9........10........11......12........13.......14............15...............16..................17..........18.....................19............20.................21......................22.........................23..................24.........
  	#............1.....2.........3.........4.........5.........6.........7.........8.........9........10.........11........12.......13........14........15.........16......17......18.......19........20........21........22.........23....
    if($_=~/^(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\n]+)/)
      {  
      	#chr,site,read count,relative location, other information                        
         assign($4,$11,$2,$15,$16);
      }
    elsif($_!~/\d/)
      {
        next;	
      }
    else
      {
        print "error $_\n";	
      }
    #my $test=<>;
  }
close IN;

sub assign
  {
    my $gme=shift;	
    my $site=shift;
    $seqbs{$gme}->{$site}->[0]+=shift;#read_count
    $seqbs{$gme}->{$site}->[1]=shift;#relative location
    $seqbs{$gme}->{$site}->[2]=shift;#other information
    $seqbs{$gme}->{$site}->[3]++;#bla_count 
    #print "$seqbs{$gme}->{$site}->[0]\t$seqbs{$gme}->{$site}->[1]\t$seqbs{$gme}->{$site}->[3]\n"
  }

open OUT, ">$outfile" or die $!;
my $num=0;
print OUT"site_id\tsite\tread_counts\tbla_counts\tgenome\trelative location\tstrain\torganelle\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tfunction_unit_function\tunit_direction\tabnormal\tfunction_unit_note\n";
foreach $m(sort{$a cmp $b}keys %seqbs)
  {
    foreach $n(keys %{$seqbs{$m}})
      {
      	$num++;
      	print OUT "$num\t$n\t$seqbs{$m}->{$n}->[0]\t$seqbs{$m}->{$n}->[3]\t$m\t$seqbs{$m}->{$n}->[1]\t$seqbs{$m}->{$n}->[2]\n";        	
      }	
  }
close OUT;