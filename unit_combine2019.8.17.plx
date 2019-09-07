#this is for dealing with the site combine file, gather the unit together.
#updated 2019.08.17: add the key site_start as there may be same relative location of different region
#        *****       there will be repeat for one iteration of hash combine. use hash combine for a second time in this script.
#updated 2019.08.16: retain the reads count and site count.
#updated 2019.07.27: fix the bug that there may be mutiple out for same unit or same site_easy 
#updated 2019.07.25: add bla_count to the result 
#updated 2019.07.24: fix the length calculation bug, add 1 to the difference
#                    add unit type in the genome unit info result
#                    change the location of the unit length to the second column
#updated 2019.7.19: replace bla_count with site_counts. This change provides the inserted site number of the region and remove the different records number matched to genome(same read sequence is taken as the same record).
use strict;
use warnings;

print "enter site combine file:\t";
my $file=<>;

chomp $file;
my @out_file=split/site_combine/,$file;
my $outf=$out_file[0]."unit_combine.re";

print "enter after/before length:\t";
my $region=<>;
chomp $region;

#######################################################################genome info
my %feature;
my %feature_n;
my $tmpl;#last site
my @tmplf;#last info
my $i=-1;
my $j=-1;
open IN,"feature.db" or die $!;
while(<IN>)
  {
  	chomp;
  	#strain	         organelle	genome	      source	function_unit_type	function_unit_id	chromosome	function_unit site_easy	function_unit site_complex	function_unit_function	function_unit_note	
    #Yarrowia lipolytica	na	CP017553.1		1..2257857	protein	YALI1_A00014g	1A	1230..1648	1230..1457,1523..1648	hypothetical protein	AOW00072.1	na
    #...........1..........2.......3.........4..........5.......6..............7.......8................9...................10...............11........12
    #...........1........2........3..........4.........5........6..........7.........8.........9.........10........11.....
    if($_=~/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+\.\.\d+)\t([^\t]+)\t([^\t]+)\t([^\n]+)/)
      {
      	#get the start and end of region
      	my @tmpn=split /\.\./,$8;#start-end
        #my $tmp=substr $3,0,8; #chromosome
        my $tmp=$3;
        #the first region
      	if(!($j+1))#initilize 
      	  {
      	  	$j++; 
      	  	$i++; #id,new record, blank before the first region
      	  	#the blank before the first region 
      	  	#chromosome num unit_name site_easy site_complex unit source
            unit_ass($tmp,$i,"na","1\.\.".($tmpn[0]-1),"na","space",$4,"space");        
           
            $i++;   #id, new record, the second one, the first region
            #the region information from the first region start to the region end
            unit_ass($tmp,$i,$6,$8,$9,"unit",$4,$5);
            $tmpl=$tmpn[1]+1;#store the value+1 of the region end (for the start of the next region)
            
            #store the information of the current region(for the last one)
            @tmplf=($tmp,$i,$6,$8,$9,"unit",$4,$5);
            next;    	    
      	  }
      	else
      	  {
      	    if($feature{$tmp}->[0][0])#same chromosome
      	      {
      	        $i++; #id,new record
                #the blank region infomation of the nst region
      	  	    #chromosome num unit_name site_easy site_complex location source
                unit_ass($tmp,$i,"na","$tmpl\.\.".($tmpn[0]-1),"na","space",$4,"space");            
                
                $i++;
                #the region information from the nst region start to the region end
                unit_ass($tmp,$i,$6,$8,$9,"unit",$4,$5);
                $tmpl=$tmpn[1]+1;#store the value+1 of the region end (for the start of the next region)
                #store the information of the current region(for the last one)
                @tmplf=($tmp,$i,$6,$8,$9,"unit",$4,$5);
                next;
      	      }
      	    elsif(!$feature{$tmp}->[0][0])#different chromosome
              {
              	#get the information of the last one
              	my @tmpv=split /\.\./,$tmplf[6];
              	unit_ass($tmplf[0],$i+1,"na","$tmpl\.\.$tmpv[1]","na","space",$tmplf[6],"space");#last chromosome the last one 
              	
                $i=0;#new chromosome
                unit_ass($tmp,$i,"na","1\.\.".($tmpn[0]-1),"na","space",$4,"space");#new first blank for the new chromosome
               
                #first region
                $i++;
                unit_ass($tmp,$i,$6,$8,$9,"unit",$4,$5);
                
                $tmpl=$tmpn[1]+1;#store last site
                @tmplf=($tmp,$i,$6,$8,$9,"unit",$4,$5);
                next;                	
              }
            else
              {
                print "error:$_\n";	
              } 
          }      	
      }	
  }
close IN;
my @tmpv=split /\.\./,$tmplf[6];#chr length
unit_ass($tmplf[0],$i+1,"na","$tmpl\.\.$tmpv[1]","na","space",$tmplf[6],"space");#last chromosome the last one

sub unit_ass
  {
  	my $subi=shift;#chromosome
  	my $subj=shift;#num
    $feature{$subi}->[$subj][0]=shift;#unit name
    $feature{$subi}->[$subj][1]=shift;#site_easy
    $feature{$subi}->[$subj][2]=shift;#site_complex
    $feature{$subi}->[$subj][3]=shift;#unit
    $feature{$subi}->[$subj][4]=shift;#source    
    $feature{$subi}->[$subj][5]=shift;#unit_type, added 2019.07.24
    $feature_n{$subi}=$subj;#number of features(including blank)
  }

open OUT, ">genome_unit_info.log" or die $!;
foreach my $key(keys %feature)
  {
    foreach $i(0..$#{$feature{$key}})	
      {
      	#print "$i\n";
      	#my $test=<>;
      	my @tmp=split /\.\./,$feature{$key}->[$i][1];
        print OUT "$key\t$feature{$key}->[$i][4]\t",$i+1,"\t$feature{$key}->[$i][0]\t$feature{$key}->[$i][1]\t",$tmp[1]-$tmp[0]+1,"\t$feature{$key}->[$i][2]\t$feature{$key}->[$i][3]\t$feature{$key}->[$i][5]\n";	
      }
  }
close OUT;



##############################################################################site_info
#read site_combine file
my %seqbs;#before sort
my $na=0;
my @re;
open IN, "$file" or die $!;
while(<IN>)
  {
  	chomp;
     #site_id	site	read_counts	bla_counts	genome	relative location	strain	organelle	function_unit_type	function_unit_id	chromosome	function_unit site_easy	function_unit site_complex	function_unit_function	function_unit_note
      #1 	1561041   	225        	38	   CP017553        	in	Yarrowia lipolytica		        protein	           YALI1_A15624g	   1A	        1560526..1562445	                      	          hypothetical protein	Compare to YALI0A15576g, similar to uniprot|Q07824 Saccharomyces cerevisiae YLL028w TPO1 polyamine transport	
  	#...1.......2........3...........4........5.........   6          ....7......8............9.....................10...........11...............12..................13........................14...........................15..................16.............).
  	#..........1.......2.........3.........4.........5......6..7..........8.........9........10.........11....12..13(12)......14
  	if($_=~/^(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\n]+)/)
  	  {
  	  	#get the value of site_combine file
  	    my @tmp=split /\t/,$_;
  	    my $val=$5;
  	    my $len;
  	    #print "$tmp[7]\t$tmp[8]\t$tmp[9]\n";
  	    #my $test=<>;
  	    my @fulsite=split /\,/,$tmp[12];
  	    my @esite=split /\.\./,$tmp[11];  
  	    	    
  	    #get the length of the region
  	    if($tmp[5] eq "in" or $tmp[5] eq "other_in")
      	  {
      	  	#site,#chromosome
      	  	@re=lookup($tmp[1],$tmp[4]);
      	  	$len=$re[0];
      	  	#print "$len\n";
      	  }
      	elsif($tmp[5] eq "out" or $tmp[5] eq "other_out")
      	  {
      	  	@re=lookup($tmp[1],$tmp[4]);
      	  	$len=$re[0]-2*$region;#exclude the length of prior and after
      	  	$val=~s/$tmp[11]/$re[1]/;#replace the unit region with the out region
      	  	#$val=~s/$tmp[12]/\t/;  	  	
      	  } 
  	    elsif($tmp[5] eq "intron" or $tmp[5] eq "exon")
  	      {
  	      	my $flen=$esite[1]-$esite[0];
  	      	my $ex=0;
  	        foreach my $fuli(0..@fulsite-1)
  	          {
  	          	my @fulsitee=split /\.\./,$fulsite[$fuli];
  	          	$ex+=$fulsitee[1]-$fulsitee[0];
  	          }
  	        if($tmp[5] eq "intron")
  	          {
  	            $len=$flen-$ex;
  	          }  
  	        else
  	          {
  	            $len=$ex;	
  	          }	
  	      }
  	    elsif($tmp[5]=~/prior|after/)
  	      {  	      	
  	      	@re=lookup($tmp[1],$tmp[4]);
  	      	$len=$re[0];
  	      	if($len>$region)
  	      	  {
  	      	    $len=$region;	
  	      	  }
  	      }
  	    
  	    #assign result to seqbs%  
  	    #....gme....unit...loc.....rc.......bc......len....other
  	    #no locus_tag, not out relative location
  	    if($tmp[9]!~/\w/ and $tmp[5] ne "out")
  	      {
  	      	#print "ok\n";
  	        assign($tmp[4],"na".$na++,$tmp[5],$esite[0],$tmp[2],$len,$val); 	  	        
  	      }
  	    #no locus_tag, out
  	    elsif($tmp[9]!~/\w/ and $tmp[5] eq "out")
  	      {
  	      	#print "ok\n";
  	        assign($tmp[4],"na".$na++,$tmp[5],$esite[0],$tmp[2],$len,$val); 	  	        
  	      }
  	    #locus_tag, out
  	    elsif($tmp[9]=~/\w/ and $tmp[5] eq "out")
  	      {
  	        assign($tmp[4],$tmp[9],$tmp[5],$esite[0],$tmp[2],$len,$val); 	
  	      }  
  	    #locus_tag, not out
  	    elsif($tmp[9]=~/\w/ and $tmp[5] ne "out")
  	      {
  	        assign($tmp[4],$tmp[9],$tmp[5],$esite[0],$tmp[2],$len,$val); 	
  	      }
  	  }
  	elsif($_=~/site_id/)
  	  {
  	    next;	
  	  }
  	else
  	  {
  	    print "error\t$_\n";	
  	  }
  }
close IN;

sub assign
  {
    my $gme=shift;	
    my $unit=shift;
    my $loc=shift;
    my $sta_end=shift;
    $seqbs{$gme}->{$unit}->{$loc}->{$sta_end}->[0]+=shift;#read_count
    #$seqbs{$gme}->{$unit}->{$loc}->{$sta_end}->[4]+=shift;#bla_count
    $seqbs{$gme}->{$unit}->{$loc}->{$sta_end}->[1]+=1;#inserted site number of the region, added 2019.7.19
    $seqbs{$gme}->{$unit}->{$loc}->{$sta_end}->[2]=shift;#length
    $seqbs{$gme}->{$unit}->{$loc}->{$sta_end}->[3]=shift;#other information 
    
    #print "$seqbs{$gme}->{$site}->[0]\t$seqbs{$gme}->{$site}->[1]\t$seqbs{$gme}->{$site}->[3]\n"
  }

sub lookup #match genome and site info;
  {
  	my $site=shift;#site
  	my $gme=shift;#genome 
  	my $min=0;
  	my $max=$feature_n{$gme};
  	my $subi=0;
  	my @tmp;
  	my $match;
  	my $si=0;
  	while(!$match)
  	  {
  	  	$si++;
  	  	#$feature{$subi}->[$subj][1]=shift;#site_easy
  	  	$subi=int(($min+$max)/2);            	
        @tmp=split /\.\./,$feature{$gme}->[$subi][1] if($feature{$gme}->[$subi][1]);
  	    if($site <= $tmp[1] and $site >= $tmp[0])
  	      {
  	      	#print "$gme\t$site\t$tmp[1]\t$tmp[0]\n";
  	      	#length,#annotation
  	      	return ($tmp[1]-$tmp[0]+1,$feature{$gme}->[$subi][1]);
  	      	$match=1;
  	      	last;
  	      }	
  	    elsif($site > $tmp[1])
  	      {
  	      	$min=$subi;
  	      }
  	    elsif($site < $tmp[0])
  	      {
  	      	$max=$subi;
  	      }
  	    if($si >= "30")
  	      {
  	        #print "$si\t$site\t$subi\t$tmp[0]\t$tmp[1]\t$max\t$min\n";
  	        $min=$max;	
  	        #my $test=<>;	
  	      } 
  	  }
  }

open OUT, ">$outf" or die $!;
my $num=0;
print OUT"id\tfrag_length\tread count\tsite_count\tgenome\trelative location\tstrain\torganelle\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tfunction_unit_function\tdirection\tfunction_unit_note\n";
foreach my $m(sort{$a cmp $b}keys %seqbs)
  {
    foreach my $n(sort keys %{$seqbs{$m}})
      {
      	foreach my $k(sort keys %{$seqbs{$m}->{$n}})
         {
         	 foreach my $o(sort keys %{$seqbs{$m}->{$n}->{$k}})
         	   {
      	       $num++;
      	       print "$m\t$n\t$k\t" if(!$seqbs{$m}->{$n}->{$k}->{$o}->[2]);
      	       print OUT "$num\t$seqbs{$m}->{$n}->{$k}->{$o}->[2]\t$seqbs{$m}->{$n}->{$k}->{$o}->[0]\t$seqbs{$m}->{$n}->{$k}->{$o}->[1]\t$seqbs{$m}->{$n}->{$k}->{$o}->[3]\n";        	
      	     }
      	 }
      }	
  }
close OUT;

#deal unit_combine file for a second time
{
  my %seqbsa;#before sort
  my $na=0;
  my @re;
  my $numt;
  open IN, "$outf" or die $!;
  while(<IN>)
    {
    	chomp;
      if($_=~/\d+\t\d+\t\d+\t\d+\t([^\n]+)/)
        {
        	#get the value
        	my @tmp=split /\t/,$_;
        	my $gme=$tmp[4];#chr information
        	my $unit=$tmp[9];#locus tag
        	my $loc=$tmp[5];#relative location
        	my $region=$tmp[11];#site easy
        	print "$_\n" if($seqbs{$gme}->{$unit}->{$loc}->{$region}->[0]);
        	$seqbsa{$gme}->{$unit}->{$loc}->{$region}->[0]=$tmp[1]; #unit length
        	$seqbsa{$gme}->{$unit}->{$loc}->{$region}->[1]+=$tmp[2]; #read count
        	$seqbsa{$gme}->{$unit}->{$loc}->{$region}->[2]+=$tmp[3]; #site count
        	$seqbsa{$gme}->{$unit}->{$loc}->{$region}->[3]=$1;#annotation information      	
        	$numt++;
        }
     }
  close IN;
  
  open OUT, ">$outf" or die $!;
  my $num=0;
  print OUT"id\tfrag_length\tread count\tsite_count\tgenome\trelative location\tstrain\torganelle\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tfunction_unit_function\tdirection\tfunction_unit_note\n";
  foreach my $m(sort{$a cmp $b}keys %seqbsa)
    {
      foreach my $n(sort keys %{$seqbsa{$m}})
        {
        	foreach my $k(sort keys %{$seqbsa{$m}->{$n}})
           {
           	 foreach my $o(sort keys %{$seqbsa{$m}->{$n}->{$k}})
           	   {
           	   	 $num++;
        	       print OUT "$num\t$seqbsa{$m}->{$n}->{$k}->{$o}->[0]\t$seqbsa{$m}->{$n}->{$k}->{$o}->[1]\t$seqbsa{$m}->{$n}->{$k}->{$o}->[2]\t$seqbsa{$m}->{$n}->{$k}->{$o}->[3]\n";        	
        	     }
        	 }
        }	
    }
  close OUT;
}