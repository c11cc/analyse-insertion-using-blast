#annotate the bla result according to their relative location to features on the genome
#updated 2019.07.23: rewrite the part for reading genebank files
#                    add abnormal tags for the features that have other features inside it, write these features into abnormal.log
#                    assign the abnormal feature to the bla result only when the site is aligned to the exon part or the abnormal feature has no exon part.
#                    assign the feature inside the abnormal feature to the bla result for other circumstance. 

use strict;
use warnings;
#use Math::BigFloat;
#use List::AllUtils qw(min max);
#use Cwd;
#use Bio::Seq;
#use Bio::Tools::Run::StandAloneBlastPlus;
#use 5.012;

print "enter genebank file:\t";
my $file=<>;

print "enter blast file:\t";
my $file1=<>;

print "enter after/before length:\t";
my $region=<>;
chomp $region;

print "take unit distance less than $region all as prior?\tYes:1\tNo:2\n";
my $asbef=<>;
chomp $asbef;

print "enter the q_start site:\t";
my $qmin=<>;

print "enter the q_end site:\t";
my $qmax=<>;

print "enter the count threshold:\t";
my $cou=<>;

#use this part for dealing from the clean data input
#require "fastq_deal.plx";
#require "sorting.plx";
#
#
#opendir my $curr, "." or die $!;
#while (readdir $curr)
#  {
#  	chomp;
#  	if($_=~/\.fa/)
#  	  {
#        my $fac= Bio::Tools::Run::StandAloneBlastPlus->new(-db_name=>"./genome.db", -create=> "1",-db_data=>"./genomic.fna",-db_type=>"nucl");
#        my $result=$fac->blastn(-query=>"$_",
#                                  -outfile=>"1.bla",
#                                  -outformat=>6,
#                                  -num_thread=>3,
#                                  #-method_args=> [ -threshold=>1,
#                                  #                -evalue=> 0.05]
#                                 );# or print " write error";
#        $fac->cleanup;   
#      } 
#  }
#

my %feature;
my $i=-1;
my $j=-1;
my $notes=0;
my $joi=0;
my @num;
my $mrn=0;

open IN, "$file" or die $!;
while(<IN>)
  {
    chomp;
    my $tmp=$_;  
    if($tmp=~/FEATURES/)
      {	
        if($j gt -1 and $i gt -1)
          {
            $num[$j]=$i;
          } 
        $j++;#chrmosome
        $i=-1;#unit
        if($notes)
          {
            $notes=0;	
          } 
        if($joi)
          {
            $joi=0;	
          }                  
      }
    elsif($tmp=~/\/organism\=\"([^\t]+)\"/)#/organism="Yarrowia lipolytica"
      {
        $feature{organism}->[$j]=$1;      
        next;
      }      
    elsif($tmp=~/\/strain\=\"([^\t]+)\"/)#/strain="CLIB89(W29)"
      {
        $feature{strain}->[$j]=$1;   
        next;  
      }      
    elsif($tmp=~/\/mol_type\=\"([^\t]+)\"/)#/mol_type="genomic DNA"
      {
        $feature{mol_type}->[$j]=$1;  
        next;   
      } 
    elsif($tmp=~/\/organelle\=\"([^\t]+)\"/)#/organelle="mitochondrion"
      {
        $feature{organelle}->[$j]=$1; 
        next;    
      }
    elsif($tmp=~/\/chromosome\=\"([^\t]+)\"/)#/chromosome="1A"
      {
        $feature{chromosome}->[$j]=$1;     
        next;
      }
    elsif($tmp=~/VERSION[^\w]+([^\t]+)/)#VERSION     CP017553.1
      {
        $feature{chromosome_id}->[$j+1]=$1;     
        next;
      }           
    elsif($tmp=~/[^\d]+(\d+\.\.\d+)/)#gene/misc_feature/cds/rna            3285..3920  /   gene            <13383..14354
      {
      	my $val=$1;
      	if($tmp=~/source\D+(\d\.\.\d+)$/)#source          1..2257857
      	  {
      	  	$feature{source}->[$j]=$1;    
            next;
      	  }
      	elsif($tmp!~/[A-Z]|[a-z]/)#19174..19184,19490..19500,20671..20806,22104..22343,records with mutiple exon
      	  {      	  	                                                               
            if($joi eq 1 and $tmp=~/\D+(\d+[^\)]+)\)/)  #  # 315569..315894,315955..315975,316049..316064)                   
              {
                $feature{start}->[$j][$i][2].=$1; 
                $joi=0;       
              }  
                                                                      
            elsif($joi eq 1 and $tmp=~/\D+(\d+[^\)]+)/)  #      # 315569..315894,315955..315975,     
                                                                    #  316049..316064,          
              {
                $feature{start}->[$j][$i][2].=$1;
              } 
            else 
              {
                die "not considered line: $tmp";
              }
      	  }
      	elsif($tmp!~/source/)
      	  {

      	    #gene/misc_feature <3285..3920;gene/misc_feature 3285..3920;
      	    if($tmp=~/gene|misc_feature/) #new records
      	      { 
      	      	$i++;     
      	      	$feature{start}->[$j][$i][1]=$val;#region
      	      	$feature{start}->[$j][$i][0]=0;#default
      	      }    
      	    if($tmp=~/joi/)
      	      {
      	      	$feature{start}->[$j][$i][0]=1;# have intron 
      	      	
      	      	if($tmp=~/[^\d+](\d+[^\)]+)\)/)#oneline
      	      	  {
      	      	  	$feature{start}->[$j][$i][2]=$1;
      	      	  }
      	        elsif($tmp=~/join\((\d+[^\)]+)/)#mutiline
      	          {   
      	            $feature{start}->[$j][$i][2]=$1;  	      
      	          }
      	        else
      	          {
      	            die $tmp;	
      	          }
      	      }       	
  
      	    $feature{type}->[$j][$i]="other" if($tmp=~/misc_feature/); #misc_feature
      	    $feature{type}->[$j][$i]="RNA" if($tmp=~/RNA/);
      	    $feature{type}->[$j][$i]="protein" if($tmp=~/mRNA/);#replace the value of RNA when it is protein
      	    $feature{type}->[$j][$i]="exon" if($tmp=~/exon/);#exon
      	    $feature{type}->[$j][$i]="intron" if($tmp=~/intron/);#exon
      	    $feature{reverse}->[$j][$i]=1 if($tmp=~/complement/);
      	    
      	    if($notes)
              {
                $notes=0;	
              }  
            $joi=1 if($tmp!~/\)/);#more than one line
            $joi=0 if($tmp=~/\)/ and $joi);# one line  	;
          }
      }   
          
    elsif($tmp=~/mRNA\s+/)#mRNA            complement(join(1230..1457,1523..1648))
      {
      	$feature{type}->[$j][$i]="protein";#product type 
        if($notes)
          {
            $notes=0;	
          }  
        if($joi)
          {
            $joi=0;	
          }    
        $mrn=1;      	    	
      }
                 
    elsif($tmp=~/\/locus_tag\=\"([^\t]+)\"/)#/locus_tag="YALI1_A00014g"
      {
        $feature{name}->[$j][$i]=$1;
        if($notes)
          {
            $notes=0;	
          }  
        if($joi)
          {
            $joi=0;	
          }    
        $mrn=0;  
      }      
    elsif($tmp=~/\/product=\"([^\t]+)\"/)#/product="hypothetical protein"
      {
        $feature{product}->[$j][$i]=$1;
        if($notes)
          {
            $notes=0;	
          }
        if($joi)
          {
            $joi=0;	
          }
      }   
    elsif($tmp=~/\/protein_id\=\"([^\t]+)\"/)#/protein_id="AOW00072.1"
      {
        $feature{protein_id}->[$j][$i]=$1;
      }   
    elsif($tmp=~/\/gene\=\"([^\t]+)\"/)#/gene="trnA(tgc)"
      {
        $feature{protein_id}->[$j][$i]=$1;
      }  
    elsif($tmp=~/pseudo/)#/gene="trnA(tgc)"
      {
        $feature{type}->[$j][$i].="pseudo";
      }
    elsif($tmp=~/\/note\=\"([^\"]+)/)# /note="Compare to YALI0A00110g, weakly similar to
      {
        $feature{note}->[$j][$i]=$1;
        $notes=1;
        $mrn=0;
        next;
      } 
    
    elsif($tmp=~/ORIGIN/)
      {
        $notes=0;	
      }
    elsif($notes and $tmp=~/\/codon_start/)
      {
        $notes=0;        
      }   
            
    elsif($notes eq 1 and $tmp=~/                     [^\t]+\s\w+/ and $tmp!~/(\/locus_tag)|(\/codon_start)|(\/product)|(\/protein_id)|(\/translation)|(\/gene)|(\/inference)/)
      { 
      	if($tmp=~/                     ([^\"]+)$/)
      	  {
            $feature{note}->[$j][$i].=" ".$1;     
          }   
        elsif($tmp=~/                     ([^\t]+)\"/)
      	  {
            $feature{note}->[$j][$i].=" ".$1;     
          }
      }                                       
  } 
close IN;
$num[$j]=$i;#the last one 
my $chro_num=$j;

open (my $outlog, ">abnormal.log") or die $!;
#get the abnormal values whose region contain the latter or the third one(there will be errors if the abnormal one is out of this range)
my $n=0;
foreach $j(0..$j)
  {
  	$n=0;
    foreach $i(0..$num[$j])	
      {
      	if($i>=1 and $i<$num[$j])
      	  { 
      	    my @val_former=split /\.\./,$feature{start}->[$j][$i-1][1] if($feature{start}->[$j][$i-1][1]);#get the start and end of the former one
      	    my @val_current=split /\.\./,$feature{start}->[$j][$i][1] if($feature{start}->[$j][$i][1]);#get the start and end of the current one
      	    my @val_latter=split /\.\./,$feature{start}->[$j][$i+1][1] if($feature{start}->[$j][$i+1][1]);#get the start and end of the latter one      	
            if(
      	      ($val_current[0]>=$val_former[0] and $val_current[1] <=$val_former[1]) or   #store the value of the abnormal one bigger than the current
      	      ($val_latter[0]>=$val_former[0] and $val_latter[1] <=$val_former[1])   	#store the value of the abnormal one bigger than the latter      	  
      	      )
      	        {
      	    	    $feature{abnormal}->[$j][$n]=$i-1;
      	    	    $n++;
      	          print $outlog "n-1:$n-1\ti:$i\tj:$j\t$feature{abnormal}->[$j][$n-1]\t$feature{chromosome_id}->[$j]\t$feature{start}->[$j][$feature{abnormal}->[$j][$n-1]][1]\t$feature{start}->[$j][$i][1]\n"; 
      	        }
      	  }  
      	#bla or not,filehand,unit,chr
      	#prin_out(0,$outlog,$i,$j);
     }
  }
close $outlog;

###############################################################################################
#deal with bla file
my %bla;
open IN, "$file1" or die $!;
while(<IN>)
  {
  	chomp;
  	#query acc.ver, sject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
  	#        0|0          CP017555.1	100.000	     35	      0	      0	   115	  149	  46203	 46237	8.92e-011	 65.8
  	#........1..2.....3.........4.........5........6..........7......8.....9.....10.....11......12......13........14
    #if($_=~/^((\d+)\|(\d+))\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\n]+)$/)#sample d series
    if($_=~/^((\d+)\|(\d+)\|\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\n]+)$/)#
      { 
      	 if($9>=$qmin and $9<=$qmax)
      	   {
      	     if($5>95)
      	       {
      	       	 #my $val=Math::BigFloat->new($13);
      	       	 my $val=$13;
                 if(!$bla{$2}->[0])#new site
                   {     
                   	  {     
                        val_ass($2,$_,$3,$4,$5,$6,$11,$val);      	
                      }
                   }
                 else #same site, matain one record
                   {
                     if($bla{$2}->[6]>$val)#based on e_value
                       { 
                       	 val_ass($2,$_,$3,$4,$5,$6,$11,$val); 
                       }
                     elsif($bla{$2}->[6] eq $val and $bla{$2}->[4]<$6)#same e_value,based on alignment length
                       {
                       	 val_ass($2,$_,$3,$4,$5,$6,$11,$val);
                       }
                     elsif($bla{$2}->[3] < $5 and $bla{$2}->[4] eq $6 and $bla{$2}->[6] eq $val)#same e_value and alignment length, based on identity
                       {
                       	 val_ass($2,$_,$3,$4,$5,$6,$11,$val);
                       }
                   }
               }
           }
      }
  }
close IN;

sub val_ass
  {
  	my $stmp=shift;#id
  	$bla{$stmp}->[0]=shift;#full records    
    $bla{$stmp}->[1]=shift;#read_counts
    $bla{$stmp}->[2]=shift;#seq_id/chromosome
    $bla{$stmp}->[3]=shift;#identity
    $bla{$stmp}->[4]=shift;#aln_length
    $bla{$stmp}->[5]=shift;#s_start
    $bla{$stmp}->[6]=shift;#e_value 
  }

#####################################################################################
#get the insertion site and annotate the result
sub annotate
  {
  	#blast record
  	my $m=shift;
  	#normal or abnormal
  	my $abnormal=shift;
  	#my $sasbef=shift;
   	my @tmp1;
    my @tmp;
    my $match;#flag for matching result
    my $subj;#tmp variable for chrmosome
    my $subi;#the average of the upper and lower region id
    my $site=$bla{$m}->[5]; 	     
    $match=0;
    my $find=0;
    foreach $subj(0..$chro_num)#db
      {  	
      	if($feature{chromosome_id}->[$subj]=~/$bla{$m}->[2]/)#match chromosome
      	  { 
      	  	#reformat the chromosome format      	  	
      	  	if($bla{$m}->[2] ne $feature{chromosome_id}->[$subj])#the two chromosome records are stored differently 
      	  	  {
      	  	  	my @bla_tmp=split /\t/,$bla{$m}->[0];#separate blast result
      	  	  	$bla_tmp[1]=$feature{chromosome_id}->[$subj];#change the chromosome format of the blast result
                $bla{$m}->[0]=join("\t", @bla_tmp);#get the full blast records
      	  	  }
      	  	
      	  	#abnormal ones, foreach every one  
      	  	if($abnormal and $feature{abnormal}->[$subj][0])
      	  	  {
      	  	    foreach my $subn(0..@{$feature{abnormal}->[$subj]}-1)
      	  	      {     	  	      	
      	  	      	$subi=$feature{abnormal}->[$subj][$subn];#get the unit id
      	  	        die "error line:subn:$subn\tsubj:$subj\t$feature{abnormal}->[$subj][$subn]\t$feature{chromosome_id}->[$subj]\t$subj\t$subi\n" if(!$feature{start}->[$subj][$subi][1]);
      	  	        @tmp=split /\.\./,$feature{start}->[$subj][$subi][1];
      	  	        if($site <= $tmp[1] and $site >= $tmp[0])#inside the region
      	  	          {
      	  	          	my $other=0;#flag for judging "other"
                        if($feature{type}->[$subj][$subi])#with type value
                          {
                            ($match,$other)=others($m,$subj,$subi);            	    	 
                          }
                        if(!$other)#not "other"
                          {
              	  	        if($feature{start}->[$subj][$subi][0])
              	  	          {
              	  	          	#chr,unit,bla_record,abnormal
                                exon($subj,$subi,$m,1);
              	  	          }
              	  	        else 
              	  	          {  
                          	 	  assign($m,"in",$subi,$subj);
              	              }
              	          }
              	        $bla{$m}->[10]="abnormal_region";#set the abnormal value
      	  	          }
      	  	        else
      	  	          {
      	  	            next;	
      	  	          }	
      	  	      }	
      	  	  }
      	  	  
      	  	#non abnormal, bisection method
      	  	my $max=$num[$subj];#initial upper region id, the total records of certain chromosome
      	  	my $min=0;#initial lower region id
      	  	#print "matching...\n";  	 
      	  	my $iteration=0; #count the iteration times,for 1024, when $iteration=10, it will be 1		    	
            while(!$match)#unit
              {  
              	$iteration++;
              	$subi=int(($min+$max)/2);   
              	@tmp=split /\.\./,$feature{start}->[$subj][$subi][1];
              	if($site <= $tmp[1] and $site >= $tmp[0])# inside feature, value first, or there may be errors when the distance of two gene are smaller than $region.  	
              	  {              	  	
              	  	my $other=0;#flag for judging "other"
                    if($feature{type}->[$subj][$subi])#with type value
                      {
                        ($match,$other)=others($m,$subj,$subi);            	    	 
                      }
                    if(!$other)#not "other"
                      {
              	  	    if($feature{start}->[$subj][$subi][0])
              	  	      {
                            exon($subj,$subi,$m,0);
              	  	      }
              	  	    else 
              	  	      {  
                      	 	  assign($m,"in",$subi,$subj);
              	          }
              	        $match=1;	
              	      }
              	  } 
                elsif($site >= $tmp[1] and $site <= $tmp[1]+$region)# region after forward    	
                  {                  	  	
                  	my $other=0;
                    if($feature{type}->[$subj][$subi] and !$find)#others
                      {
                        ($match,$other)=others($m,$subj,$subi);
                      }
                    if(!$other and !$find)
                      {
                        assign($m,"prior $region"."bp",$subi,$subj) if($feature{reverse}->[$subj][$subi]);
                        assign($m,"after $region"."bp",$subi,$subj) if(!$feature{reverse}->[$subj][$subi]);
                        $find=1;
                        $min=$subi;
                      }
                    elsif(!$other and $find)
                      {
                      	if($asbef eq 1)
                      	  {
                            assign($m,"prior $region"."bp",$subi,$subj) if($feature{reverse}->[$subj][$subi]);
                          }
                        elsif($asbef eq 2)
                      	  {
                            assign($m,"after $region"."bp",$subi,$subj) if(!$feature{reverse}->[$subj][$subi]);
                          }
                        $find=0;
                        $match=1;	                          
                      }                    
                  }
                elsif($site <= $tmp[0] and $site >= $tmp[0]-$region)# 100bp prior forward    	
                  {                 	  	
                  	 my $other=0;
                    if($feature{type}->[$subj][$subi] and !$find)#others
                      {
                        ($match,$other)=others($m,$subj,$subi);
                      }
                    if(!$other and !$find)
                      {
                        assign($m,"after $region"."bp",$subi,$subj) if($feature{reverse}->[$subj][$subi]);
                        assign($m,"prior $region"."bp",$subi,$subj) if(!$feature{reverse}->[$subj][$subi]);
                        $find=1;
                        $max=$subi;
                      }
                    elsif(!$other and $find)
                      {
                        if($asbef eq 1)
                      	  {
                            assign($m,"prior $region"."bp",$subi,$subj) if(!$feature{reverse}->[$subj][$subi]);
                          }
                        elsif($asbef eq 2)
                      	  {
                            assign($m,"after $region"."bp",$subi,$subj) if($feature{reverse}->[$subj][$subi]);
                          }
                        $match=1;	
                        $find=1;
                      } 
                  }  	
              	elsif($site > $tmp[1]+$region)# bigger part     	
              	  {                 	  	
              	    $min=$subi;	
              	    if($subi eq int(($subi+$max)/2) and $iteration>20) 
              	      {
              	        if(!$bla{$m}->[7])
              	          {              	        	
              	  	        assign($m,"out",$subi,$subj);
              	  	      }
                	      $match=1;	
              	      }
              	  }      
              	elsif($site < $tmp[0]-$region)	#smaller part
              	  {
              	  	$max=$subi;
              	    if($subi eq int(($subi+$min)/2) and $iteration>20)
              	      {
               	        if(!$bla{$m}->[7])
              	          {
                  	  	    assign($m,"out",$subi,$subj);
                  	  	  }
                	      $match=1;	
              	      }                	  	       	    
              	  }     
              	print "$site\t$bla{$m}->[7]\n" if($feature{start}->[$subj][$subi][2] and $feature{start}->[$subj][$subi][2]=~/"19174\."/);          	
              } 
            last;
          }
      }  
  }

sub others
  {
  	my $m=shift;#blast record id
    my $chr=shift;#chromosome id
    my $id=shift;	#unit id
    if($feature{type}->[$chr][$id]=~/other|pseudo/)
      {
        assign($m,"other",$chr,$id);
        return(1,1);
      }
    else
      {
        return(0,0);
      }
  }

sub exon
  {
    my $chr=shift;#chromosome id
    my $id=shift;	#unit id
    my $inexo=0;  
    my $m=shift;#blast record id  
    my $abnormal=shift;#set value for abnomal ones  
    my $site=$bla{$m}->[5];#

    my @exonm=split /\,/,$feature{start}->[$chr][$id][2];#get each exon part
    EXON: foreach my $exo(0..@exonm-1)#deal with each exon part	
      {
        my @exon=split /\.\./,$exonm[$exo];#get the start/end of the exon
        print "exon:$exon[0]\t$exon[1]" if($exonm[$exo]=~/"19174"/);
        if($exon[1] and $exon[0])
          {
            if($site <= $exon[1] and $site >= $exon[0])#exon
              {
                $inexo=1;
                last EXON;
              }
          }
        elsif(!$exon[1] or !$exon[0])
          {
            assign($m,"other",$chr,$id) if(!$abnormal);
          }
      }
    if($inexo)
      {
        assign($m,"exon",$chr,$id);
      }
    else
      {
        assign($m,"intron",$chr,$id) if (!$abnormal);        	
      }	
  }

sub assign
 {
   my $m=shift;#blast id
   $bla{$m}->[7]=shift;#relative location
   $bla{$m}->[8]=shift;#function unit id
   $bla{$m}->[9]=shift;#chromosome id
 }
 
open (my $out, ">annotation.re ") or die $!;
open (my $outt, ">annotation_all_count.re") or die $!;
print $out "read_id\tcount\tread_id_full\tgenome\tidentity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\trelative location\tstrain\torganelle\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tfunction_unit_function\tstrand\tfunction_unit_note\n";
print $outt "read_id\tcount\tread_id_full\tgenome\tidentity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\trelative location\tstrain\torganelle\tfunction_unit_type\tfunction_unit_id\tchromosome\tfunction_unit site_easy\tfunction_unit site_complex\tfunction_unit_function\tstrand\tfunction_unit_note\n";
foreach my $key(sort {$a<=>$b} keys %bla)
  {
  	#annotate the blast file
  	annotate($key,0);
  	#abnormal one
  	annotate($key,1); 
  	
    if($bla{$key}->[1]>=$cou)
  	  {
  	  	prin_out(1,$out,$key);
  	  }
    prin_out(1,$outt,$key);      
  }
close $out;
close $outt;
  
sub prin_out
  {
  	my $for_bla=shift;
  	my $filehand=shift;
  	my $key;#blast id
  	$key=shift if($for_bla);
  	
  	my $unit;#unit id
  	$unit=shift if(!$for_bla);  	
  	$unit=$bla{$key}->[8] if($for_bla);
  	
  	my $chr;#chromosome id
  	$chr=shift if(!$for_bla);  	
  	$chr=$bla{$key}->[9] if($for_bla);  
	
    print $filehand "$key\t$bla{$key}->[1]\t$bla{$key}->[0]\t$bla{$key}->[7]\t" if($for_bla); #bla information 
    if($feature{organism}->[$chr])
      {
      	print $filehand "$feature{organism}->[$chr]\t";
      }
    else
      {
        print $filehand "\t";
      }      
    if($feature{organelle}->[$chr])
      {
      	print $filehand "$feature{organelle}->[$chr]\t";
      }
    else
      {
        print $filehand "\t";
      }   
    if($feature{type}->[$chr][$unit])#product
      {
      	print $filehand "$feature{type}->[$chr][$unit]\t";
      }
    else
      {
        print $filehand "mit_protein\t";
      }                   
    if($feature{name}->[$chr][$unit])#product_id
      {
      	print $filehand "$feature{name}->[$chr][$unit]\t";
      }
    else
      {
        print $filehand "\t";
      }           
    #if($feature{strain}->[$chr])
    #  {
    #  	print $filehand "$feature{strain}->[$chr]\t";
    #  }
    #else
    #  {
    #    print $filehand "\t";
    #  }
    #if($feature{mol_type}->[$chr])
    #  {
    #  	print $filehand "$feature{mol_type}->[$chr]\t";
    #  }
    #else
    #  {
    #    print $filehand "\t";
    #  }
    if($feature{chromosome}->[$chr])#chromosome
      {
      	print $filehand "$feature{chromosome}->[$chr]\t";
      }
    else
      {
        print $filehand "\t";
      }
    if(exists $feature{start}->[$chr][$unit][0])#site_easy/tsite complex
      {
      	print $filehand "$feature{start}->[$chr][$unit][1]\t$feature{start}->[$chr][$unit][2]\t" if($feature{start}->[$chr][$unit][0]);
        print $filehand "$feature{start}->[$chr][$unit][1]\t\t" if(!$feature{start}->[$chr][$unit][0]);
        die "$chr\t$unit\n" if(!$feature{start}->[$chr][$unit][1]);
      } 
    else
      {
        print $filehand "\t\t";	
      }
    if($feature{product}->[$chr][$unit])
      {
      	print $filehand "$feature{product}->[$chr][$unit]\t";
      }
    else
      {
        print $filehand "\t";
      }   
    if($feature{reverse}->[$chr][$unit]) 
      {
        print $filehand "-\t";	
      }
    else
      {
        print $filehand "+\t";	
      }
    #if($feature{protein_id}->[$chr][$unit])
    #  {
    #  	print $filehand "$feature{protein_id}->[$chr][$unit]\t";
    #  }
    #else
    #  {
    #    print $filehand "\t";
    #  }      
    if($for_bla)
      {
        print $filehand "$bla{$key}->[10]\t" if($bla{$key}->[10]);
        print $filehand "\t" if (!$bla{$key}->[10]);	
      }
    if($feature{note}->[$chr][$unit])
      {
      	print $filehand "$feature{note}->[$chr][$unit]\t";
      }
    else
      {
        print $filehand "\t";
      }
    print $filehand "\n";
  }

