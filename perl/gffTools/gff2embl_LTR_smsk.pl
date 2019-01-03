use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#embl format:
#
#
#FH       key         Location/Qualifiers
#
#FT       feature_key   regions
#FT                     /qualifier1="  "
#FT                     /qualifier2="  "
#   ...
#
#SQ     sequence 1000 Bp ; ...
#       CGTAGCTACGAT   GTAGTTTCGATG  ACGTAGTCACAG  36
#          ...

#//

# valid feature_keys/qualifiers groups see  ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html


# usage: perl thisfile.pl  test.LTR.gff > test.LTR.gff.embl


my ($gff_idx,$repeat_idx)=&gffRead_repeat($ARGV[0]);

#print Dumper $gff_idx,$repeat_idx;
#exit;


printf("%-5s%-18s%s\n", "FH","Key","Location/Qualifier");
printf("%-23s\n", "FH");


#=pod
foreach my $chr(keys %{$gff_idx}){
  foreach my $repeatID(sort keys %{$gff_idx->{$chr}}){
         my ($chr,$start,$end,$ID,$strand,$note)=split/\t/,$repeat_idx->{$repeatID};
         my $element=$gff_idx->{$chr}->{$repeatID};
         if(exists $element->{"repeat_fragment"}){
            my @regions=@{$element->{"repeat_fragment"}};
            my $coord=&getCoord(\@regions,$strand);
            printf("%-5s%-18s%s\n","FT","LTR",$coord);         
            printf("%-23s%s\n","FT","/ID=\"$ID\"");         
            #printf("%-23s%s\n","FT","/note=\"$note\"");         
            foreach(@{&segSort_start(\@regions)}){
              my @box=split/-/,$_;
              shift @box;
              shift @box;
              my $str=join("-",@box);
              printf("%-23s%s\n","FT","/note=\"$str\"");
            }
         }else{
                my @regions=("$start-$end");
                my $coord=&getCoord(\@regions,$strand);
                printf("%-5s%-18s%s\n","FT","LTR",$coord);          
                printf("%-23s%s\n","FT","/ID=\"$ID\"");
                printf("%-23s%s\n","FT","/note=\"$note\"");
              }
  }#foreach Transposon
}#foreach chr

#=cut





###sub##


sub gffRead_repeat(){
  ##modified standard nest gff3 parsing code for repeat gff 
  my $file=shift;
  my %gff;#main structure table:  chr->Transposon->repeat_fragment
  my %repeat;#repeat info
  #1,read&fill without order consideration
  open GFF, $file or die "$!";
  while(<GFF>){
    chomp;
    next if($_ eq "" || $_=~/^#/);
    my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
    die"empty value at $_\n" if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
    next if($chr eq "chrUn" || $chr eq "chrSy" || $chr eq "Pt" || $chr eq "Mt"); #get rid of other chrs
    next if($ftype eq "chromosome" || $ftype eq "Chromosome");
    ($start,$end)=($start<$end)?($start,$end):($end,$start);
    my ($ID,$pID,$note);
    if($ftype eq "Transposon"){
      $attr=~/ID=([^;]+);Note=([^;]+)/;
      ($ID,$note)=($1,$2);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $repeat{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$note";
      }else{die "dup Transposon name $ID\n"}
     }elsif($ftype eq "repeat_fragment" ){
           $attr=~/Parent=([^;]+);Note=([^;]+)/;
           ($pID,$note)=($1,$2);
           if(exists $repeat{$pID}){}else{die "elements comes first, I can't assign it to Transposon:$pID\n"}
           if(!exists $gff{$chr}->{$pID}->{$ftype}){
              my @temp;
              push @temp,"$start-$end-$note";
              $gff{$chr}->{$pID}->{$ftype}=\@temp;
           }else{ push @{$gff{$chr}->{$pID}->{$ftype}},"$start-$end-$note" }
         }else{die "unknow feature type\n"}
  }#while end
  close GFF;
  return (\%gff,\%repeat);

}


sub getCoord(){
   my ($region,$strand)=@_;
   my $string;
   my @regions_sorted=@{&segSort_start($region)};
   if($strand eq "+"){
     if(scalar @regions_sorted == 1){
       my ($s,$e)=split/-/,$regions_sorted[0];
       $string="$s..$e";
     }elsif(scalar @regions_sorted > 1){
       my @str;
       foreach my $seg(@regions_sorted){
         my ($s,$e)=split/-/,$seg;
         push @str,"$s..$e";
       }
       $string=join ",",@str;
       $string="join($string)";
     }else{die "epmty segs at @regions_sorted"}
   }elsif($strand eq "-"){
       if(scalar @regions_sorted == 1){
         my ($s,$e)=split/-/,$regions_sorted[0];
         $string="complement($s..$e)";
       }elsif(scalar @regions_sorted > 1){
          my @str;
          foreach my $seg(@regions_sorted){
            my ($s,$e)=split/-/,$seg;
            push @str,"$s..$e";
          }
          $string=join ",",@str;
          $string="complement(join($string))";
        }else{die "epmty segs at @regions_sorted"}
     }else{die "unknow strand $strand at @{$region}"}
   return $string;
}


sub segSort_start(){ #simple sort by start
     #start < end
     my @segs=@{shift @_};
     my @segs_sorted=sort{
       my ($s1,$e1)=split/-/,$a;
       my ($s2,$e2)=split/-/,$b;
       die "s > e at @segs" if($s1 > $e1 || $s2 > $e2);
       if($s1<$s2){return -1}elsif($s1 > $s2){return 1}else{return 0}
     } @segs;
     return \@segs_sorted;
}




