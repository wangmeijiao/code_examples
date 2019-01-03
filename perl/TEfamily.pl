use strict;
use warnings;
use Data::Dumper;
use File::Path;

#classify te family fasta, link LTR to full length LTR retrotransposons
#repeatmasker lib header format id#family/subfamily
#export file structure:  
#
#                 family1_dir-----subfamily1.fa (multiple seqs, LTR retrotransposon will be linked as LTR-int-LTR)
#                            -----subfamily2.fa
#                            ...
#                 family2_dir-----subfamily1.fa
#                            -----subfamily2.fa
#                            -----subfamily3.fa

my %tefamily;
my %LTR;


$/=">";
while(<stdin>){
   chomp;
   next if($_ eq "" || $_ =~/^#/ || $_ =~/^\s+$/);    
   my @box=split/\n+/,$_;
   my $head=shift @box;
   my @ids=split/ +/,$head;   
   my ($id,$familyID,$subfamID)=split/#|\//,$ids[0];
   if(!defined $id || $id eq ""){die "no id found at $_"}   #id must not empty
   if(!defined $familyID || $familyID eq ""){die "no Familyname found at $_"}  #familyname must not empty
   if(!defined $subfamID || $subfamID eq ""){$subfamID="noSubFamname"}
   my $seq=join"\n",@box;
   $seq=uc($seq);  #upcased

   if($familyID eq "LTR"){
      my $name=$id;
      $name=~s/_OS$//i;
      $name=~s/[-_]LTR$|[-_]I$|[-_]I[-_]int$|[-_]int$//i;
      #$name=~s/-LTR$|-I$|_LTR$|_I$|-int$|_int$|_I-int$//i;
      $LTR{$familyID}->{$subfamID}->{$name}->{$id}=$seq;
   }else{  
         if(exists $tefamily{$familyID}->{$subfamID}->{$id}){die "dup id $id at $_"}else{$tefamily{$familyID}->{$subfamID}->{$id}=$seq}

    }


}#while end
$/="\n";


#print Dumper \%LTR,\%tefamily; 
#exit;




####report LTR, both ori and linked##########################################################
foreach my $familyID(sort keys %LTR){
  #if(-d $familyID ){}else{if(mkdir($familyID) == 1){}else{die "can't mkdir $familyID"}}
  print "$familyID\n";
  foreach my $subfamID(sort keys %{$LTR{$familyID}}){
    print "=======$subfamID\n";
    mkpath("$familyID/$subfamID"); 
    foreach my $id(sort keys %{$LTR{$familyID}->{$subfamID}}){
       print "           $id\n";
       open OUT, ">$familyID/$subfamID/$id.fa" or die "$!"; #ori multi seqs, int and LTR       
       my %full_ltr;
       my $full_ltr_seq;
       foreach my $name(keys %{$LTR{$familyID}->{$subfamID}->{$id}}){ 
         print OUT ">$name\n$LTR{$familyID}->{$subfamID}->{$id}->{$name}\n";
         my $string=$name;
         $string=~s/_OS//i;
         if($string=~/[-_]LTR$/i ){
           if(!exists $full_ltr{"ltr"}){$full_ltr{"ltr"}=$LTR{$familyID}->{$subfamID}->{$id}->{$name}}else{die "dup ltr seq at $name"}
           open LTR, ">$familyID/$subfamID/${id}_LTR.fa" or die "$!"; # LTR seq only       
           print LTR ">$id\n$LTR{$familyID}->{$subfamID}->{$id}->{$name}\n";
         }elsif($string=~/[-_]I$|[-_]I[-_]int$|[-_]int$/i){
           if(!exists $full_ltr{"int"}){$full_ltr{"int"}=$LTR{$familyID}->{$subfamID}->{$id}->{$name}}else{die "dup int seq at $name"}

           }else{ if(!exists $full_ltr{"unknown"}){$full_ltr{"unknown"}=$LTR{$familyID}->{$subfamID}->{$id}->{$name}}else{die "dup not ltr/int seq at $name"}}
       }
 
       #link to produce full ltr retrotransposon
       if(exists $full_ltr{"ltr"} && exists $full_ltr{"int"} && !exists $full_ltr{"unknown"}){
          $full_ltr_seq=$full_ltr{"ltr"}."\n".$full_ltr{"int"}."\n".$full_ltr{"ltr"};
          $full_ltr_seq=~s/\n+//g;
          $full_ltr_seq=&faFormat($full_ltr_seq,60);
          open FULL, ">$familyID/$subfamID/${id}_full.fa" or die "$!"; #linked full LTR 
          print FULL ">$id\n$full_ltr_seq";
       }else{print STDERR "$id#$familyID/$subfamID is not full_LTR_retrotransposon: ",join(",",keys %full_ltr),"\n"}
       #}elsif(!exists $full_ltr{"ltr"} && !exists $full_ltr{"int"} && exists $full_ltr{"unknown"}){
           #$full_ltr_seq=$full_ltr{"unknown"};
           #print FULL ">$id\n$full_ltr_seq\n";
        # }elsif(exists $full_ltr{"ltr"} && !exists $full_ltr{"int"} && !exists $full_ltr{"unknown"}){
            # $full_ltr_seq=$full_ltr{"ltr"};
            # print FULL ">${id}_LTR\n$full_ltr_seq\n";
         #  }elsif(!exists $full_ltr{"ltr"} && exists $full_ltr{"int"} && !exists $full_ltr{"unknown"}){
               # $full_ltr_seq=$full_ltr{"int"};    
                #print FULL ">${id}-I\n$full_ltr_seq\n";
          #    }else{ #print STDERR "unknown circurence at $id:\n";print STDERR join" ",keys( %full_ltr),"\n"; system("rm $familyID/$subfamID/${id}_full.fa");print STDERR "$familyID/$subfamID/${id}_full.fa removed\n"
             #      }
       close OUT;
       close FULL;
       close LTR;
    }#foreach id    
    #last;
  }#foreach subfamiID 
  #last;
}#foreach familyID


######report nonLTR###############################################
foreach my $familyID(sort keys %tefamily){
  #if(-d $familyID ){}else{if(mkdir($familyID) == 1){}else{die "can't mkdir $familyID"}}
  print "$familyID\n";
  foreach my $subfamID(sort keys %{$tefamily{$familyID}}){
    print "=======$subfamID\n";
    mkpath("$familyID/$subfamID");
    foreach my $id(sort keys %{$tefamily{$familyID}->{$subfamID}}){
       print "           $id\n";
       open OUT, ">$familyID/$subfamID/$id.fa" or die "$!";
       #foreach my $name(keys %{$tefamily{$familyID}->{$subfamID}->{$id}}){
       print OUT ">$id\n$tefamily{$familyID}->{$subfamID}->{$id}\n";

       #}
       close OUT;
    }
    #last;
  }#foreach subfamiID 
  #last;

}




##sub#
sub faFormat(){
  my ($seq,$num)=@_;
  my @seq=split//,$seq;
  my ($string,$cnt);
  foreach(@seq){
   $cnt++;
   if($cnt%$num==0){$string.=$_."\n"}else{$string.=$_}
  }
  if($cnt%$num!=0){return $string."\n"}else{return $string}
}





