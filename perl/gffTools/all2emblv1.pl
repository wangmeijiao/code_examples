use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#desc: transform geneGff (default), teGff and fa to a single embl file; one seq only
#embl contents: gene records
#               te  records
#               fa seqs
#usage: perl thisfile.pl -geneGff gene.gff -teGff te.gff -fa in.fa > out.embl


my ($geneGff,$teGff,$fa,$gapGff);
GetOptions("geneGff=s",\$geneGff,"teGff:s",\$teGff,"fa:s",\$fa,"gapGff:s",\$gapGff);


#read and combine

if(defined $geneGff){
   #print Dumper &readGff_gene($geneGff); exit;
   &geneGff2embl(&readGff_gene($geneGff));
}else{die "gene gff file can't be empty\n"}


if(defined $teGff){ 
   my @files = split/,/,$teGff; 
   foreach my $file (@files){
     &teGff2embl(&readGff_te($file)) 
   }

 }

if(defined $gapGff){ &gapGff2embl(&readGff_gap($gapGff))   }

if(defined $fa){ &formatFaEMBL(&readFa($fa),60)  }




###sub###

sub readGff_gene(){
   #only capture mrna and CDS and or exon lines
   my $file = shift;
   die "file $file not exists\n" if(! -f $file);
   my %gff;
   open GFF, $file or die "$!";
   while(<GFF>){  
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my ($chr_in,$foo1,$ftype,$start,$end,$foo2,$strand,$fshift,$attr)=split/\t/,$_;
     die"empty value at $_\n" if(!defined $chr_in || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
     ($start,$end)=($start<$end)?($start,$end):($end,$start);
     my ($id,$pid);
     if($ftype eq "mRNA"){
        $attr=~/ID=([^;]+)/;
        $id=$1;
        if(!exists $gff{$id}){
          $gff{$id}->{"range"}="$start-$end";
          $gff{$id}->{"strand"}=$strand;
        }else{die "CDS comes before mrna, quit\n"}
     }elsif($ftype eq "CDS"){
        $attr=~/Parent=([^;]+)/;
        $pid=$1;
        if(exists $gff{$pid}){
          push @{$gff{$pid}->{"CDS"}},"$start-$end";
        }else{die "CDS comes before mrna, quit\n"}
      }elsif($ftype eq "exon"){
         $attr=~/Parent=([^;]+)/;
         $pid=$1;
         if(exists $gff{$pid}){
            push @{$gff{$pid}->{"exon"}},"$start-$end";
         }else{die "exon comes before mrna, quit\n"}
       }
   }#while end
   close GFF;
   return \%gff;
}


sub readGff_te(){
   #single line 10 columns records
   my $file = shift;
   die "file $file not exists\n" if(! -f $file);
   my %gff;
   open GFF, $file or die "$!";
   while(<GFF>){ 
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my ($chr_in,$foo1,$ftype,$start,$end,$foo2,$strand,$fshift,$attr)=split/\t/,$_;
     die"empty value at $_\n" if(!defined $chr_in || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
     ($start,$end)=($start<$end)?($start,$end):($end,$start);
     my $id;
     $attr=~s/;Note=/\|/;
     $attr=~/ID=([^;]+)/;
     $id=$1;
     if(!exists $gff{$id}){ $gff{$id}="$id\t$ftype\t$start\t$end\t$strand"  }else{die "dup $id in te gff $file\n"};
   }
   close GFF;
   return \%gff;
}


sub readGff_gap(){
   #single line 10 columns records
   my $file = shift;
   die "file $file not exists\n" if(! -f $file);
   my %gff;
   open GFF, $file or die "$!";
   while(<GFF>){ 
     chomp;
     next if($_ eq "" || $_=~/^#/);
     my ($chr_in,$foo1,$ftype,$start,$end,$foo2,$strand,$fshift,$attr)=split/\t/,$_;
     die"empty value at $_\n" if(!defined $chr_in || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr);
     ($start,$end)=($start<$end)?($start,$end):($end,$start);
     my $id;
     $attr=~/ID=([^;]+)/;
     $id=$1;
     if(!exists $gff{$id}){ $gff{$id}="$id\t$ftype\t$start\t$end\t$strand"  }else{die "dup $id in te gff $file\n"};
   }
   close GFF;
   return \%gff;
}






sub readFa(){
  #single fa only
  my $file=shift;
  my %fa;
  $/=">";
  open FA, $file or die "$!";
  while(<FA>){
    chomp;
    next if($_ eq ""|| $_=~/^#/);
    my @box=split/\n+/,$_;
    my $id=shift @box;
    my ($chr,$s,$e)=split/:|-/,$id;
    my $seq=join"",@box;
    $fa{$chr}=$seq;
  }
  close FA;
  $/="\n";
  return \%fa;
}


sub formatFaEMBL(){
  #embl sequence embeded
  #one seq only
  #SQ   Sequence 324123 BP;
  #     ttcacacaac aatatacgat acaacaatgg agttccattg cgagaagaaa cagtaatctt    871200     # 10*6 
  #     tcagaaaaca aaacgactca ataatggaaa tccatcacaa gtagaaacat caatctttca    871260
  #     cacc                                                                 871264
  #//
  my ($refseq,$width)=@_;
  my @head=keys %$refseq;
  my $seq=$refseq->{$head[0]};
  die "uncomfortable line width $width, please use number > 10 and number % 10 ==0\n" if($width <= 10 || $width % 10 != 0);
  my @seqs=split//,$seq;
  my $len=length $seq;
  die "empty sequence $seq\n" if($len == 0);
  my $cnt;
  print  "SQ   Sequence $len BP;\n     ";
  foreach(@seqs){
     next if($_ eq " " || !defined $_ || $_ eq "");
     $cnt++;
     if($cnt % 10 ==0){print "$_ "}else{print "$_"}
     if($len % $width == 0 && $cnt == $len){
        if($cnt % $width ==0){printf("%8i\n",$cnt)} 
     }else{if($cnt % $width ==0){printf("%8i\n     ",$cnt)}}

  }
  die "detected illegal char in your sequence $seq\n" if($len != $cnt);
  #handle end space with care
  if($len % $width != 0){
     my $tail_len=$len % $width;
     my $space_num= $width+($width/10-1) - $tail_len;
     for (1..$space_num){print " "}
     printf("%8i",$len);
     print "\n\/\/\n";
  }else{print "\/\/\n"}
  return 1;
}




sub fa2embl(){
  #embl sequence embeded
  #one seq only
  #SQ   Sequence 324123 BP;
  #     ttcacacaac aatatacgat acaacaatgg agttccattg cgagaagaaa cagtaatctt    871200     # 10*6 
  #     tcagaaaaca aaacgactca ataatggaaa tccatcacaa gtagaaacat caatctttca    871260
  #     cacc                                                                 871264
  #//
  my $refseq=shift; 
  my @head=keys %$refseq;
  my $seq=$refseq->{$head[0]};
  my $len = length ($seq);
  my $mod = $len % 60;
  my $times = ($len - $mod) / 60;
  print  "SQ   Sequence $len BP;\n";
  my ($i,$j,$seg);
  for ($i = 0; $i < $times; $i++) {
          $seg = substr ($seq, $i * 60, 60);
          print  "    ";
          for ($j = 0; $j < 6; $j++) {
                  print " ", substr ($seg, $j * 10, 10);
          }
          printf  ("%10d\n", $i * 60 + 60);
  }
  #handle the tail carefully
  my $tail = substr ($seq, $len - $mod);
  print  "    ";
  my $l="";
  for ( $j = 0; $j < 6; $j++) {
          print " ", substr ($tail, $j * 10, 10);
          $l = $l." ".substr ($tail, $j * 10, 10);
  }
  $l = length ($l) + 4;
  $l = 70 - $l;
  for ( $j = 0; $j < $l; $j++) {
          print  " ";
  }
  printf  ("%10d\n", $i * 60 + $mod);
  print  "\/\/\n";
  return 1;
}

sub geneGff2embl(){
   #only include mRNA regions(exons) and CDS regions
   #
   #FT   mRNA            join(480444..482052,489311..489830,489934..491323) #exons
   #FT                   /gene="B6-4"
   #FT   CDS             join(481224..482052,489311..489830,489934..490522)
   #FT                   /gene="B6-4"
   #FT                   /translation=""
   my $table=shift;
   foreach my $id(sort keys %{$table}){
     die "trunct data in gff of $id: range or strand not found: $table->{$id}->{'range'},$table->{$id}->{'strand'}" if(!exists $table->{$id}->{"range"} || !exists $table->{$id}->{"strand"});
     if(!exists $table->{$id}->{"CDS"}){
       print STDERR "trunct data in gff of $id: CDS not found, use mrna range to add \n";
       push @{$table->{$id}->{"CDS"}},$table->{$id}->{"range"} 
     }
     if(!exists $table->{$id}->{"exon"}){
       print STDERR "trunct data in gff of $id: exon not found, use mrna range to add \n";
       push @{$table->{$id}->{"exon"}},$table->{$id}->{"range"} 
     }

     my ($start_mrna,$end_mrna)=split/-/,$table->{$id}->{"range"};
     die "start >= end\n" if($start_mrna >= $end_mrna);
     my $strand=$table->{$id}->{"strand"};
     my @cds=@{&segSort($table->{$id}->{"CDS"})};
     my @exon;
     if(!exists $table->{$id}->{"exon"}){ # use mrna range + CDS to calculate exon if not exist exon
        if(scalar @cds == 1){ #single cds
           my ($s,$e)=split/-/,$cds[0];
           if($s > $start_mrna){ $s=$start_mrna }
           if($e < $end_mrna){ $e=$end_mrna }
           push @exon,"$s-$e";
        }elsif(scalar @cds > 1){ #more than two cds
             foreach(0..$#cds){
               my ($s,$e)=split/-/,$cds[$_];
               if($_ == 0){  
                    if($s > $start_mrna){ $s=$start_mrna };
                    push @exon,"$s-$e";
                 }elsif( $_ == $#cds ){
                      if($e < $end_mrna){ $e=$end_mrna };
                      push @exon, "$s-$e";
                    }else{push @exon, "$s-$e"}
             }
          }else{die "empty cds of $id: @cds\n"} #empty record
     }else{ @exon=@{&segSort($table->{$id}->{"exon"})} } #exon does exist

     #start to print embl records
     my @regions_cds;
     my @regions_exon;
     foreach(@cds){my ($s,$e)=split/-/,$_; push @regions_cds,"$s..$e"}     
     foreach(@exon){my ($s,$e)=split/-/,$_; push @regions_exon,"$s..$e"}     
     my $string_cds;
     my $string_exon;
     if($strand eq "+" || !defined $strand || $strand eq "."){ # plus strand or unknown strand
        #cds
        if(scalar @regions_cds == 1){
           $string_cds=$regions_cds[0];
        }elsif(scalar @regions_cds > 1){
            $string_cds=join(",",@regions_cds);
            $string_cds="join($string_cds)";
          }else{die "empty regions_cds \n"}
       #exon
       if(scalar @regions_exon == 1){
           $string_exon=$regions_exon[0];
        }elsif(scalar @regions_exon > 1){
            $string_exon=join(",",@regions_exon);
            $string_exon="join($string_exon)";
          }else{die "empty regions_exon \n"}
     }elsif($strand eq "-"){ #minus strand
             #cds
             if(scalar @regions_cds == 1){
                $string_cds="complement($regions_cds[0])";
             }elsif(scalar @regions_cds > 1){
                  $string_cds=join(",",@regions_cds);
                  $string_cds="complement(join($string_cds))";
                }else{die "empty regions_cds \n"}
            #exon
            if(scalar @regions_exon == 1){
                $string_exon="complement($regions_exon[0])";
              }elsif(scalar @regions_exon > 1){
                  $string_exon=join(",",@regions_exon);
                  $string_exon="complement(join($string_exon))";
                }else{die "empty regions_exon \n"}
         }else{die "unknown strand in $id\n"}

       print  "FT   mRNA            $string_exon\n";
       print  "FT                   \/gene\=\"$id\"\n";
       print  "FT   CDS             $string_cds\n";
       print  "FT                   \/gene\=\"$id\"\n";

  }#foreach id
  return 1;
}

sub segSort(){
    my $segs=shift;
    my (@segs,@segs_sorted);
    @segs=@{$segs};
    @segs_sorted=sort{ 
                   my ($s1,$e1)=split/-/,$a;
                   my ($s2,$e2)=split/-/,$b;
                   die "start > end in cds of $a,$b\n" if($s1 > $e1 || $s2 > $e2);
                   if($e2<=$s1){1}elsif($e1<=$s2){-1}else{print STDERR "unknown circumstance at $s1-$e1:$s2-$e2\n"}
                 }@segs;

    return \@segs_sorted;
}




sub teGff2embl(){
    #for LTR or DNA te
    #FT(3 space)LTR             complement(487351..487464)
    #FT                   /note="SZ-36_LTR|LTR"
    #FT                   /note="truncated"
    #FT(3space)repeat_region   complement(488191..488359)
    #FT                   /note="Explorer|DNA"
    #FT                   /note="truncated"
    my $table=shift;    
    foreach my $id(sort keys %{$table}){
        my ($id,$ftype,$start,$end,$strand)=split/\t/,$table->{$id};
        if($strand eq "+" || !defined $strand || $strand eq "."){
            printf ("FT   %3s             %i..%i\n",$ftype,$start,$end);
            print "FT                    /note=\"$id\"\n";
        }elsif($strand eq "-"){
             printf ("FT   %3s             complement(%i..%i)\n",$ftype,$start,$end);
             print "FT                    /note=\"$id\"\n";
          }else{die "unknown strand $strand\n"}
    }
    return 1;
}



sub gapGff2embl(){
    #for gap color = 11
    #FT   gap             complement(487351..487464)
    #FT                   /note="gapxxx"
    my $table=shift;    
    foreach my $id(sort keys %{$table}){
        my ($id,$ftype,$start,$end,$strand)=split/\t/,$table->{$id};
        if($strand eq "+" || !defined $strand || $strand eq "."){
            printf ("FT   %3s             %i..%i\n",$ftype,$start,$end);
            print "FT                    /note=\"$id\"\n";
            print "FT                    /colour=11\n";
        }elsif($strand eq "-"){
             printf ("FT   %3s             complement(%i..%i)\n",$ftype,$start,$end);
             print "FT                    /note=\"$id\"\n";
            print "FT                    /colour=11\n";
          }else{die "unknown strand $strand\n"}
    }
    return 1;
}




