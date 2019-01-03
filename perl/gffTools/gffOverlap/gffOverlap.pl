use strict;
use warnings;
use Data::Dumper;

#usage : cat oryza_nivara.maker.gff |perl gffOverlap.pl OnivaChr04:10000-500000 > OnivaChr04-10000-500000.gff

my ($chr,$s,$e)=split/:|-/,$ARGV[0];
my @gff;
#read, extract and store
while (<stdin>){
  chomp;
  next if($_ eq "" || $_=~/^#/);
  my @box=split/\t/,$_;
  die "start >= end at @box\n" if($box[3] > $box[4]);
  next if($box[3] == $box[4]);
  next if($box[2] eq "chromosome");
  next if($box[0] ne $chr || $box[4] <= $s || $box[3] >= $e );
  push @gff, $_;
}

#check gff validation and output

if(1){ foreach (@gff){print "$_\n"} }else{ print STDERR "\nerr in:\n",join("\n",@gff)  }
#if(! &gffCheck(\@gff)){ foreach (@gff){print "$_\n"} }else{ print STDERR "\nerr in:\n",join("\n",@gff)  }


##sub##
sub gffCheck(){
  my $gff=shift;
  my $flag=0;  
  my %gff;
  my %gene;
  my %mrna;
  foreach(@{$gff}){
    my ($chr,undef,$ftype,$start,$end,undef,$strand,$fshift,$attr)=split/\t/,$_;
    if(!defined $chr || !defined $ftype || !defined $start || !defined $end || !defined $strand || !defined $attr){$flag=1;print STDERR "undefined columns found \n"}
    if($chr eq "chrUn" || $chr eq "chrSy"){$flag=1; print STDERR "chromosome err\n"} #get rid of other chrs
    if($ftype eq "chromosome" || $ftype eq "Chromosome"){$flag=1; print STDERR "chromosome err\n"}
    if($start >= $end){$flag=1; print STDERR "start >= end\n"}
    my ($ID,$pID,$note,$name);
    if($ftype eq "gene"){
      $attr=~/ID=([^;]+);Name=([^;]+);biotype=([^;]+)/;
      ($ID,$name,$note)=($1,$2,$3);
      if(!exists $gff{$chr}->{$ID}){
        my %temp;
        $gff{$chr}->{$ID}=\%temp;
        $gene{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$note";
      }else{$flag=1;print STDERR "repeat gene name $ID\n"}
     }elsif($ftype eq "mRNA"){
        $attr=~/ID=([^;]+);Parent=([^;]+)/;
        ($ID,$pID)=($1,$2);
        if(!exists $gff{$chr}->{$pID}->{$ID}){
          my %temp;
          $gff{$chr}->{$pID}->{$ID}=\%temp;
          $mrna{$ID}="$chr\t$start\t$end\t$ID\t$strand\t$pID";
        }else{$flag=1; print STDERR "repeat mRNA ID at $ID\n"}
       }elsif($ftype eq "CDS" || $ftype eq "exon" ||$ftype eq "intron"){
           $attr=~/Parent=([^;]+);/;
           $pID=$1;
           my $geneID;
           if(exists $mrna{$pID}){
            my @box=split/\t/,$mrna{$pID};
            if(defined $box[5]){$geneID=$box[5]}else{$flag=1; print STDERR "pID empty at $pID\n"}
           }else{$flag=1; print STDERR "elements comes first, I can't assign it to gene\n"}
  if(!exists $gff{$chr}->{$geneID}->{$pID}->{$ftype}){
              my @temp;
              if($ftype eq "CDS"){push @temp,"$start-$end-$fshift"}else{push @temp,"$start-$end"};
              $gff{$chr}->{$geneID}->{$pID}->{$ftype}=\@temp;
           }else{if($ftype eq "CDS"){push @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end-$fshift"}else{push  @{$gff{$chr}->{$geneID}->{$pID}->{$ftype}},"$start-$end"}}
         }else{$flag = 1 ;print STDERR "unknow feature type\n"}
  }#foreach end 
  return $flag;
}


