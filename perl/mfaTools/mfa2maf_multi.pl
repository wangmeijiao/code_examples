use strict;
use warnings;
use Data::Dumper;


my ($info,undef) = &readInfo($ARGV[0]);
my @order = ("japo","glab","punc","brac","lper","sorg");
my $len = &getLenAll($info,\@order);

#print Dumper $len,$info; exit;


$/="\n\n";
my $cnt=0;

print "##maf version=1 scoring=mfa2maf\n";
while(<stdin>){
    chomp;
    next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
    my @box=split/>/,$_;
    die "empty line at $_" if(scalar @box == 0);
    my %aln;
    print "a score=$cnt\n";
    foreach my $str(@box){
       next if($str eq "" || $str eq "\n");
       my @temp = split/\n+/,$str;
       die "empty line at |$str|" if(scalar @temp == 0);
       my ($cord,$anno) = split/[\t ]+/,shift( @temp);
       my ($spec,$chr,$gene_id) = split/\|/,$cord;
       my ($tss,$strand,$extract) = split/\|/,$anno;
       my (undef,$region_s,$region_e) = split/:|-/,$extract;
       my $seq = join"",@temp;
       my $len_aln = length $seq;
       my $seq_real = $seq;
       $seq_real =~s/-+//g;
       my $len_real = length $seq_real;
       if(!exists $aln{$spec}){ 
          $aln{$spec}->{"chr"} = $chr; 
          $aln{$spec}->{"gene_id"} = $gene_id; 
          $aln{$spec}->{"strand"} = $strand; 
          $aln{$spec}->{"region_s"} = $region_s; 
          $aln{$spec}->{"region_e"} = $region_e; 
          $aln{$spec}->{"anno"} = $anno; 
          $aln{$spec}->{"aln"}=$seq;
          $aln{$spec}->{"seq"}=$seq_real;
          $aln{$spec}->{"alnLen"}=$len_aln;
          $aln{$spec}->{"realLen"}=$len_real;
       }else{die "dup spec $spec at $_"}
    }
    #print Dumper \%aln;
    foreach my $spec(@order){
       if(!exists $aln{$spec}){}else{
         my $chr = $aln{$spec}->{"chr"};
         #printf("s %-15s %-8s %-5s %1s %-8s %s\n", $spec.".".$chr, $aln{$spec}->{'region_s'}, $aln{$spec}->{'realLen'},$aln{$spec}->{'strand'},$len->{$spec}->{$chr},$aln{$spec}->{'aln'});
         printf("s %s %s %s %s %s %s\n", $spec.".".$chr, $aln{$spec}->{'region_s'}, $aln{$spec}->{'realLen'},$aln{$spec}->{'strand'},$len->{$spec}->{$chr},$aln{$spec}->{'aln'});
       }
    }
    print "\n";

    $cnt++;
    #if($cnt == 1){last}

}

$/="\n";



##sub

sub readInfo(){ #genome info file with spec order
  my %info;
  my @order;
  my $file = shift @_;
  open INF, $file or die"$!";
  $/="#";
  <INF>;
  while(<INF>){
   chomp;
   my @box=split/\n+/,$_;
   my $species=shift @box;
   push @order,$species;
   foreach(@box){
     my ($type,$file)=split/[\t ]+/,$_;
     next if($type =~/^"/);
     if(!exists $info{$species}->{$type}){ $info{$species}->{$type}=$file }else{die "dup $type in species $species\n"}
   }
  }
  close INF;
  $/="\n";
  return (\%info,\@order);
}

sub readLen(){
    my $file=shift @_;
    my %len;
    open LEN, $file or die "$!";
    while(<LEN>){
      chomp;
      next if ($_ eq "" || $_=~/^#/);
      my ($chr,$len)=split/\t/,$_;
      if(!exists $len{$chr}){$len{$chr}=$len}else{die "dup chr $chr\n"}
    }
    close LEN;
    return \%len;
}


sub getLenAll(){
  my ($info,$archi)=@_;
  my %len;
  foreach my $species(@{$archi}){
   my $sizef=$info->{$species}->{"faSize"};
   if ( -e $sizef ){}else{print "\n$species\t$sizef\n"; die "not exists $sizef\n"}
   my $len=&readLen($sizef);  
   $len{$species} =  $len;
  }
  return \%len;
}



