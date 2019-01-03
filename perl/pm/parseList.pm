
use strict;
use warnings;
use Data::Dumper;


my ($mark_o,$mark_p,$mark_f) = &parseList_new($ARGV[0]);
print Dumper $mark_o,$mark_p,$mark_f;



sub parseList_new(){
#list format:
#9       input_k271      /home/mjwang/pwdex/epi-oryza-data/6libs_BGI_1122/1115-input-rice/input-k271_tigr6_bowtiev0M1/hits_input-k271_tigr6_bowtiev0M1.density.gz
#9       k271    /home/mjwang/pwdex/epi-oryza-data/6libs_BGI_1122/1115-k271-rice/k271_tigr6_bowtiev0M1/hits_k271_tigr6_bowtiev0M1.density.gz 
  my $list=shift;
  #read list
  my %mark_list;
  my %mark_file;
  open LIST, $list or die "$!";
  while(<LIST>){
   chomp;
   next if($_ eq "" || $_=~/^#/);
   my($num,$mark,$file)=split/[\t ]+/,$_;
   next if ($file eq "NA");
   if($mark=~/input/){
     if(!exists $mark_list{$num}->{'control'}){$mark_list{$num}->{'control'}=$mark}else{die "dup mark control $mark"}
   }else{
      if(!exists $mark_list{$num}->{'signal'}){$mark_list{$num}->{'signal'}=$mark}else{die "dup mark signal $mark"}
    }
   if(!exists $mark_file{$mark}){ $mark_file{$mark} = $file }else{die "dup mark file $mark:$file"}
  }
  close LIST;
  #pair & order them
  my @mark_order;
  my %pair;
  foreach my $i(sort {$a<=>$b} keys %mark_list){  
     if(exists $mark_list{$i}->{'signal'} && exists $mark_list{$i}->{'control'}){
        push @mark_order,$mark_list{$i}->{'signal'};
        $pair{$mark_list{$i}->{'signal'}} = $mark_list{$i}->{'control'}; 
     }elsif(exists $mark_list{$i}->{'signal'} && !exists $mark_list{$i}->{'control'}){
          push @mark_order,$mark_list{$i}->{'signal'};
          $pair{$mark_list{$i}->{'signal'}} = undef;
       }else{die "not valid signal-control pair at group $i"}
  }

  return \@mark_order,\%pair,\%mark_file;
}



