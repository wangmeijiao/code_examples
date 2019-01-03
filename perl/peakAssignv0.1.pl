use strict;
use warnings;
use Data::Dumper;


my $bed_idx=&readBed_segTree($ARGV[0]);
my $prefix=$ARGV[1];
$prefix||="out";
#print Dumper $bed_idx;

#exit;

open LOST, ">$prefix.lost.list" or die "$!";
my $cnt;
while(<stdin>){
   chomp;
   next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
   $cnt++;
   print STDERR "#" if($cnt % 1000 == 0);  
   my ($chr,$start,$strand,$geneID)=split/[\t ]+/,$_;
   my $end=$start+1;
   my $idx=int 0.5*($start+$end)*1e-5;  #100kb range
   my $hits;
   ##in situ
   if(exists $bed_idx->{$chr}->{$idx}){  
     $hits=&segQuery_hash($bed_idx->{$chr}->{$idx},"$start-$end");
   }
   if(keys %{$hits} == 0){
       if($strand eq "+"){
       ##inward  step 200bp, 1kb most, stop when hits found
         for my $i(1..5){
            my $end_ext = $start+$i*200;
            my $idx_ext = int 0.5*($start+$end_ext)*1e-5;         
            if(!exists $bed_idx->{$chr}->{$idx_ext}){next}
            $hits=&segQuery_hash($bed_idx->{$chr}->{$idx_ext},"$start-$end_ext");
            if(keys %{$hits} != 0){last}
         } 
         ##outward step 500bp, 5kb most, stop when hits found 
         if(keys %{$hits} == 0){
            for my $i(1..10){
              my $start_ext = $start-$i*500; 
              if($start_ext <= 0){$start_ext=1} 
              my $idx_ext = int 0.5*($start_ext+$start)*1e-5;
              if(!exists $bed_idx->{$chr}->{$idx_ext}){next}
              $hits=&segQuery_hash($bed_idx->{$chr}->{$idx_ext},"$start_ext-$start");
              if(keys %{$hits} != 0){last}
            }
    
         }
       }elsif($strand eq "-"){
            ##inward  step 200bp, 1kb most, stop when hits found
             for my $i(1..5){
                my $start_ext = $start-$i*200;
                if($start_ext <= 0){$start_ext=1}
                my $idx_ext = int 0.5*($start_ext+$start)*1e-5;         
                if(!exists $bed_idx->{$chr}->{$idx_ext}){next}
                $hits=&segQuery_hash($bed_idx->{$chr}->{$idx_ext},"$start_ext-$start");
                if(keys %{$hits} != 0){last}
             } 
             ##outward step 500bp, 5kb most, stop when hits found 
             if(keys %{$hits} == 0){
                for my $i(1..10){
                  my $end_ext = $start+$i*500;  
                  my $idx_ext = int 0.5*($end_ext+$start)*1e-5;
                  if(!exists $bed_idx->{$chr}->{$idx_ext}){next}
                  $hits=&segQuery_hash($bed_idx->{$chr}->{$idx_ext},"$start-$end_ext");
                  if(keys %{$hits} != 0){last}
                }
      
             }  
         }else{die "unknown strand type $strand"}
   }#in situ not found then extend inward and outward in a stepwise manner

   ##report
   if(keys %{$hits} != 0){
     foreach(keys %{$hits}){ 
       my @temp=split/\t/,$hits->{$_};
       $temp[3]=$geneID."_".$temp[3];
       print join("\t",@temp),"\n";
     }
   }else{print LOST "not found for $geneID\n"}

}#while end

print STDERR "\n";
close LOST;



#########sub#

sub readBed(){
   my $file=shift;
   my %bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my @box=split/[\t ]+/,$_;
     #die "trunct bed file $file at $_\n" if(scalar @box < 6);
     die "start >= end at @box" if($box[1] >= $box[2]);
     if(!exists $bed{$box[0]}->{$box[3]}){$bed{$box[0]}->{$box[3]}=join"\t",@box   }else{die "dup $box[3] at $_\n"}
   }
   close BED;
   return \%bed;
}

sub readBed_segTree(){
   my $file=shift;
   my %bed;
   open BED, $file or die "$!";
   while(<BED>){
     chomp;
     next if($_ eq "" || $_=~/^#/ || $_=~/^\s+$/);
     my @box=split/[\t ]+/,$_;
     die "start >= end at @box" if($box[1] >= $box[2]);
     my $idx=int 0.5*($box[1]+$box[2])*1e-5; #100kb range indexed
     if(!exists $bed{$box[0]}->{$idx}->{$box[3]}){$bed{$box[0]}->{$idx}->{$box[3]}=join"\t",@box   }else{die "dup $box[3] at $_\n"}
   }
   close BED;
   return \%bed;
}




sub segQuery_hash(){  #query a given segment in a table, report all overlap ones
   my ($table,$seg)=@_;
   my %segs_hit;
   my ($s,$e)=split/-/,$seg;
   die "$s >= $e at $seg" if($s >= $e);
   foreach my $id(sort keys %{$table}){
     my (undef,$start,$end)=split/\t/,$table->{$id};
=pod
     if($start >= $s && $end <= $e){
       $segs_hit{$id}=$table->{$id};
      }elsif($start < $s && $end > $s && $end <= $e){
        $segs_hit{$id}=$table->{$id}
       }elsif($end > $e && $start > $s && $start < $e){
          $segs_hit{$id}=$table->{$id}
         }elsif($start <= $s && $end >= $e){$segs_hit{$id}=$table->{$id}}
=cut
     if($e < $start || $end < $s){}else{$segs_hit{$id}=$table->{$id} }  #segment overlap 
  
   }
   return \%segs_hit;
}


