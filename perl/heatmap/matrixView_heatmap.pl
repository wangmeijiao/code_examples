use strict;
use warnings;
use Getopt::Long;
use SVG;
use Data::Dumper;

#data matrix format:
#  a    b     c
#  0.2  0.12  0.9
#  0.89 7.1   1.9

# row order must be clustered or ordered beforehand (draw as is)
# header=T/F
# check and make sure that not have NA and empty values, all datapoint must be numeric (rownames col are not allowed)
# try to guess matrix nrow and ncol after readin the data and assign the best width x height

my($data,$header,$dimXY,$zlim,);
GetOptions("data=s"=>\$data,"header=s",\$header,"dimXY=s"=>\$dimXY,"zlim=s"=>\$zlim);
my ($width,$height)=split/,/,$dimXY; #not include x/y margin
my ($zmin,$zmax)=split/,/,$zlim;

my ($data_idx,$header_idx,$nrow,$ncol,$min,$max)=&readTable($data,$header);
#print Dumper $header_idx,$nrow,$ncol,$min,$max,$data_idx;
#exit;

my ($box_h,$box_w);
my $y_start=50;
my $edge=10;
$box_w=$width/$ncol;
$box_h=$height/$nrow;
my $x_left=$edge;
my $x_right=$edge+$width;
my @color=@{&greyMono256()};
#my @color=@{&redMono256()};
#my @color=@{&greenMono256()};
my $ncolor=$#color;

my $svg=SVG->new("width",$width+2*$edge+10*3,"height",$height+$y_start*2);

my $cnt_line;
my $y=$y_start;
foreach (@{$data_idx}){
 $cnt_line++;
 #if($cnt_line >2000){last}
 if($cnt_line%1000==0){print STDERR "#"}
 my @box=split/[\t ]+/,$_;
 my $x=0;
 foreach my $data(@box){   
   my $colN; #color index
   if($data<=$zmin){$colN=0}elsif($data>$zmin && $data<$zmax){$colN=&round ($ncolor*($data-$zmin)/($zmax-$zmin))}else{$colN=$ncolor} #cut and map datapoint to color index
   #$svg->rect("x",$x_left+$x,"y",$y,"width",$box_w,"height",$box_h,"fill",$color[$colN],"stroke","none","stroke-opacity",0.8,"stroke-width",0);
   $svg->rect("x",$x_left+$x,"y",$y,"width",$box_w,"height",$box_h,"fill",$color[$colN],"stroke",$color[$colN],"stroke-width",6);#stroke-width, the bigger, the shapper when deal with huge data. But not too big
   $x+=$box_w+10;
 }
 $y+=$box_h; 
}
print STDERR "\n";

#draw legends
my $text_x=$x_left+0.5*$box_w;
my $text_y=$y_start-10;
foreach(@{$header_idx}){
   $svg->text("x",$text_x,"y",$text_y,"width",$box_w,"height",10,"font-family","Arial", "font-weight","bold","text-anchor","middle","font-size",10, "-cdata", $_);
   $text_x+=$box_w+10;  
}

print $svg->xmlify();

#####subs######

sub readTable(){
   my $file=shift;
   my $head=shift;
   my @header;
   my @data;
   my $cnt;
   my $ncol;
   my ($min,$max);
   open DATA, $file or die "$!";
   while(<DATA>){
     chomp;
     next if($_ eq "" || $_ =~/^#/ || $_=~/^\s+$/);
     $cnt++;
     if($cnt==1){
       if($head eq "T"){
           @header=split/[\t ]+/,$_;
           $ncol=scalar @header;
           next;
       }elsif($head eq "F"){
           my @box=split/[\t ]+/,$_;
           $ncol=scalar @box;
           if(&notNA(\@box) == 0){die "has NA at $_"}
           if(&isNumeric(\@box) == 0){die "not numeric at $_"}
           push @data,$_;
           ($min,$max)=($box[0],$box[0]);
         }else{die "unknow header swith $head"}
     }
     my @box_line=split/[\t ]+/,$_;
     my $ncol_line=scalar @box_line;
     if(!defined $ncol || $ncol_line != $ncol){die "ncol err $ncol != $ncol_line at $_"}
     if(&notNA(\@box_line) == 0){die "has NA at $_"}
     if(&isNumeric(\@box_line) == 0){die "not numeric at $_"}
     push @data,$_;
     if(! defined $min && !defined $max){
       ($min,$max)=($box_line[0],$box_line[0])
     }elsif(defined $min && defined $max){
         foreach my $ctrl(@box_line){ 
           if($ctrl < $min){$min = $ctrl}
           if($ctrl > $max){$max = $ctrl}
         }
       }else{die "err min max $min , $max"}    
   }#while end
   close DATA;
   if($head eq "T"){$cnt-=1}
   return (\@data,\@header,$cnt,$ncol,$min,$max);

}

sub isNumeric(){
  my $idx=shift;
  my $flag=1;
  foreach my $ctrl (@{$idx}){
    #signed    
    $ctrl=~s/^-//;
    #decimals
     if($ctrl =~/[^0-9\.]/){$flag = 0}

    #scientific notation

  }
  return $flag;
}

sub notNA(){
    my $idx=shift;
    my $flag=1;
    foreach my $ctrl (@{$idx}){if($ctrl eq "NA" || $ctrl eq "NULL" || $ctrl eq "Na" || $ctrl eq "n/a" || $ctrl eq "N/A"){$flag=0}}
    return $flag;
}

sub redMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb(255,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub redMono9(){
   my @color=("#FFFFFF", "#FFDFDF", "#FFBFBF", "#FF9F9F" ,"#FF7F7F" ,"#FF5F5F", "#FF3F3F" ,"#FF1F1F", "#FF0000");
   return \@color;  
}

sub greenMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,255,$i)";
   }
   #print "@color\n";
   return \@color;
}

sub greyMono256(){
   my @color;
   for(my $i=255;$i>=0;$i--){
     push @color,"rgb($i,$i,$i)";
   }
   #print "@color\n";
   return \@color;
}


sub round(){
   my $n=shift;
   my $i=int $n;
   my $d=$n-$i;
   if($d>0.6){return $i+1}else{return $i}
}

sub minMax(){
   my @data=@{shift @_};
   my ($min,$max)=(0,0);
   foreach(@data){
     if($min>$_){$min=$_}
     if($max<$_){$max=$_}
   }
   return($min,$max);
}
