
#####color spectrum code start###
#1,full-spectrum visuable natural light (1019 colors, blue to yellow to red )
sub fullSpecColor(){
  my ( @color_arr, @arr );
  @arr= (0,0,128);
  for ( my $i=128; $i<255; $i++ ){
     push @color_arr, "rgb($arr[0],$arr[1],$arr[2])";
     $arr[2]++;
  }
  for ( my $i=0; $i<255; $i++ ){
     $arr[1]++;
     push @color_arr, "rgb($arr[0],$arr[1],$arr[2])";
  }
  for ( my $i=0; $i<255; $i++ ){
     $arr[0]++;
     $arr[2]--;
     push @color_arr, "rgb($arr[0],$arr[1],$arr[2])";
  }
  for ( my $i=0; $i<255; $i++ ){
     $arr[1]--;
     push @color_arr, "rgb($arr[0],$arr[1],$arr[2])";
  }
  for ( my $i=0; $i<127; $i++ ){
     $arr[0]--;
     push @color_arr, "rgb($arr[0],$arr[1],$arr[2])";
  }
  return \@color_arr;
}
#2, hot/cold two-color spectrum (510 colors, blue to white to red)
sub hotCold510(){
  my @color;
  #left to right side:  blue(0,0,255) -> white(255,255,255)
  for(my $i=0;$i<=255;$i++){
    push @color,"rgb($i,$i,255)";
  }
  #left to right side:  white(255,255,255) -> red(255,0,0)
  for(my $i=255;$i>=0;$i--){
    push @color,"rgb(255,$i,$i)";
  }
  #print "@color\n";
  return \@color;
}

#3,mono-color gradient spectrum (from colorbrewer2.org)
sub mono9(){
  #9-class Purples:
  #fcfbfd
  #efedf5
  #dadaeb
  #bcbddc
  #9e9ac8
  #807dba
  #6a51a3
  #54278f
  #3f007d
    
  #9-class Greys:
  #ffffff
  #f0f0f0
  #d9d9d9
  #bdbdbd
  #969696
  #737373
  #525252
  #252525
  #000000

  #9-class Greens
  #f7fcf5
  #e5f5e0
  #c7e9c0
  #a1d99b
  #74c476
  #41ab5d
  #238b45
  #006d2c
  #00441b
}


#########################

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



