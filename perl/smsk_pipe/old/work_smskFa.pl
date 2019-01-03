use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
##use repeatMasker softmask only
## $prefix.fa -> gffDir + $prefix.smsk.fa
##usage: perl thisfile.pl -fa in.fa -species rice -keep
#        perl work_smskFa.pl -fa H1_brac.fa -lib brac

my ($fa,$species,$lib,$keep);
GetOptions("fa=s",\$fa,"species:s",\$species,"lib:s",\$lib,"keep!",\$keep);
my $prefix=`basename $fa .fa`;
chomp $prefix;

my %lib=( #for custom libs
          "brac"=>"/home/mjwang/progs/repeatMasker4.0.3/Libraries/customed_lib_ff/TElib_ff.lib",
          "tigr_repeatv3.3" => "/home/mjwang/progs/repeatMasker4.0.3/Libraries/customed_lib_tigrrepeat/tigr_oryza_repeatsv3.3.lib",

);

#1,repeatmasker use species or custom lib, only one per time
  my $outdir="out_".$prefix."-smsk";
  mkdir $outdir;
  #in fact, you can let RepeatMasker do some trf job if remove -nolow . Perhaps it will not do as good as trfBig ?
   my $CMD="RepeatMasker -a -para 14 -dir $outdir -no_is -xsmall -engine wublast "; 
   #my $CMD="RepeatMasker -a -para 14 -dir $outdir -no_is -xsmall -nolow -engine wublast ";
   #my $CMD="RepeatMasker -a  -para 14 -dir $outdir -no_is -xsmall -norna -nolow -engine wublast ";

  if(defined $species && !defined $lib ){
        $CMD.="-species '$species' "." $fa "; #print $CMD,"\n";
        print "start to run repeatMasker with build in lib #$species# ..\n";
        system($CMD);
      }elsif( !defined $species && defined $lib){
          $CMD.="-lib $lib{$lib} "." $fa "; #print $CMD,"\n";
          print "start to run repeatMasker with custom lib #$lib# ..\n";
          system($CMD);
         }elsif( !defined $species && !defined $lib ){
              die "err: empty paras\n"
             }elsif(defined $species && defined $lib){die "err: species and lib confused\n"}

#2, repeat2gff/bed
print "repeatMasker.out -> gff/bed.. \n";

`bash rpmkout2bedgff_filter.sh $prefix bed`;
`bash rpmkout2bedgff_filter.sh $prefix gff`;

#`cat $outdir/${prefix}.fa.out |perl bin/rpmk2gffv2.1.pl -filter 0.3 -bed -prefix $prefix -outdir bed_${prefix} 2> out2bed.err`;
#`cat bed_${prefix}/${prefix}_*.bed > bed_${prefix}/${prefix}.fa.repeats.bed`; # not include lost.bed
#`mv  out2bed.err bed_${prefix}`;
#`mv  ${prefix}_lost.bed bed_${prefix}`;

#`cat $outdir/${prefix}.fa.out |perl bin/rpmk2gffv2.1.pl -filter 0.3 -prefix $prefix -outdir gff_${prefix} 2> out2gff.err`;
#`cat gff_${prefix}/${prefix}_*.gff > gff_${prefix}/${prefix}.fa.repeats.gff`; # not include lost.gff
#`mv   out2gff.err gff_${prefix}`;
#`mv  ${prefix}_lost.gff gff_${prefix}`;



#3,cleanup
print "clean up ..";
`mv $outdir/${prefix}.fa.out .`;
`mv $outdir/${prefix}.fa.masked ${prefix}.smsk.fa `;
if($keep){}else{`rm -r $outdir`}

print "all done.\n";

