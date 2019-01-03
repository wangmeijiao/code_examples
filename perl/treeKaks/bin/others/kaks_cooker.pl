#!/usr/bin/perl

# Compute number of synonymous and nonsynonymous sites, synonymous and nonsynonymous mutations from pairwise-aligned nucleotide sequences
# By Yi Xing, UCLA May 2004

# Use R=0.5 ( transition vs transversion = 1 : 2)



%codon_table= ("TTT",'F',"TTC",'F',"TTA",'L',"TTG",'L',"CTT",'L',"CTC",'L',"CTA",'L',"CTG",'L',"ATT",'I',"ATC",'I',"ATA",'I',"ATG",'M',"GTT",'V',"GTC",'V',"GTA",'V',"GTG",'V',"TCT",'S',"TCC",'S',"TCA",'S',"TCG",'S',"CCT",'P',"CCC",'P',"CCA",'P',"CCG",'P',"ACT",'T',"ACC",'T',"ACA",'T',"ACG",'T',"GCT",'A',"GCC",'A',"GCA",'A',"GCG",'A',"TAT",'Y',"TAC",'Y',"TAA",'*',"TAG",'*',"CAT",'H',"CAC",'H',"CAA",'Q',"CAG",'Q',"AAT",'N',"AAC",'N',"AAA",'K',"AAG",'K',"GAT",'D',"GAC",'D',"GAA",'E',"GAG",'E',"TGT",'C',"TGC",'C',"TGA",'*',"TGG",'W',"CGT",'R',"CGC",'R',"CGA",'R',"CGG",'R',"AGT",'S',"AGC",'S',"AGA",'R',"AGG",'R',"GGT",'G',"GGC",'G',"GGA",'G',"GGG",'G');





=cut

# compute number of synonymous and nonsynonymous site per codon

foreach $codon (keys (%codon_table)) {
    $aa=$codon_table{$codon};

    $syn_site=0;
    $nosyn_site=0;
    @nt=('A','G','C','T');
    for $m(0..2) {
        $codons[$m]=substr $codon,$m,1;
    }
    for $m(0..2) {
        for $n(0..3) {
            $#new_codon=-1;
            @new_codons=@codons;
            if ($nt[$n] ne $codons[$m]) {
                $new_codons[$m]=$nt[$n];
            
                $new_codon=$new_codons[0].$new_codons[1].$new_codons[2];

                if ($codon_table{$new_codon} eq $aa) {
                    $syn_site++;
                }
                else {
                    $nosyn_site++;
                }
            }
    
        }
    }
    print $codon,"\t",$syn_site,"\t",$nosyn_site,"\n";

}

=cut

# read codon synonymous & nonsynonymous sites info
    
open(In,"codon_info") || die " can not open codon syn vs nonsyn sites table!\n";

while(<In>) {
    chomp $_;
    @codon_info=split(/\t/);
    
    $syn{$codon_info[0]}=$codon_info[1];
    $nosyn{$codon_info[0]}=$codon_info[2];
}

close In;


# Read two sequences from FASTA file
open(FASTA,"$ARGV[0]") || die "can not open FASTA file for two sequences!\n";

$n_of_seq=-1;
while(<FASTA>) {
    if (/>/) {
        if ($n_of_seq >=0) {
            $sequence[$n_of_seq]=$seq;
        }
        $n_of_seq+=1;
        chomp $_;
        $id[$n_of_seq]=substr $_,1;
        $seq="";
    }
    else {
        chomp $_;
        $seq=$seq.$_;
    }
}

$sequence[$n_of_seq]=$seq;

close FASTA;




$n_of_nt=length $sequence[0] ;
$n_of_codons=$n_of_nt / 3;



$Ls=0;
$Ln=0;
$Ss=0;
$Sn=0;


# Use sequence 1 as template
for $i(1..$n_of_codons) {
    $current_codon=substr $sequence[0],($i-1)*3,3;
    
    for $k(0..2) {
        $current_nt[$k]=substr $current_codon,$k,1
        }
    $other_codon=substr $sequence[1],($i-1)*3,3;

    
    for $k(0..2) {
        $other_nt[$k]=substr $other_codon,$k,1;
    }

    # Compute # of synonymous and nonsynonymous sites
    
    $Ls=$Ls+$syn{$current_codon};
    $Ln=$Ln+$nosyn{$current_codon};


    # Compute # of synonymous and nonsynonymous mutations
    for $k(0..2) {
        if ( $current_nt[$k] ne $other_nt[$k]) {
            $new_codon="";
            for $p(0..2) {
                if ($p == $k) {
                    $new_codon=$new_codon.$other_nt[$p];
                }  
                else {
                    $new_codon=$new_codon.$current_nt[$p];
                }
            }

            if ($codon_table{$current_codon} ne $codon_table{$new_codon}) {
                $Sn++;
            }
            else {
                $Ss++;
            }
        }
    }
}


# switch order , use sequence 2 as template
for $i(1..$n_of_codons) {
    $current_codon=substr $sequence[1],($i-1)*3,3;
    
    for $k(0..2) {
        $current_nt[$k]=substr $current_codon,$k,1
        }
    $other_codon=substr $sequence[0],($i-1)*3,3;
    
    
    for $k(0..2) {
        $other_nt[$k]=substr $other_codon,$k,1;
    }
    
    
    $Ls=$Ls+$syn{$current_codon};
    $Ln=$Ln+$nosyn{$current_codon};


    
    for $k(0..2) {
        if ( $current_nt[$k] ne $other_nt[$k]) {
            $new_codon="";
            for $p(0..2) {
                if ($p == $k) {
                    $new_codon=$new_codon.$other_nt[$p];
                }  
                else {
                    $new_codon=$new_codon.$current_nt[$p];
                }
            }
        
            if ($codon_table{$current_codon} ne $codon_table{$new_codon}) {
                $Sn++;
            }
            else {
                $Ss++;
            }
        }
    }
}

$Sn=$Sn/2;
$Ln=$Ln/2;
$Ss=$Ss/2;
$Ls=$Ls/2;


$ka=$Sn/$Ln;
$ks=$Ss/$Ls;

print $Sn,"\t",$Ln,"\t",$ka,"\t",$Ss,"\t",$Ls,"\t",$ks,"\t",$ka/$ks,"\n";





