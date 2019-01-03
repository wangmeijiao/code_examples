
   cat chr03.japo-VS-zs97-VS-glab-VS-meri.mfa |perl mfaView_chrbig.pl  > chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.svg && display chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.svg



     bigWigToBedGraph /home/mjwang/pwdexx/oryza_epiCompara_sixSpecies/tss_gBGC_tigr6/ancestral_infer/ancestral_alignment/chr03/pecan_ortheus_giventree_four/pecan_out/post_process/conserveScore/phylofit/../../../post_process/conserveScore/phyloP/chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.phastP.gerp.bw chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.phastP.gerp.bedgraph

    cat chr03.japo-VS-zs97-VS-glab-VS-meri.refDegap.mfa |perl mfaView_chrbig.pl chr03.japo-VS-zs97-VS-glab-VS-meri.mfa.phastP.gerp.bedgraph > chr03.japo-VS-zs97-VS-glab-VS-meri.refDegap.mfa.svg
  
    convert chr03.japo-VS-zs97-VS-glab-VS-meri.refDegap.mfa.svg chr03.japo-VS-zs97-VS-glab-VS-meri.refDegap.mfa.svg.png #small png file  
