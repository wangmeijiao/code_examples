

sub lineFit(){
 #result resemble to R lm()
 my ($x,$y) = @_;
 my $lineFit = Statistics::LineFit->new();
 $lineFit->setData ($x, $y) or die "Invalid data";
 my ($intercept, $slope) = $lineFit->coefficients();
 defined $intercept or die "Can't fit line if x values are all equal";

 #print "$intercept $slope\n";
 return ($intercept,$slope);
=pod
 my $rSquared = $lineFit->rSquared();
 my $meanSquaredError = $lineFit->meanSqError();
 my $durbinWatson = $lineFit->durbinWatson();
 my $sigma = $lineFit->sigma();
 my ($tStatIntercept, $tStatSlope) = $lineFit->tStatistics();
 my @predictedYs = $lineFit->predictedYs();
 my @residuals = $lineFit->residuals();
 my ($varianceIntercept, $varianceSlope) = $lineFit->varianceOfEstimates();
=cut


}

