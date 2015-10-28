package Spiders::PeakSimilarity;

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Search database	                      ****    #
#       ****					                      ****    #
#       ****Copyright (C) 2014 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

use Storable;
use Statistics::Basic qw(:all nofill);
use Statistics::Distributions; 

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(findTag readTags exportMatches);
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.01;


sub new{
	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self
}

sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}

sub get_2D_peaks
{
	my ($self,$mass_3D,$peaks) = @_;
	my $peak_int;
	my $RT = $peaks->{'2Drt'};
	my $int = $peaks->{'2Dintensity'};

	for(my $i=0;$i<=$#$RT;$i++)
	{
		$peak_int->{$RT->[$i]} = $int->[$i];

	}

	return $peak_int;
}

sub peak_similarity
{
	my ($self,$peak1_int,$peak2_int)=@_;
	my $params = $self->get_parameter();
	
### find the missing peaks in the both C12 and C13
	my @peak1_array = ();
	foreach (sort {$a<=>$b} keys %$peak1_int) {
		if(exists $peak2_int->{$_})
		{
			next if(($peak2_int->{$_}/$peak1_int->{$_})>3 or ($peak1_int->{$_}/$peak2_int->{$_})>3); 	
			push(@peak1_array, $peak1_int->{$_});
		}
	}
	if($#peak1_array<2)
	{
		return 0;
	}
	my $per_peak1 = $#peak1_array / (scalar keys %$peak1_int);
	
	my @peak2_array = ();	
	foreach (sort {$a<=>$b} keys %$peak2_int) {
		if(exists $peak1_int->{$_})
		{
			next if(($peak2_int->{$_}/$peak1_int->{$_})>3 or ($peak1_int->{$_}/$peak2_int->{$_})>3); 
			push(@peak2_array, $peak2_int->{$_});		
		}
	}
	if($#peak2_array<2)
	{
		return 0;
	}	
	my $per_peak2 = $#peak2_array / (scalar keys %$peak2_int);

########  peak overlap 	
	if ($per_peak1< 0.5 or $per_peak2<0.5)
	{
		return 0;
	}
	my $corr=1;
	$corr = correlation(\@peak1_array,\@peak2_array);

	$corr = 1 if($corr eq 'n/a');
	if(defined($corr))
	{
		my $n_num = 1;
		if($#peak1_array>2)
		{
			$n_num = $#peak1_array-2;
		}
		$tstat=$corr/sqrt((1-$corr*$corr)/$n_num); 

		my $p_value= Statistics::Distributions::tprob($n_num,abs($tstat));
		if($p_value == 0)
		{
			return 30;
		}
		my $sim_score=-log($p_value)/log(10);
		return $sim_score;

	}
	else
	{
		return 0;
	}
}

1;
