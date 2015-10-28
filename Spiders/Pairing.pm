#!/usr/bin/perl 

######### Pairing #################################################
#                                                             #
#       **************************************************    #  
#       **** Pairing                    	          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2014 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Pairing;

use warnings;
use lib 'lib';
use Spiders::Dta;
use Spiders::PeakSimilarity;
use Statistics::R;
use Storable;
use Spiders::MathUtils;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.0;

@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
    };
    bless $self, $class;
    return $self;
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



sub set_CC_ratio
{
	my ($self,$CC_ratio)=@_;
	$self->{'_CC_ratio'}=$CC_ratio;
}

sub get_CC_ratio
{
	my $self=shift;
	return $self->{'_CC_ratio'};
}

sub set_NC_ratio
{
	my ($self,$NC_ratio)=@_;
	$self->{'_NC_ratio'}=$NC_ratio;
}

sub get_NC_ratio
{
	my $self=shift;
	return $self->{'_NC_ratio'};
}

sub set_NC_std
{
	my ($self,$NC_std)=@_;
	$self->{'_NC_std'}=$NC_std;
}

sub get_NC_std
{
	my ($self)=@_;
	return $self->{'_NC_std'};
}


sub set_CC_std
{
	my ($self,$CC_std)=@_;
	$self->{'_CC_std'}=$CC_std;
}

sub get_CC_std
{
	my ($self)=@_;
	return $self->{'_CC_std'};
}

sub set_NC_defect_loc
{
	my ($self,$NC_defect_loc)=@_;
	$self->{'_NC_defect_loc'}=$NC_defect_loc;
}

sub get_NC_defect_loc
{
	my ($self)=@_;
	return $self->{'_NC_defect_loc'};
}

sub set_CC_defect_loc
{
	my ($self,$CC_defect_loc)=@_;
	$self->{'_CC_defect_loc'}=$CC_defect_loc;
}

sub get_CC_defect_loc
{
	my ($self)=@_;
	return $self->{'_CC_defect_loc'};
}


sub set_NC_defect_scale
{
	my ($self,$NC_defect_scale)=@_;
	$self->{'_NC_defect_scale'}=$NC_defect_scale;
}

sub get_NC_defect_scale
{
	my ($self)=@_;
	return $self->{'_NC_defect_scale'};
}

sub set_CC_defect_scale
{
	my ($self,$CC_defect_scale)=@_;
	$self->{'_CC_defect_scale'}=$CC_defect_scale;
}	

sub get_CC_defect_scale
{
	my ($self)=@_;
	return $self->{'_CC_defect_scale'};
}	

sub clustering
{
	my $self=shift;
	my $params = $self->get_parameter();
	
	my $neutron = 1.00335;
	my $tolerance_ppm = $params->{cluster_tolerance};

	
	my ($cluster,$mz_hash,$select_mz) = @_;


# search the following peaks 

	my $search_loop = 1;
	my $flag=1;
	my $previous_int = $mz_hash->{$select_mz};
	my $selected_int = $mz_hash->{$select_mz};
	
	
# search backward C13 case
	while($search_loop && $flag)
	{	
		my $tolerance = $tolerance_ppm*$select_mz/1000000;
		my ($peak_mz_low,$peak_mz_high) =  ($select_mz - ($neutron)*$search_loop-$tolerance, $select_mz - ($neutron)*$search_loop +$tolerance);
		$flag = 0;
		foreach my $mz (keys %$mz_hash)
		{
			if($mz>$peak_mz_low && $mz<$peak_mz_high)
			{

				$flag=1;

				if(($mz_hash->{$mz} > $previous_int))  # C12
				{

					$cluster->{$select_mz}->{-$search_loop}->{'mz'} = $select_mz;
					$cluster->{$select_mz}->{-$search_loop}->{'int'} += $mz_hash->{$mz};
										
					$search_loop++;
#					$previous_int = $mz_hash->{$mz};
					
################## Always keep the strongest peaks ############
## especially true for C13 
					#delete $mz_hash->{$mz};
					delete $mz_hash->{$select_mz};	
##############################################						

				}
				elsif(($mz_hash->{$mz} < $previous_int))  # C13
				{

					$cluster->{$select_mz}->{-$search_loop}->{'mz'} = $select_mz;
					$cluster->{$select_mz}->{-$search_loop}->{'int'} += $mz_hash->{$mz};
										
					$search_loop++;
#					$previous_int = $mz_hash->{$mz};
					
################## Always keep the strongest peaks ############
## especially true for C13 
					delete $mz_hash->{$mz};
					#delete $mz_hash->{$select_mz};	
##############################################						

				}											
				else
				{
					$search_loop=0;				
				}
			}
		}
	}
	

# searching forward 		
	$flag=1;
	$search_loop = 1;

	while($search_loop && $flag)
	{	
		my $tolerance = $tolerance_ppm*$select_mz/1000000;
		my ($peak_mz_low,$peak_mz_high) =  ($select_mz + ($neutron)*$search_loop-$tolerance, $select_mz + ($neutron)*$search_loop +$tolerance);

		$flag = 0;
		foreach my $mz (keys %$mz_hash)
		{
			if($mz>$peak_mz_low && $mz<$peak_mz_high)
			{

				$flag=1;
#				print $select_mz,"\t",$mz,"\t",$mz_hash->{$mz},"\t",$previous_int,"\t00\n";	
				if($mz_hash->{$mz} < $previous_int) # C12
				{
#					print $select_mz,"\t",$mz,"\tcc\n";			
					$cluster->{$select_mz}->{$search_loop}->{'mz'} = $mz;
					$cluster->{$select_mz}->{$search_loop}->{'int'} += $mz_hash->{$mz};									
					$search_loop++;
					#$previous_int = $mz_hash->{$mz};

					delete $mz_hash->{$mz};
				}
				elsif($mz_hash->{$mz} > $previous_int) ## C13
				{
#						print $select_mz,"\t",$mz,"\t",$mz_hash->{$mz},"\t",$previous_int,"\tdd\n";	
					$cluster->{$select_mz}->{$search_loop}->{'mz'} = $mz;
					$cluster->{$select_mz}->{$search_loop}->{'int'} += $mz_hash->{$mz};									
					$search_loop++;
					#$previous_int = $mz_hash->{$mz};

					delete $mz_hash->{$select_mz};
				}
				
				else
				{
					$search_loop=0;
				}
			}
		
		}
	}
}


sub select_mono_peaks
{
	my ($self,$cluster) =@_;
	my $mono_cluster;

	foreach my $select_mz (keys %$cluster)
	{
## it must have isotopes
#		print $select_mz,"\t",scalar keys (%{$cluster->{$select_mz}}),"\n";
		#next if((scalar keys (%{$cluster->{$select_mz}}))<2);
		foreach my $mono_order (reverse sort {$a<=>$b} keys %{$cluster->{$select_mz}})
		{
			#$mono_cluster->{$select_mz}->{'mz'} = $cluster->{$select_mz}->{$mono_order}->{'mz'};
			$mono_cluster->{$select_mz}->{'mz'} = $cluster->{$select_mz}->{0}->{'mz'};
			$mono_cluster->{$select_mz}->{'int'} = $cluster->{$select_mz}->{0}->{'int'};
			
			#last;

		}		
	}
	return ($mono_cluster);
}

sub init_pairing
{
	my ($self,$mono_cluster,$dir) =@_;
	my $pairing_results;
	my $neutron = 1.00335;
	my $C_N_diff = 0.99703;

	my $params = $self->get_parameter();
	my $tolerance_C = $params->{'c12_c13_tolerance'};
	my $tolerance_N = $params->{'c12_n15_tolerance'}; 
	my $min_pair_correlation = $params->{min_pair_correlation};
	my $relative_isotopes_intensity = $params->{relative_isotopes_intensity};
	

	foreach my $mz1 (reverse sort {$mono_cluster->{$a}->{'int'}<=>$mono_cluster->{$b}->{'int'}} keys %$mono_cluster)
	{
		my @C_array=();
		my @N_array=();
		my @comb_N_C=();

		foreach my $mz2 (reverse sort {$mono_cluster->{$a}->{'int'}<=>$mono_cluster->{$b}->{'int'}} keys %$mono_cluster)
		{

			next if($mz2<=($mz1+0.5));

			next if($mz2>2*$mz1);

			my $diff = abs($mz2 - $mz1);
			next if (int($diff+0.1)==0);
			my $N_or_C = $diff/int($diff+0.1);
			my $max_C_num = $mz1 / 12;

			next if($diff>$max_C_num);
					
			next if($mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'} < $relative_isotopes_intensity or $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'} < $relative_isotopes_intensity);

			if($N_or_C<($neutron+$tolerance_C) and $N_or_C>($neutron-$tolerance_C))
			{

				$paired->{$mz1}->{'C'} = int($diff);				
				my $min_relative_int = $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'};
				push(@C_array,"C$paired->{$mz1}->{'C'}:$mz2|$N_or_C|$min_relative_int");
				$diff = 0;
 
			} 			
			elsif($N_or_C<($C_N_diff+$tolerance_N) and $N_or_C>($C_N_diff-$tolerance_N))
			{	

				$paired->{$mz1}->{'N'} = int($diff+0.1);
				my $min_relative_int = $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'};
				push(@N_array,"N$paired->{$mz1}->{'N'}:$mz2|$N_or_C|$min_relative_int");						
				$diff = 0;	
			} 			
		}
				
### generate all possible combination of C and N
		my %temp_score;
		for $carbon (@C_array) 
		{
		  my ($c_formula,$c_diff_score, $c_int_score, $c_similarity) = split(/\|/,$carbon);
		  for $nitrogen (@N_array) 
		  {
			my $temp = "$carbon,$nitrogen";
			push(@comb_N_C,$temp_score{$score});			
		  }
		}
		
		$pairing_results->{$mz1}->{'number'} = \@comb_N_C;
		$pairing_results->{$mz1}->{'intensity'} = $mono_cluster->{$mz1}->{'int'};		
	}
	
	return $pairing_results;
}

sub pairing
{
	my ($self,$mono_cluster,$dir) =@_;
	my $pairing_results;
	my $neutron = 1.00335;
	my $C_N_diff = 0.99703;
	#my $tolerance_ppm = 5;
######## +/-10 times of SD (see table 2) for tolerance	
#	my $tolerance_C = 0.001;
#	my $tolerance_N = 0.0045; 

	my $params = $self->get_parameter();
	my $tolerance_C = $params->{'c12_c13_tolerance'};
	my $tolerance_N = $params->{'c12_n15_tolerance'}; 
	my $min_pair_correlation = $params->{min_pair_correlation};
	my $relative_isotopes_intensity = $params->{relative_isotopes_intensity};
	
	my $peaksim = new Spiders::PeakSimilarity;
#	my %used_peaks;
    my $R = Statistics::R->new(r_bin=>"/home/xwang4/software/R-3.1.0/bin/R");

	foreach my $mz1 (reverse sort {$mono_cluster->{$a}->{'int'}<=>$mono_cluster->{$b}->{'int'}} keys %$mono_cluster)
	{
		my @C_array=();
		my @N_array=();
		my @comb_N_C=();
## skip the peak if the peak has been used for once		
		#next if(defined($used_peaks{$mz1}));
		

		#my $tolerance = $tolerance_ppm*$mz1/1000000;
		foreach my $mz2 (reverse sort {$mono_cluster->{$a}->{'int'}<=>$mono_cluster->{$b}->{'int'}} keys %$mono_cluster)
		{
## skip the peak if the peak has been used for once		
#			next if(defined($used_peaks{$mz2}));


			next if($mz2<=($mz1+0.5));

## limit the size the of the second peak			
			next if($mz2>2*$mz1);

			my $diff = abs($mz2 - $mz1);
			next if (int($diff+0.1)==0);
			my $N_or_C = $diff/int($diff+0.1);
### limit the C number			
			my $max_C_num = $mz1 / 12;

			next if($diff>$max_C_num);
					
			next if($mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'} < $relative_isotopes_intensity or $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'} < $relative_isotopes_intensity);

			if($N_or_C<($neutron+$tolerance_C) and $N_or_C>($neutron-$tolerance_C))
			{
			#	$used_peaks{$mz1} = 1;
			#	$used_peaks{$mz2}=1;
				
				$paired->{$mz1}->{'C'} = int($diff);


				my $mz1_peaks =retrieve("$dir/$mz1");
				my $mz2_peaks =retrieve("$dir/$mz2");
				
				my $mz1_2D = $peaksim->get_2D_peaks($mz1,$mz1_peaks);
				my $mz2_2D = $peaksim->get_2D_peaks($mz2,$mz2_peaks);
				my $peak_similarity = 0;
				$peak_similarity = $peaksim->peak_similarity($mz1_2D,$mz2_2D);

			#	$peak_similarity =~ /^\d+$/
				next if(!(defined($peak_similarity)) and $peak_similarity eq '');				
			#	next if($peak_similarity<$min_pair_correlation);				

				
				my $min_relative_int = 1;
#				if($mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'} < $min_relative_int)
#				{
#					$min_relative_int = $mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'};
#				}
#				else
#				{
					$min_relative_int = $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'};
#				}				
				my ($diff_score,$int_score) = $self->pair_scoring($N_or_C,$min_relative_int,"CC");
				push(@C_array,"C$paired->{$mz1}->{'C'}:$mz2|$diff_score|$int_score|$peak_similarity");
#				print $mz1,"\t",$mz2,"\t",$diff,"\t",$diff/int($diff+0.5),"\t",$paired->{$mz1}->{'C'},"C\t",$peak_similarity,"\t",$score,"\n";
			
				$diff = 0;
 
			} 
			
			
			
			
			elsif($N_or_C<($C_N_diff+$tolerance_N) and $N_or_C>($C_N_diff-$tolerance_N))
			{	
			#	$used_peaks{$mz1} = 1;
			#	$used_peaks{$mz2}=1;
				$paired->{$mz1}->{'N'} = int($diff+0.1);


				my $mz1_peaks =retrieve("$dir/$mz1");
				my $mz2_peaks =retrieve("$dir/$mz2");
				
				my $mz1_2D = $peaksim->get_2D_peaks($mz1,$mz1_peaks);
				my $mz2_2D = $peaksim->get_2D_peaks($mz2,$mz2_peaks);
				my $peak_similarity = 0;
				$peak_similarity = $peaksim->peak_similarity($mz1_2D,$mz2_2D);
				
				next if(!defined($peak_similarity));
#				next if($peak_similarity<$min_pair_correlation);

				
				my $min_relative_int = 1;
#				if($mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'} < $min_relative_int)
#				{
#					$min_relative_int = $mono_cluster->{$mz1}->{'int'}/$mono_cluster->{$mz2}->{'int'};
#				}
#				else
#				{
					$min_relative_int = $mono_cluster->{$mz2}->{'int'}/$mono_cluster->{$mz1}->{'int'};
#				}				
				
				my ($diff_score,$int_score) = $self->pair_scoring($N_or_C,$min_relative_int,"NC");
				push(@N_array,"N$paired->{$mz1}->{'N'}:$mz2|$diff_score|$int_score|$peak_similarity");				
	#			print $mz1,"\t",$mz2,"\t",$diff,"\t",$diff/int($diff+0.5),"\t",$paired->{$mz1}->{'N'},"N\t",$peak_similarity,"\t",$score ,"\n";			
				$diff = 0;	
				
			} 			
		}
				
### generate all possible combination of C and N
		my %temp_score;
		for $carbon (@C_array) 
		{
		  my ($c_formula,$c_diff_score, $c_int_score, $c_similarity) = split(/\|/,$carbon);
		  for $nitrogen (@N_array) 
		  {
			#print "$number,$letter\n";
			my ($n_formula,$n_diff_score, $n_int_score, $n_similarity) = split(/\|/,$nitrogen);
			my $total_score = ($c_diff_score + $c_int_score + $c_similarity + $n_diff_score + $n_int_score + $n_similarity);

			$R->set('x', $total_score);				
	
			$R->run(q`y<-pchisq(2*x,12,lower.tail=FALSE)`);
			my $final_p = $R->get('y');
			my $final_score = -log($final_p)/log(10);
			print $c_formula,"\t",$n_formula,"\t",$total_score,"\t",$final_score,"\n";
			$temp_score{$final_p} = "$carbon,$nitrogen,$final_score";
			#push (@comb_N_C, "$carbon,$nitrogen|");
		  }
		}

		foreach my $score (sort {$a<=>$b} keys %temp_score)
		{
			push(@comb_N_C,$temp_score{$score});
			last;
		}			
		
		$pairing_results->{$mz1}->{'number'} = \@comb_N_C;
		$pairing_results->{$mz1}->{'intensity'} = $mono_cluster->{$mz1}->{'int'};		
	}
	
	return $pairing_results;
}

sub pair_scoring
{
	my ($self,$diff,$relative_int,$NC_or_CC) = @_;

	my $nc_ratio = $self->get_NC_ratio;	
	my $cc_ratio = $self->get_CC_ratio;		
	my $NC_std=$self->get_NC_std();
	my $CC_std=$self->get_CC_std();
	
	my $NC_defect_loc=$self->get_NC_defect_loc();
	my $CC_defect_loc=$self->get_CC_defect_loc();
	my $NC_defect_scale=$self->get_NC_defect_scale();
	my $CC_defect_scale=$self->get_CC_defect_scale();
	
#	my $percentage_diff = $diff / $tolerance;
    my $R = Statistics::R->new(r_bin=>"/home/xwang4/software/R-3.1.0/bin/R");
	my $diff_p = 0.5;
	my $rel_int_p = 0.5;
	
	$relative_int=log($relative_int)/log(2);	
	$R->set('defect_diff', $diff);
	$R->set('rel_int', $relative_int);
	$R->set('CC_ratio', $cc_ratio);
	$R->set('NC_ratio', $nc_ratio);		
	$R->set('NC_ratio_std', $NC_std);
	$R->set('CC_ratio_std', $CC_std);
	$R->set('NC_defect_loc', $NC_defect_loc);
	$R->set('CC_defect_loc', $CC_defect_loc);	
	$R->set('NC_defect_scale', $NC_defect_scale);
	$R->set('CC_defect_scale', $CC_defect_scale);

	my $defect;
	if($NC_or_CC eq "NC")
	{
		$defect = 0.99703;
	}
	elsif($NC_or_CC eq "CC")
	{
		$defect = 1.00335;
	}
	else
	{
		print "please set a right defect\n";
		exit;
	}
	if($NC_or_CC eq "NC")
	{	
			$R->run(q`y<-pcauchy(defect_diff,NC_defect_loc,NC_defect_scale,FALSE)`);
			$R->run(q`y<-min(y,1-y)`);
			$diff_p = 1-2*$R->get('y');
	}
	elsif($NC_or_CC eq "CC")
	{	
			$R->run(q`y<-pcauchy(defect_diff,CC_defect_loc,CC_defect_scale,TRUE)`);
			$R->run(q`y<-min(y,1-y)`);
			$diff_p = 1-2*$R->get('y');	
	}

	if($NC_or_CC eq "NC")
	{		
			$R->run(q`y<-pnorm(rel_int,NC_ratio,NC_ratio_std,FALSE)`);
			$R->run(q`y<-min(y,1-y)`);			
			$rel_int_p = 1-2*$R->get('y');	
	}
	elsif($NC_or_CC eq "CC")
	{			
			$R->run(q`y<-pnorm(rel_int,CC_ratio,CC_ratio_std,FALSE)`);
			$R->run(q`y<-min(y,1-y)`);			
			$rel_int_p = 1-2*$R->get('y');	
	}	

=head	
	if($peak_similarity==0)
	{
		$peak_similarity=1;
	}
=cut	
	if($diff_p==0)
	{
		$diff_p=1;
	}
	if($rel_int_p==0)
	{
		$rel_int_p=1;
	}

	my $diff_score = -log($diff_p)/log(10);
	my $int_score = -log($rel_int_p)/log(10);
#	my $sim_score = $peak_similarity;

	return ($diff_score,$int_score);
}


1;

		



