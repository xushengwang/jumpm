#!/usr/bin/perl -I /home/xwang4/scripts/dev/

package Spiders::PeakDetection;
        
use strict;
use warnings;
use List::Util qw(min max);
use Statistics::Lite qw(mean median);
use Spiders::Params;
use File::Basename;
use MIME::Base64;
use Storable;
use vars qw($VERSION @ISA @EXPORT);


$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

#####################
## Read parameters ##
#####################

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}

sub Gen_3Dpeaks
{
	my ($self,$params,$mzXML,$dtaPath) = @_;
=head	
	my $paramFile = shift (@ARGV);
	my $params = Spiders::Params -> new('-path' => $paramFile);
	$params = $params -> parse_param();

	my $mzXML = shift (@ARGV);
=cut	
	my $inputFile = $mzXML;

	####################
	## Get MS spectra ##
	####################

	my %ms = getMsHash($inputFile, $$params{'first_scan_extraction'}, $$params{'last_scan_extraction'});
	print "\n  MS data conversion finished\n";

	###############################
	## Mass ranges for each scan ##
	###############################

	my %massRanges;
	my @minMass;
	my @maxMass;
	my $msCount = scalar(keys %ms);
	for (my $i = 0; $i < $msCount; $i++) {
		push (@minMass, min(@{$ms{$i}{'mass'}}));
		push (@maxMass, max(@{$ms{$i}{'mass'}}));
	}
	%massRanges = ("minMass" => \@minMass, "maxMass" => \@maxMass);

	###########################
	## MaxQuant main routine ##
	###########################

	my @peaks;

	my $peaks3DCount = -1;
	my @cache;
	my $oldMinInd = -1;
	my $oldMaxInd = -1;
	my $gap = $$params{'scanWindow'} - 1;
	if ($gap < 0) {
		die "\"scanWindow\" parameter should be greater than 0\n";
	}
	my $intensityThresholdScale2D = $$params{'intensityThresholdScale2D'};
	my $intensityThreshold2D = $$params{'intensityThreshold2D'};

	for (my $i = 0; $i < $msCount; $i++) {
		## Put former and latter scans to cache (temporary) array 
		my $minInd = max(0, $i - $gap - 1);
		my $maxInd = min($msCount - 1, $i + $gap + 1);
		if ($i == 0) {
			for (my $j = 0; $j <= $maxInd; $j++) {
				my $msCentroidList = DetectPeaks($$params{'numCentroidPoints'}, $intensityThresholdScale2D, $intensityThreshold2D, $ms{$j});
				$$msCentroidList{'scanNumber'} = $ms{$j}{'scanNumber'};
				$$msCentroidList{'scanIndex'} = $j;
				$$msCentroidList{'RT'} = $ms{$j}{'RT'};
				push @cache, $msCentroidList;
			}
		} else {
			for (my $j = $oldMinInd; $j < $minInd; $j++) {
				shift @cache;
			}
			for (my $j = $oldMaxInd + 1; $j <= $maxInd; $j++) {
				my $msCentroidList = DetectPeaks($$params{'numCentroidPoints'}, $intensityThresholdScale2D, $intensityThreshold2D, $ms{$j});
				$$msCentroidList{'scanNumber'} = $ms{$j}{'scanNumber'};
				$$msCentroidList{'scanIndex'} = $j;
				$$msCentroidList{'RT'} = $ms{$j}{'RT'};
				push @cache, $msCentroidList;
			}
		}
		
		##################################################################################
		## Reduction step								##
		## For each scan, take the MS centroid list and select (i.e. reduce) the peaks	##
		## which can form a 3D-peak with other peaks in neighboring scans		##
		##										##
		## Reduction is done in forward and backward direction				## 
		##################################################################################
		
		my %p = %{$cache[$i - $minInd]};	## A MS centroid list under consideration
		my $pCount = @{$p{'peakCenterMass'}};
		my @valids;
		my $count = 0;
		for (my $j = 0; $j < $pCount; $j++) {
			my $cm = $p{'peakCenterMass'}[$j];
			my $match = 0;
			my $ntries = 0;
			for (my $k = $i - 1; $k >= $minInd; $k--) {
				if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
					my %q = %{$cache[$k - $minInd]};
					my $ind;
					($match, $ind) = GetClosestValue(\%q, $cm, $$params{'matchPpm'});
					if ($match) {
						last;
					}
					$ntries++;
					if ($ntries > $$params{'scanWindow'}) {
						last;
					}
				}
			}
			if (!$match) {
				$ntries = 0;
				for (my $k = $i + 1; $k <= $maxInd; $k++) {
					if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
						my %q = %{$cache[$k - $minInd]};
						my $ind;
						($match, $ind) = GetClosestValue(\%q, $cm, $$params{'matchPpm'});
						if ($match) {
							last;
						}
						$ntries++;
						if ($ntries > $$params{'scanWindow'}) {
							last;
						}
					}
				}
			}
			if ($match) {
				push (@valids, $j);
			}
		}
		my %reduced = ExtractMsCentroidList(\%p, \@valids);
		
		##########################################################################################
		## Generation of 3D peaks by combining peaks within the specified mass tolerance range 	##
		## in the former scans (+/- scanWindow) 						##
		##											##
		## 3D-peak generation is done in only backward direction				##
		##########################################################################################

		$cache[$i - $minInd] = \%reduced;
		my $reducedCount = @{$reduced{'peakCenterMass'}};
		for (my $j = 0; $j < $reducedCount; $j++) {
			my $cm = $reduced{'peakCenterMass'}[$j];
			my $match = 0;
			my $ntries = 0;
			my @matchedPeaks3DInds;
			for (my $k = $i - 1; $k >= $minInd; $k--) {	## Search for previous scans stored in a cache array
				if ($cm >= ($massRanges{'minMass'}[$k]) && $cm <= ($massRanges{'maxMass'}[$k])) {
					my %q = %{$cache[$k - $minInd]};
					my ($matchIndicator, $ind) = GetClosestValue(\%q, $cm, $$params{'matchPpm'});
					## $matchIndicator = 1 means that the j-th (reduced) peak in the i-th scan
					## can form a 3D-peak with $ind-th (reduced) peak in the previous scan (%q)
					if ($matchIndicator) {
						push (@matchedPeaks3DInds, @{$q{'peaks3DIndex'}}[@{$ind}]);
						$match = 1;
					}
				}
			}
			if ($match) {
				@matchedPeaks3DInds = Unique(@matchedPeaks3DInds);
				my $peaks3DInd;
				if ($#matchedPeaks3DInds > 0) {	## Match to the peaks in mutliple previous scans
					$peaks3DInd = min(@matchedPeaks3DInds);
					foreach my $matchedPeaks3DInd (@matchedPeaks3DInds) {
						# Processing of singleton 3D peaks
						if ($matchedPeaks3DInd ne $peaks3DInd) {
							my $peakTobeRemoved = $peaks[$matchedPeaks3DInd];
							## comment out by xusheng ######################################################
							push (@{$peaks[$peaks3DInd]{'centerMz'}}, @{$$peakTobeRemoved{'centerMz'}});
							push (@{$peaks[$peaks3DInd]{'intensity'}}, @{$$peakTobeRemoved{'intensity'}});
							push (@{$peaks[$peaks3DInd]{'scanIndex'}}, @{$$peakTobeRemoved{'scanIndex'}});
							push (@{$peaks[$peaks3DInd]{'scanNumber'}}, @{$$peakTobeRemoved{'scanNumber'}});
							push (@{$peaks[$peaks3DInd]{'scanRT'}}, @{$$peakTobeRemoved{'scanRT'}});
							################################################################################
						
							# Revise the information in cache
							foreach my $scanInd (@{$peaks[$matchedPeaks3DInd]{'scanIndex'}}) {
								for (my $l = 0; $l <=$#cache; $l++) {
									if ($cache[$l]{'scanIndex'} == $scanInd) {
										for (my $m = 0; $m <= $#{$cache[$l]{'peaks3DIndex'}}; $m++) {
											if ($cache[$l]{'peaks3DIndex'}[$m] == $matchedPeaks3DInd) {
												$cache[$l]{'peaks3DIndex'}[$m] = $peaks3DInd;
											}
										}
									}
								}
							}
							delete $peaks[$matchedPeaks3DInd];
						}
					}
				} else {
					$peaks3DInd = shift @matchedPeaks3DInds;
				}
				push (@{$cache[$i - $minInd]{'peaks3DIndex'}}, $peaks3DInd);	# dummy value for keeping array size
				########################################################	
				push (@{$peaks[$peaks3DInd]{'centerMz'}}, $reduced{'peakCenterMass'}[$j]);
				push (@{$peaks[$peaks3DInd]{'minMz'}}, $reduced{'peakMinMass'}[$j]);
				push (@{$peaks[$peaks3DInd]{'maxMz'}}, $reduced{'peakMaxMass'}[$j]);
				push (@{$peaks[$peaks3DInd]{'intensity'}}, $reduced{'peakIntensity'}[$j]);
				push (@{$peaks[$peaks3DInd]{'scanIndex'}}, $reduced{'scanIndex'});
				#########################################################
				push (@{$peaks[$peaks3DInd]{'scanNumber'}}, $reduced{'scanNumber'});
				push (@{$peaks[$peaks3DInd]{'scanRT'}}, $reduced{'RT'});
			}
			if (!$match) {
				if ($i < $msCount) {
					$peaks3DCount++;
					push (@{$cache[$i - $minInd]{'peaks3DIndex'}}, $peaks3DCount);
					push (@{$peaks[$peaks3DCount]{'centerMz'}}, $reduced{'peakCenterMass'}[$j]);
					push (@{$peaks[$peaks3DCount]{'minMz'}}, $reduced{'peakMinMass'}[$j]);
					push (@{$peaks[$peaks3DCount]{'maxMz'}}, $reduced{'peakMaxMass'}[$j]);
					push (@{$peaks[$peaks3DCount]{'intensity'}}, $reduced{'peakIntensity'}[$j]);
					push (@{$peaks[$peaks3DCount]{'scanIndex'}}, $reduced{'scanIndex'});
					push (@{$peaks[$peaks3DCount]{'scanNumber'}}, $reduced{'scanNumber'});
					push (@{$peaks[$peaks3DCount]{'scanRT'}}, $reduced{'RT'});
				}
			}
		}	
		$oldMinInd = $minInd;
		$oldMaxInd = $maxInd;
		my $currentCount = $i + 1;
		print "  Generation of 3D-peaks: traversing scan#$currentCount out of $msCount\r";	
	}

	## Peak re-numbering (because some peaks (i.e. singleton 3D peaks) are 'undef')
	my @tempPeaks;
	for (my $i = 0; $i <= $#peaks; $i++) {
		if (defined $peaks[$i]) {
			push (@tempPeaks, $peaks[$i]);
		}
	}
	@peaks = @tempPeaks;
	undef @tempPeaks;

	##################################################################
	## Filtering step						##
	## A 3D-peaks may contain multiple peaks from one scan		##
	## In this case, choose only one with the largest intensity	##
	##################################################################	
	
	for (my $i = 0; $i < scalar(@peaks); $i++) {
		if (scalar(@{$peaks[$i]{'scanNumber'}}) != scalar(Unique(@{$peaks[$i]{'scanNumber'}}))) {
			my %tempHash;
			for (my $j = 0; $j < scalar(@{$peaks[$i]{'scanNumber'}}); $j++) {
				if (defined $tempHash{$peaks[$i]{'scanNumber'}[$j]}) {
					my $currentIntensity = $peaks[$i]{'intensity'}[$j];
					if ($currentIntensity > $tempHash{$peaks[$i]{'scanNumber'}[$j]}{'intensity'}) {
						$tempHash{$peaks[$i]{'scanNumber'}[$j]}{'intensity'} = $currentIntensity;
						$tempHash{$peaks[$i]{'scanNumber'}[$j]}{'index'} = $j;
					}
				} else {
					$tempHash{$peaks[$i]{'scanNumber'}[$j]}{'intensity'} = $peaks[$i]{'intensity'}[$j];
					$tempHash{$peaks[$i]{'scanNumber'}[$j]}{'index'} = $j;
				}
			}
			my @uniqueIndex;
			foreach my $key (sort {$a <=> $b} keys %tempHash) {
				push (@uniqueIndex, $tempHash{$key}{'index'});
			}
			@{$peaks[$i]{'centerMz'}} = @{$peaks[$i]{'centerMz'}}[@uniqueIndex];
			@{$peaks[$i]{'intensity'}} = @{$peaks[$i]{'intensity'}}[@uniqueIndex];
			@{$peaks[$i]{'maxMz'}} = @{$peaks[$i]{'maxMz'}}[@uniqueIndex];
			@{$peaks[$i]{'minMz'}} = @{$peaks[$i]{'minMz'}}[@uniqueIndex];
			@{$peaks[$i]{'scanIndex'}} = @{$peaks[$i]{'scanIndex'}}[@uniqueIndex];
			@{$peaks[$i]{'scanNumber'}} = @{$peaks[$i]{'scanNumber'}}[@uniqueIndex];
			@{$peaks[$i]{'scanRT'}} = @{$peaks[$i]{'scanRT'}}[@uniqueIndex];
		}
	}

	##########################################################
	## Smoothing of the intensity profile of each 3D-peak	##
	## Split a 3D-peak into multiple ones, if necessary 	##
	##########################################################

	$peaks3DCount = $#peaks;
	for (my $i = 0; $i <= $peaks3DCount; $i++) {
		## 1. Peak smoothing
		my @intensity = @{$peaks[$i]{'intensity'}};
		my ($smoothIntensity) = PeakSmoothing(\@intensity);	
		@{$peaks[$i]{'intensity'}} = @{$smoothIntensity};
		
		## 2. Peak split when local minimum is found
		my @splitCandidates;
		my @splitResults;
		push (@splitCandidates, $peaks[$i]);
		while ($#splitCandidates > -1) {
			my %maybeSplit = %{$splitCandidates[0]};
			shift (@splitCandidates);
			my ($needToSplit, $split1, $split2) = SplitIntoPair($$params{'valleyFactor'}, %maybeSplit);
			if ($needToSplit) {
				unshift (@splitCandidates, $split1, $split2);
				%maybeSplit = ();
			} else {
				push (@splitResults, \%maybeSplit);
			}
		}
		## Post-processing of peak split
		if ($#splitResults > 0) {
			for (my $j = 0; $j <= $#splitResults; $j++) {
				if ($j == 0) {
					$peaks[$i] = $splitResults[$j];
				} else {
					$peaks3DCount++;
					$peaks[$peaks3DCount] = $splitResults[$j];
				}
			}
		}
	}

	###################################################################################################
	## Calculation of "MEDIAN" of all smoothed intensities of 3D-peaks 				 ##
	## It will be used to update %msHash by excluding some low-intensity peaks even though they 	 ##
	## contributed to the generateion of 3D-peak 							 ##
	################################################################################################### 

	my @tempIntensity;
	foreach (@peaks) {
		push (@tempIntensity, @{$$_{'intensity'}});
	}
	my $medianIntensity3D = median(@tempIntensity);
	my $noisePercentage3D = $$params{'noisePercentage3D'};
	my $intensityThreshold3D = ($noisePercentage3D * $medianIntensity3D) / 100; 

	print "\n  Postprocessing of 3D-peaks is finished\n";

	##################################################
	## Estimation of 3D-peak masses and errors	##
	##################################################

	my $tot2DpeaksFor3Dpeaks = 0;
	my $removed2Dpeaks = 0;
#	my @tempRemoved;
	for (my $i = 0; $i <= $peaks3DCount; $i++) {
		my $nPoints = @{$peaks[$i]{'centerMz'}};
		my $numerator = 0;
		my $denominator = 0;
		for (my $j = 0; $j < $nPoints; $j++) {
			$numerator += $peaks[$i]{'centerMz'}[$j] * $peaks[$i]{'intensity'}[$j];
			$denominator += $peaks[$i]{'intensity'}[$j];
		}
		my $mass = $numerator / $denominator;	## Mass estimate calculated by the weighted mean of smoothed 3D-peak intensities
		my $minMass = min(@{$peaks[$i]{'centerMz'}});
		my $maxMass = max(@{$peaks[$i]{'centerMz'}});
		my $mErr = ($maxMass - $minMass) / $mass * 1e+6;
		$peaks[$i]{'3Dmass'} = $mass;
		#$peaks[$i]{'3DmassErr'} = $mErr;
		
		## Update ms hash #########################################################
		for (my $j = 0; $j < scalar(@{$peaks[$i]{'scanIndex'}}); $j++) {
			$tot2DpeaksFor3Dpeaks++;
			my $msScanIndex = $peaks[$i]{'scanIndex'}[$j];
			my $msSmoothedIntensity = $peaks[$i]{'intensity'}[$j];
			if ($msSmoothedIntensity >= $intensityThreshold3D) 
			{
				
				push (@{$ms{$msScanIndex}{'3Dmass'}}, $mass);
				push (@{$ms{$msScanIndex}{'3Dintensity'}}, $msSmoothedIntensity);
				#push (@{$ms{$msScanIndex}{'3DmassErrPpm'}}, $mErr);
				#push (@{$ms{$msScanIndex}{'3Dmass'}}, $mass);
				
			} else {
				$removed2Dpeaks++;
				my $entity = $i."_".$peaks[$i]{'scanIndex'}[$j];
#				push (@tempRemoved, $entity);
			}
		}
		###########################################################################
		print "  Updating MS peaks based on the information of 3D-peaks ($i out of $peaks3DCount)\r";
	}
	print "\n";
	print "  Out of $tot2DpeaksFor3Dpeaks (2D) peaks used to form 3D-peaks, $removed2Dpeaks have been EXCLUDED from the updated MS peaks due to their small intensities\n";

	##############################################################
	## Export the updated ms hash information to an output file	##
	##############################################################
	print "  Generating dta files\n";
	my @scanNumberArray;
	my %used_3Dpeaks;
	foreach my $index (keys %ms) {
		next if (!defined $ms{$index}{'3Dmass'});
		next if(defined($used_3Dpeaks{$ms{$index}{'3Dmass'}}));
		
		my $scanNumber = $ms{$index}{'scanNumber'};
		push (@scanNumberArray,$scanNumber);
		
		my $iDtaFile = "$dtaPath/$scanNumber" . ".MS1";

		open (DTAFILE, ">", $iDtaFile);
		my $sizeArray = scalar(@{$ms{$index}{'3Dmass'}});
		## Sort 3Dmass and 3Dintensity
		my @sortedIndex = sort {$ms{$index}{'3Dmass'}[$a]<=>$ms{$index}{'3Dmass'}[$b]} 0..$#{$ms{$index}{'3Dmass'}};
		@{$ms{$index}{'3Dmass'}} = @{$ms{$index}{'3Dmass'}}[@sortedIndex];
		@{$ms{$index}{'3Dintensity'}} = @{$ms{$index}{'3Dintensity'}}[@sortedIndex];
		print DTAFILE "1\t2\n";
		#@{$ms{$index}{'3DmassErrPpm'}} = @{$ms{$index}{'3DmassErrPpm'}}[@sortedIndex];
		for (my $i = 0; $i < $sizeArray; $i++) {
			my $mz = $ms{$index}{'3Dmass'}[$i];
			my $int = $ms{$index}{'3Dintensity'}[$i];
			#print $mz,"\t",$int,"\t",$index,"\t",$i,"\n";
			print DTAFILE "$mz\t$int\n";			
		}
	}
	undef %ms;

	print "  Generating 2D hash for each 3D mass\n";
	##################################################
	## Generation of %msToPeaks hash		##
	## Added on 10/15/2014 per Xusheng's request	##
	##################################################
	
	for (my $i = 0; $i < scalar(@peaks); $i++) 
	{
		print "\r  store peak $i   ";
		my %info_2D;
		@{$info_2D{'2Dmass'}} = @{$peaks[$i]{'centerMz'}};
		@{$info_2D{'2Dintensity'}} = @{$peaks[$i]{'intensity'}};
		@{$info_2D{'2DscanNumber'}} = @{$peaks[$i]{'scanNumber'}};
		@{$info_2D{'2Drt'}} = @{$peaks[$i]{'scanRT'}};
		my $file = $peaks[$i]{'3Dmass'}; # . "_" . shift(@{$peaks[$i]{'scanNumber'}});
		store(\%info_2D, "$dtaPath/$file");
	}

#	store(\@peaks, "$dtaPath/peaks_hash");
#	store (\%ms, $$params{'outputFile'});
	print "\n  Finished 3D peaks generation\n\n";
	##########################################
	## Modified on 10/15/2014 by JCho	##
	##########################################

	@peaks=();

	#return (\%ms,\@peaks, \%msToPeaks);
}

##################
## Subroutines	##
##################

sub Unique {
	my %seen;
	my @uniq =  grep { !$seen{$_}++ } @_;
	return (@uniq);
}

sub SplitIntoPair {
	my ($valleyFactor, %peak) = @_;
	my $needToSplit = 0;
	my %split1;
	my %split2;
	my @intensity = @{$peak{'intensity'}};
	my @minPos = CalcLocalMinPositions(@intensity);	
	foreach my $pos (@minPos) {
		my $leftMax = GetLeftMax($pos, @intensity);
		my $rightMax = GetRightMax($pos, @intensity);
		my $smallMax = min($leftMax, $rightMax);
		if ($smallMax / $intensity[$pos] > $valleyFactor) {
			foreach my $key (keys %peak) {
				@{$split1{$key}} = @{$peak{$key}}[0..$pos];
				@{$split2{$key}} = @{$peak{$key}}[$pos + 1..$#intensity];
			}
			$needToSplit = 1;
			return ($needToSplit, \%split1, \%split2);
		}
	}
	return ($needToSplit, \%split1, \%split2);
} 

sub GetLeftMax {
	my ($pos, @y) = @_;
	my $max = -1;
	for (my $i = 0; $i < $pos; $i++) {
		if ($y[$i] > $max) {
			$max = $y[$i];
		}
	}
	return $max;
}

sub GetRightMax {
	my ($pos, @y) = @_;
	my $length = @y;
	my $max = -1;
	for (my $i = $pos + 1; $i < $length; $i++) {
		if ($y[$i] > $max) {
			$max = $y[$i];
		}
	}
	return $max;	
}

sub CalcLocalMinPositions {
	my @y = @_; ## smoothed intensity
	my @minPos;
	my $length = @y;
	for (my $i = 2; $i < $length - 2; $i++) {
		my $b2 = $y[$i - 2];
		my $b1 = $y[$i - 1];
		my $x = $y[$i];
		my $a1 = $y[$i + 1];
		my $a2 = $y[$i + 2];
		if (IsMin($b2, $b1, $x, $a1, $a2)) {
			push (@minPos, $i);
		}
	}
	my @minY = @y[@minPos];
	my @yOrderedIndex = sort {$minY[$a] <=> $minY[$b]} 0..$#minY;
	@minPos = @minPos[@yOrderedIndex];
	return (@minPos);
}

sub IsMin {
	my ($b2, $b1, $x, $a1, $a2) = @_;
	my $IsMin = 0;
	if ($x < $b1 && $x < $a1) {
		$IsMin = 1;
	} elsif ($x == $b1 && $x < $b2 && $x < $a1) {
		$IsMin = 1;
	} elsif ($x == $a1 && $x < $a2 && $x < $b1) {
		$IsMin = 1;
	} elsif ($x < $b2 && $x == $b1 && $x == $a1 && $x < $a2) {
		$IsMin = 1;
	} else {
		$IsMin = 0;
	}
	return $IsMin;	
}

sub PeakSmoothing {
	my ($origIntensity) = @_;
	my @smoothIntensity = SmoothMedian(1, @{$origIntensity});
	@smoothIntensity = SmoothMean(2, @smoothIntensity);
	return (\@smoothIntensity);
}

sub SmoothMean {
	my ($width, @m) = @_;
	my $length = @m;
	my @result;
	for (my $i = 0; $i <$length; $i++) {
		my $min = max(0, $i - $width);
		my $max = min($length - 1, $i + $width);
		$result[$i] = mean(@m[$min..$max]);
	}
	return (@result);
}

sub SmoothMedian {
	my ($width, @m) = @_;
	my $length = @m;
	my @result;
	for (my $i = 0; $i <$length; $i++) {
		my $min = max(0, $i - $width);
		my $max = min($length - 1, $i + $width);
		$result[$i] = median(@m[$min..$max]);
	}
	return (@result);	
}

sub ExtractMsCentroidList {
	my ($hash, $index) = @_;
	my %reduced;
	my $subInd = 0;
	foreach my $key (keys %{$hash}) {
		if ($key eq "scanIndex" || $key eq "scanNumber" || $key eq "RT") {
			$reduced{$key} = $$hash{$key};
		} elsif ($key eq "peakCenterMassHash") {
			next;
		} elsif ($key eq "peakCenterMassIndexHash") {
			next;
		} elsif ($key eq "peakCenterMass") {
			my @sub = @{$$hash{$key}}[@$index];
			$reduced{$key} = \@sub;
			$reduced{'peakCenterMassHash'} = {};
			$reduced{'peakCenterMassIndexHash'} = {};
			foreach my $subMass (@sub) {
				push (@{$reduced{'peakCenterMassHash'}{int($subMass)}}, $subMass);
				$reduced{'peakCenterMassIndexHash'}{$subMass} = $subInd;
				$subInd++;
			}
		} else {
			my @sub = @{$$hash{$key}}[@$index];
			$reduced{$key} = \@sub;
		}
	}
	return (%reduced);
}

sub GetClosestValue {
	my ($hash, $mass, $tol) = @_;
	my $lL = $mass - ($mass * $tol / 1e6);
	my $uL = $mass + ($mass * $tol / 1e6);
	my @searchKeys = (int($lL)..int($uL));
	
	my @returnIndex;
	my $intensity = 0;
	my $match = 0;
	
	foreach my $key (@searchKeys) {
		next if (!defined $$hash{'peakCenterMassHash'}{$key});
		for (my $i = 0; $i < scalar(@{$$hash{'peakCenterMassHash'}{$key}}); $i++) {
			my $value = $$hash{'peakCenterMassHash'}{$key}[$i];			
			if ($value >= $lL && $value <= $uL) {
				## returnIndex can be an array (there might be multiple indices in one scan)
				push (@returnIndex, $$hash{'peakCenterMassIndexHash'}{$value});
				$match = 1;
			}
		}
	}
	return ($match, \@returnIndex);
}

sub DetectPeaks {
	## Input is a single MS spectrum (mz values and corresponding intensities)
	my ($numCentroidPoints, $thresholdScale, $threshold, $msHash) = @_;
	my @peakIntensity;
	my @peakCenterMass;
	my @peakMinMass;
	my @peakMaxMass;
	my @peakIndex;
	my %peakCenterMassHash;
	my %peakCenterMassIndexHash;
	my $peakCount = 0;
	my $intensityThreshold = 0;
	
	if ($thresholdScale == 1) {	## Absolute intensity threshold for 2D-peaks
		$intensityThreshold = $threshold;
	} elsif ($thresholdScale == 2) {	## Relative intensity threshold for 2D-peaks
		if ($threshold > 0) {
			my @nonzeroIntensity;
			foreach (@{$$msHash{'intensity'}}) {
				if ($_ > 0) {
					push (@nonzeroIntensity, $_);
				}
			}
		$intensityThreshold = ($threshold * median(@nonzeroIntensity)) / 100;
		}
	}
	my $counts = @{$$msHash{'intensity'}};
	
	for (my $i = 2; $i < $counts - 2; $i++) {
		if ($$msHash{'intensity'}[$i] > 0) {
			# Consider 2 points before and after the point of interest (x), i.e. 5 point-window
			my $b2 = $$msHash{'intensity'}[$i - 2];
			my $b1 = $$msHash{'intensity'}[$i - 1];
			my $x = $$msHash{'intensity'}[$i];
			my $a1 = $$msHash{'intensity'}[$i + 1];
			my $a2 = $$msHash{'intensity'}[$i + 2];
			if ($x >= $intensityThreshold) {
				if (IsMax($b2, $b1, $x, $a1, $a2)) {
					# If x is the local maximum within a 5-point window, lower and upper bounds for a peak will be explored
					# Refer Figure 1a and b in the paper (Cox and Mann, Nature Biotechnology, 2008, 26: 1367-72)
					my $minInd = CalcMinPeakIndex($i, $$msHash{'intensity'});	# Search for a lower bound of the peak
					my $maxInd = CalcMaxPeakIndex($i, $$msHash{'intensity'});	# Search for a upper bound of the peak
					if (($maxInd - $minInd) > 2) {
						my ($calculatedIntensity, $calculatedMass) = CalcCenterMass($numCentroidPoints, $minInd, $i, $maxInd, $$msHash{'mass'}, $$msHash{'intensity'});
						if ($calculatedMass =~ /^?\d+\.?\d*$/) {
							$peakCenterMass[$peakCount] = $calculatedMass;
							$peakIntensity[$peakCount] = $calculatedIntensity;
							push (@{$peakCenterMassHash{int($calculatedMass)}}, $calculatedMass);
							$peakCenterMassIndexHash{$calculatedMass} = $peakCount;
							$peakIndex[$peakCount] = [$minInd..$maxInd];
							$peakCount++;
						}
					}
				}
			}
		}
	}
	my %msCentroidList = ('peakCenterMass' => \@peakCenterMass, 'peakCenterMassHash' => \%peakCenterMassHash, 
		'peakCenterMassIndexHash' => \%peakCenterMassIndexHash, 'peakIntensity' => \@peakIntensity, 'peakIndex' => \@peakIndex);
	return(\%msCentroidList);
}

sub IsMax {	## IsMax determines whether the middle point reaches a local maximum
	my ($b2, $b1, $x, $a1, $a2) = @_;
	my $IsMax = 0;
	if ($x > $b1 && $x > $a1) {
		$IsMax = 1;		
	} elsif ($x == $b1 && $x > $b2 && $x > $a1) {
		$IsMax = 1;
	} else {
		$IsMax = 0;
	}
	return $IsMax;
}

sub CalcMinPeakIndex {	## CalcMinPeakIndex determines the lower bound of a peak
	my ($ind, $intensity) = @_;
	while ($ind > 0 &&  $$intensity[$ind] != 0 && $$intensity[$ind - 1] <= $$intensity[$ind] ) {
		$ind--;
	}
	return $ind + 1;
}

sub CalcMaxPeakIndex {	## CalcMaxPeakIndex determines the upper bound of a peak
	my ($ind, $intensity) = @_;
	my $count = @$intensity;
	while ($ind < $count && $$intensity[$ind] != 0 && $$intensity[$ind + 1] <= $$intensity[$ind] ) { 
		$ind++;
	}
	return $ind - 1;
}

sub CalcCenterMass {	## CalcCenterMass estimates the mass (m/z) value at the center of a peak through least square fitting
	my ($nPoints, $minInd, $centerInd, $maxInd, $mass, $intensity) = @_;
	my $peakIntensity;
	my $peakCenterMass;
	for (my $j = $minInd; $j <= $maxInd; $j++) {
		$peakIntensity += $$intensity[$j];	# peakIntesntiy = Sum of intensities (or maximum intensity)
	}
	## There's a maximum point, but other points are zeros
	if ($minInd == $maxInd) {
		$peakCenterMass = $$intensity[$maxInd];
		return ($peakIntensity, $peakCenterMass);
	}
	## Right angled triangle-shaped peak
	if ($minInd == $centerInd) {
		$peakCenterMass = Estimate2($$mass[$centerInd], $$mass[$centerInd + 1], $$intensity[$centerInd], $$intensity[$centerInd + 1]);
		return ($peakIntensity, $peakCenterMass);
	}
	## Left angled triangle-shaped peak
	if ($maxInd == $centerInd) {
		$peakCenterMass = Estimate2($$mass[$centerInd - 1], $$mass[$centerInd], $$intensity[$centerInd - 1], $$intensity[$centerInd]);
		return ($peakIntensity, $peakCenterMass);
	}
	## Typical bell(triangle)-shaped peak
	if ($nPoints <= 3) {
		$peakCenterMass = Estimate3($$mass[$centerInd - 1], $$mass[$centerInd], $$mass[$centerInd + 1], $$intensity[$centerInd - 1], $$intensity[$centerInd], $$intensity[$centerInd + 1]);
		return ($peakIntensity, $peakCenterMass);
	}
	
	## Typical bell(triangle)-shaped peak, use more than three data points	
	if ($nPoints > 3) {
		my $nleft;
		my $nright;
		my $d;
			
		## Select data points to be used in the least square fitting
		if ($nPoints % 2 == 1) {	## Simple case when $npoints is an odd number
			$d = int($nPoints / 2);
			$nleft = max($centerInd - $d, $minInd);
			$nright = min($centerInd + $d, $maxInd);
		} else {	## If $npoints is an even number, it would be a bit complicated
			$d = $nPoints / 2 - 1;
			$nleft = max($centerInd - $d, $minInd);
			$nright = min($centerInd + $d, $maxInd);
			if ($nleft != $minInd && $nright != $maxInd) {
				if ($$intensity[$nleft - 1] > $$intensity[$nright + 1]) {
					$nleft--;
				} else {
					$nright++;
				}
			} elsif ($nleft != $minInd) {
				$nleft--;
			} elsif ($nright != $maxInd) {
				$nright++;
			}
		}
		
		my @x;
		my @y;
		for (my $i = 0; $i < $nright - $nleft + 1; $i++){
			$x[$i] = $$mass[$nleft + $i];
			$y[$i] = $$intensity[$nleft + $i];
		}
		$peakCenterMass = EstimateN(\@x, \@y, $centerInd - $nleft);
		return ($peakIntensity, $peakCenterMass);	
	}	
}

sub Estimate2 {	## Estimate2 gives a peakCenterMass by calculating an intensity-weighted average of mass
	my ($m1, $m2, $i1, $i2) = @_;
	my $peakCenterMass = ($m1 * $i1 + $m2 * $i2) / ($i1 + $i2);
	return $peakCenterMass;
}

sub Estimate3 {	## Estimate3 gives a peakCenterMass using the closed form of least-square fitting to a Gaussian-shaped peak
	my ($m1, $m2, $m3, $i1, $i2, $i3) = @_;
	my $l1 = log($i1);
	my $l2 = log($i2);
	my $l3 = log($i3);
	my $peakCenterMass = 0.5 * (($l1 - $l2) * ($m3 * $m3 - $m1 * $m1) - ($l1 - $l3) * ($m2 * $m2 - $m1 * $m1)) / (($l1 - $l2) * ($m3 - $m1) - ($l1 - $l3) * ($m2 - $m1));
	return $peakCenterMass;
}

sub EstimateN {	## EstimateN gives a peakCenterMass using least-square fitting to a Gaussian-shaped peak
	my ($x, $y, $ic) = @_;
	my @x = @{$x};
	my @y = @{$y};
	my $dm = $x[$ic];	## Intensity of centerInd
	my $regInput;
	for (my $i = 0; $i < @x; $i++) {
		$regInput -> [$i] = [log($y[$i]), ($x[$i]-$dm), ($x[$i]-$dm)**2];	## Mean-centering for x (singularity prevention???)		
	}
	my ($coeff, $Rsq) = linear_regression($regInput);	## Fitting to a Gaussian peak shape -> equivalent to the fitting of x and log(y) to a quadratic curve
	my $peakCenterMass = $dm - ($coeff -> [1]) / ($coeff -> [2]) * 0.5;
	return $peakCenterMass;
}

sub getMsHash {	
	my ($mzxml, $startScan, $endScan) = @_;	# Name of input .mzXML file
	my %ms;
	open (XML, "<", $mzxml) or die "Cannot open $mzxml file\n";
	while (<XML>) {
		if ($_ =~ /centroided="1"/) {
			print "  Make sure that .mzXML file has been generated in PROFILE mode\n";
			exit;
		}
		last if ($_ =~ /\<\/dataProcessing/);
	}
	my $indexOffset = getIndexOffset(*XML);
	my ($indexArray, $lastScan) = getIndexArray(*XML, $indexOffset);
	my $scanNumber = 0;
	my $scanIndex = 0;
	foreach (@$indexArray) {		
		$scanNumber++;
		if ($scanNumber < $startScan) {
			next;
		} elsif ($scanNumber > $endScan) {
			last;
		}
		print "  Gathering information of scan #$scanNumber out of $lastScan\r";
		my $index = $_;
		my $mslevel = getMSLevel(*XML, $index);
		next if ($mslevel != 1);
		my $rt = getRT(*XML, $index);
		my @peakArray;
		getPeak(*XML, \@peakArray, $index);	
		my $sizeArray = scalar(@peakArray);
		my @mz;
		my @int;
		for (my $i = 0; $i < $sizeArray; $i++) {
			$i++;
			push (@mz, $peakArray[$i - 1]);
			push (@int, $peakArray[$i]);
		}
		$ms{$scanIndex}{'RT'} = $rt;
		$ms{$scanIndex}{'scanNumber'} = $scanNumber;
		$ms{$scanIndex}{'mass'} = \@mz;
		$ms{$scanIndex}{'intensity'} = \@int;
		$scanIndex++;
	}
	close XML;
	return %ms;
}

sub getIndexOffset{
	(*XML) = @_;
	seek (XML, -120, 2);
	my ($indexOffset);
	while (<XML>){
	  next if (!/<indexOffset>/);
	  chomp;
	  $indexOffset = $_;
	  $indexOffset =~ s/\s+<indexOffset>(\d+)<\/indexOffset>.*/$1/o;
	  last;
	}
	return $indexOffset;
}

sub getIndexArray{	
	(*XML, my $indexOffset) = @_;
	my @indexArray;
	my $lastScan = 0;
	seek (XML, $indexOffset, 0);
	while (<XML>){
	  next if (/^\s+<index/);
	  last if (/^\s+<\/index>.*/);
	  chomp;
		next if (/scan/);
	  $lastScan++;
	  my $index = $_;
	  $index =~ s/[\s\t]+\<offset id="\d+"[\s]*>(\d+)<\/offset>.*/$1/o;
	  push (@indexArray, $index);
	}
	return (\@indexArray, $lastScan);
}

sub getRT{
	(*XML, my $scanIndex) = @_;
	my $rt;
	seek (XML, $scanIndex, 0);
	while (<XML>){
	  next if (!/retentionTime=/);
	  chomp;
	  $rt = $_;
	  $rt =~ s/.*retentionTime="PT([\d+\.]+)S".*/$1/o;
	  last;
	}	
	return $rt;
}

sub getMSLevel{
	(*XML, my $scan_index) = @_;
	my $msLevel;	
	seek (XML, $scan_index, 0);
	while (<XML>){
	  next if (!/msLevel=/);
	  chomp;
	  $msLevel = $_;
	  $msLevel =~ s/\s+msLevel="(\d)".*/$1/o;
	  last;
	}	
	return $msLevel;
}

sub getPeak{
	(*XML, my $peakArray, my $scanIndex) = @_;
	my ($peak, $peakLine);
	seek (XML, $scanIndex, 0);
	while(<XML>){
		if (/<peaks\sprecision="32"/ || /pairOrder="m\/z-int"[\s]*>/){
			chomp;
			next if ($_ =~ /<peaks precision="32"\Z/);
			$peakLine = $_;
			if (/<peaks precision="32">/){
				$peakLine =~ s/\s+<peaks precision="32"[\s\w\W\d\=\"]+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			} else {
				$peakLine =~ s/\s+pairOrder="m\/z-int"[\s]*>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			}
			last;
		} elsif (/compressedLen/){
			chomp;
			$peakLine = $_;
			$peakLine =~ s/\s+compressedLen="[\d]"\s+>([A-Za-z0-9\/\+\=]+)<\/peaks>.*/$1/o;
			last;
		} else {
			next;
		}
	}
	#Base64 Decode
	$peak = decode_base64($peakLine);
	my @hostOrder32 = unpack("N*", $peak);
	for (@hostOrder32){
		my $float = unpack("f", pack("I", $_));
		push (@$peakArray, $float);
	}
}

sub IndMax{
	my (@array) = @_;
	my $valMax = max(@array);
	my @indMax;
	for (my $i = 0; $i < scalar(@array); $i++) {
		if ($array[$i] == $valMax) {
			push @indMax, $i;
		}
	}
	return (@indMax);
}

1;
