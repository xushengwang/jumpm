#!/usr/bin/perl

######### Job #################################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Job;


use strict;
use warnings;
use File::Basename;
use Storable;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.03;


@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_path => undef,
        _sim_path  => undef,
    };
    bless $self, $class;
    return $self;
}

sub set_dta_path
{
	my ($self,$dta_path)=@_;
	$self->{_dta_path}=$dta_path;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{_dta_path};
}

sub set_pip
{
	my ($self,$pip)=@_;
	$self->{_pip}=$pip;	
}

sub get_pip
{
	my $self=shift;
	return $self->{_pip};
}

sub set_library_path
{
	my ($self,$lib)=@_;
	$self->{_lib}=$lib;	
}

sub get_library_path
{
	my ($self)=@_;
	return $self->{_lib};	
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


sub create_frag_script
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	my $parameter = $self->get_parameter();
	
	my $neutralLossFile = $parameter->{'neutralLossFile'};
	my $bondEnergyFile =  $parameter->{'bondEnergyFile'};

	
	open(RUNSHELL,">$dir/frag_shell.pl");
=head
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::MS2_scoring;
use Spiders::Dta;

my ($help,$dtafile,$smile,$mass,$depth);
GetOptions('-help|h'=>\\\$help,
		'-dtafile=s'=>\\\$dtafile,
		'-smile=s'=>\\\$smile,
		'-mass=s'=>\\\$mass,
		'-depth=s'=>\\\$depth,	
		);
		
		SmileFragmentation(\$smile, \$mass, \$depth, \$neutralLossFile, \$bondEnergyFile);
		my \@frag_mz_values=();
		for(my \$i=2;\$i<\$\#frag_return;\$i++) 
		{
			my \@frag_data = split(/\\s+/,\$frag_return[\$i]);
			next if(\$frag_data[-1]=~/[CHNOPS]/);
			push (\@frag_mz_values,\$frag_data[-1]);
		}
		my \$dta=new Spiders::Dta;  
		\$dta->set_dta_file(\$dtafile) ;
		\$dta->process_dtafile(\$dtafile);

		my \%mz_hash = %{\$dta->{'_mz_hash'}};
	
		my \$scoring = new Spiders::MS2_scoring;
		\$scoring->compare_theoritical_experiment(\\\%mz_hash,\\\@frag_mz_values);
		\$scoring->write_result();
EOF
=cut
	close(RUNSHELL);			
}

sub create_pscript
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(PAIRSHELL,">$dir/pair_shell.pl");
print PAIRSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	
my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

my \@dtafiles = \@ARGV;


my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;

    my \$dtafile = abs_path(\$_);


	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(MRESULT,">\${dtafile}.pairs");	
    my \$dta=new Spiders::Dta;  
    \$dta->set_dta_file(\$dtafile) ;
    \$dta->process_dtafile(\$dtafile);

    my \%mz_hash = %{\$dta->{'_mz_int_hash'}};

	
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_NC_ratio(0);
	\$pair->set_CC_ratio(0);	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
            next if (\!defined(\$mz_hash{\$mz}));
			\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
            \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

            \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	my \$dirname  = dirname(\$dtafile);
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);
    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname);
	
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>2500);
		next if(\$#missile_formula<0);
		\$mono_mass = \$mono_mass - 1.007277;
		my \$query = new Spiders::DatabaseQuery();		
		my (\$return_all,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_pairing}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "hmdb");
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "pubchem");
		#my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all", "kegg");


		foreach my \$missile_formula (\@missile_formula)
		{
			my \@N_C = split(/\\,/,\$missile_formula);
			my (\$C_num_mz,\$C_mass_defect,\$C_relative_int) = split(/\\|/,\$N_C[0]);
			my (\$N_num_mz,\$N_mass_defect,\$N_relative_int) = split(/\\|/,\$N_C[1]);
				
			my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
			my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	
					
			foreach my \$theoretical_formula (\@\$return_norm) 
			{
				chomp \$theoretical_formula;
				my \$carbon_match = 0;
				my \$nitrogen_match = 0;

				\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
				\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
				if(\$carbon_match == 1 and \$nitrogen_match ==1)
				{				
					print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$C_mass_defect,"\\t", \$C_relative_int,"\\t", \$N_mass_defect,"\\t",\$N_relative_int,"\\n";					
				}
			}
		}					
	}
}	
EOF
	close(PAIRSHELL);	
}
	
sub create_mscript_unlabel
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();
	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;

GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,		
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	
my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

my \@dtafiles = \@ARGV;

my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;

    my \$dtafile = abs_path(\$_);
	
	print "\\nProcessing scan: \$dtafile\\n";		
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	
    print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tMass\\tMISSILE\\tFormula\\tIUPAC name\\tScore\\tStructure (SMILES)\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	\$mz_hash_ref = \$deisotope->MS1_deisotope();
	\%mz_hash = \%\$mz_hash_ref;
	
    print scalar(keys \%mz_hash)," features were detected\\n";	
	
    foreach my \$mono_mass (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        print "performing mass formula database search for mass: \$mono_mass\\n";
		
		\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		my \$query = new Spiders::DatabaseQuery();
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);
		
		# my \$return = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm = \$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		\$mono_mass = \$mono_mass + 1.007277;
		
        print "found ",scalar(\@\$return)," formulas, including ",scalar(\@\$return_norm)," targets and ",scalar(\@\$return_decoy)," decoys\\n";	
		
		foreach my \$theoretical_formula (\@\$return) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\n";
					
			\$index++;					
		}
					
		foreach my \$theoretical_formula (\@\$return_norm) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t";						
			print MRESULT \$formula_db,"\t",\$mass_db,"\\tnorm\\n";
			my \$query = new Spiders::DatabaseQuery();
								
			\$query->setDatabase(\$params->{structure_database});
			\$query->setFormula(\$formula_db);
			\$query->setMass(\$mass_db);
			print "querying structure database\\n";
			\$str_db_return = \$query->QueryStructureDatabase();				
            print "formula ",\$formula_db, " found ",scalar(\@\$str_db_return)," structures\\n";

			if(\$\#\$str_db_return>0)
			{
				foreach (\@\$str_db_return)
				{
					my \$smiles="";
					my \$IUPAC = "";
					(\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
					print RESULT "norm\\t",\$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_db,"\\t",\$IUPAC,"\\t",\$smiles,"\\n";												
				}
			}
			else
			{
				print RESULT "norm\\t",\$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_db,"\\n";						
					
			}
					
			\$index++;					
		}


		foreach my \$theoretical_formula (\@\$return_decoy) 
		{				
			chomp \$theoretical_formula;
			my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
			print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t",\$formula_db,"\t",\$mass_db,"\\tdecoy\\n";
			\$index++;					
		}
	}				
}
			
EOF
	close(RUNSHELL);
}	
	
sub create_mscript
{
	my ($self)=@_;
	my $dir = $self->get_dta_path();
	my $lib = $self->get_library_path();

	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;
#!/usr/bin/perl  -I $lib
use Getopt::Long;
use Spiders::Deisotope;
use Spiders::Params;
use Spiders::Dta;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::DatabaseQuery;
use Spiders::Pairing;


GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-peaks=s'=>\\\$peaks,
		'-NC=f'=>\\\$NC,
		'-CC=f'=>\\\$CC,
		'-NC_std=f'=>\\\$NC_std,
		'-CC_std=f'=>\\\$CC_std,
		'-NC_defect_loc=f'=>\\\$NC_defect_loc,
		'-CC_defect_loc=f'=>\\\$CC_defect_loc,
		'-NC_defect_scale=f'=>\\\$NC_defect_scale,
		'-CC_defect_scale=f'=>\\\$CC_defect_scale,
		);
		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	
my \$database=\$params->{'database'};
my \$decoy_database=\$database . "_decoy";

my \@dtafiles = \@ARGV;

my \$ms_hash_mol;
foreach(\@dtafiles)
{
    my \$index=1;
    my \$dtafile = abs_path(\$_);
	print "\\nProcessing scan: \$dtafile\\n";	
	
    my \$scan = \$dtafile;
    \$scan =~s\/\(\.\*\)\\\.\(\\d\+\)\\\.dta\/\$2\/;	
    open(RESULT,">\${dtafile}.smout");
    open(MRESULT,">\${dtafile}.missile");	
	print MRESULT "Index","\\t","C12 mass","\\t","N15 mass","\\t","C13 mass","\\t","pair_score","\\t","C12_C13_int_diff","\\t", "C12_C13_int_sim","\\t", "C12_N15_int_diff","\\t", "C12_N15_int_sim","\\t", "C12_C13 score", "\\t", "C12_N15 score", "\\t", "Scan number", "\\t","C12 intensity","\\tType\\n";
		
    print RESULT "MS1 scan: \$dtafile\\n";
    print RESULT "Index\\tMass\\tMISSILE\\tFormula\\tIUPAC name\\tScore\\tStructure (SMILES)\\n";
	
	print "MS1 Deisotoping\\n";			
	my \$deisotope = new Spiders::Deisotope();  
	\$deisotope->set_parameter(\$params);	
                        
    \$deisotope->set_dta(\$dtafile);
	\$mz_hash_ref = \$deisotope->MS1_deisotope();
    #my \$dta=new Spiders::Dta;  
    #\$dta->set_dta_file(\$dtafile) ;
    #\$dta->process_dtafile(\$dtafile);

    #my \%mz_hash = %{\$dta->{'_mz_int_hash'}};
	\%mz_hash = \%\$mz_hash_ref;	
    print scalar(keys \%mz_hash)," features were detected\\n";		
	
	print "Peaks pairing\\n";			
    my \$cluster;
	my \$pair = new Spiders::Pairing();
	\$pair->set_parameter(\$params);
	\$pair->set_NC_ratio(\$NC);
	\$pair->set_CC_ratio(\$CC);
	\$pair->set_NC_std(\$NC_std);
	\$pair->set_CC_std(\$CC_std);
	\$pair->set_NC_defect_loc(\$NC_defect_loc);
	\$pair->set_CC_defect_loc(\$CC_defect_loc);
	\$pair->set_NC_defect_scale(\$NC_defect_scale);
	\$pair->set_CC_defect_scale(\$CC_defect_scale);
	
    foreach my \$mz (reverse sort {\$mz_hash{\$a}<=>\$mz_hash{\$b}} keys \%mz_hash)
    {
        next if (\!defined(\$mz_hash{\$mz}));
		\$cluster->{\$mz}->{0}->{'mz'} = \$mz;
        \$cluster->{\$mz}->{0}->{'int'} = \$mz_hash{\$mz};

        \$pair->clustering(\$cluster,\\\%mz_hash,\$mz);
    }
	
	
############### only for summary purpose ########
    my \$dta=new Spiders::Dta;                  #
	\$dtafile .= ".iso";                         #
	\$dta->set_prec_mz("1");                     #
    \$dta->set_charge("0");	                     #
    \$dta->set_dta_file(\$dtafile) ;	        #
	\$dta->set_mz_int_hash(\\\%mz_hash);	
	\$dta->write_dta_file();                    #
#################################################
	
	
	
	my \$dirname  = dirname(\$dtafile);
	
    my \$cluster_mono = \$pair->select_mono_peaks(\$cluster);

    my (\$formula_comb) = \$pair->pairing(\$cluster_mono,\$dirname);
	
	\$index=1;	
	
	print scalar(keys \%\$formula_comb)," MISSILES were detected\\n";
    foreach my \$mono_mass (sort {\$a<=>\$b} keys \%\$formula_comb)
    {
        print "performing mass formula database search for mass: \$mono_mass\\n";
		
        my \@missile_formula = \@{\$formula_comb->{\$mono_mass}->{'number'}};
        next if(\$mono_mass>2500);
		next if(\$#missile_formula<0);
###### For decoy testing		
		my \$shift = 0;
		\$mono_mass = \$mono_mass - 1.007277 + \$shift;
		my \$query = new Spiders::DatabaseQuery();
		my (\$return,\$return_norm,\$return_decoy) = \$query->QueryMassWithDatabaseFilter(\$mono_mass, \$params->{formula_mass_tolerance_searching}, \$params->{mass_formula_database}, "yes",\$database,\$decoy_database);		
		# my \$return =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{mass_formula_database}, "yes");
		# my \$return_decoy =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{decoy_str_mass_formula_database}, "all");
		# my \$return_norm =\$query->QueryMassWithDatabase(\$mono_mass, \$params->{formula_mass_tolerance}, \$params->{norm_str_mass_formula_database}, "all");
		\$mono_mass = \$mono_mass + 1.007277;
##### if the case of "mass only" 
		if(\$params->{'labeled_ID_method'} == 1)
		{
			foreach my \$theoretical_formula (\@\$return) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t";						
				print MRESULT \$formula_db,"\t",\$mass_db,"\\n";
						
				\$index++;					
			}
									
			foreach my \$theoretical_formula (\@\$return_norm) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t";						
				print MRESULT \$formula_db,"\t",\$mass_db,"\\tnorm\\n";

									
				\$query->setDatabase(\$params->{structure_database});
				\$query->setFormula(\$formula_db);
				\$query->setMass(\$mass_db);
				\$str_db_return = \$query->QueryStructureDatabase();				
	
				if(\$\#\$str_db_return>0)
				{
					foreach (\@\$str_db_return)
					{
						my (\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
						print RESULT "Norm\\t",\$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_db;						
						print RESULT "\\t",\$IUPAC,"\\t",\$smiles,"\\n";								
					}
				}
				else
				{
					print RESULT "Norm\\t",\$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_db,"\\n";						
						
				}
						
				\$index++;					
			}


			foreach my \$theoretical_formula (\@\$return_decoy) 
			{				
				chomp \$theoretical_formula;
				my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$scan,"\\t";						
				print MRESULT \$formula_db,"\t",\$mass_db,"\\tdecoy\\n";
						
				\$index++;					
			}			
			
			
		}
		else
		{		
			foreach my \$missile_formula (\@missile_formula)
			{
				my \@N_C = split(/\\,/,\$missile_formula);
				my (\$C_num_mz,\$c_diff_score,\$c_int_score, \$c_similarity) = split(/\\|/,\$N_C[0]);
				my (\$N_num_mz,\$n_diff_score,\$n_int_score, \$n_similarity) = split(/\\|/,\$N_C[1]);
				my \$pair_score = \$N_C[2];
				my (\$C_num,\$C_mz) = split(/\\:/,\$C_num_mz);
				my (\$N_num,\$N_mz) = split(/\\:/,\$N_num_mz);	
				print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity, "\\t", \$scan, "\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\tAll\\n";
				
				
			
				foreach my \$theoretical_formula (\@\$return) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = 1;
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan, "\\t", \$formula_comb->{\$mono_mass}->{'intensity'}, "\\t","pass", "\\t", \$formula_db, "\\n";			
					}

				}
			


####### generating structure and scoring #######################
			
				foreach my \$theoretical_formula (\@\$return_norm) 
				{
				
				
					chomp \$theoretical_formula;				
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = 1;
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\\\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{
					
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);

						my \$query = new Spiders::DatabaseQuery();
										
						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$formula_db);
						\$query->setMass(\$mass_db);
						\$str_db_return = \$query->QueryStructureDatabase();				

						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t","Norm","\\t",\$formula_db,"\\n";					
						
						if(\$\#\$str_db_return>0)
						{
							foreach (\@\$str_db_return)
							{
								my (\$InChI,\$smiles,\$IUPAC)=split(\/\\t\/,\$_);
								print RESULT "NORM\\t",\$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t",\$formula_db;						
								print RESULT "\\t",\$IUPAC,"\\t",\$smiles,"\\t",\$InChI,"\\n";	
								my \$smiles_orig = \$smiles;
								\$index++;
								\$smiles=\~s/\\\(/\\\\\(/g;
								\$smiles=\~s/\\\)/\\\\\)/g;
								
								\$scan = \$dtafile;
								\$scan =~s/\(\.\*\)\\\.\(\\d+\)\\\.dta/\$2/;
							}
						}
						
						\$index++;		
					}	
				}
	################################################

				
				
	###### generate decoys 
				foreach my \$theoretical_formula (\@\$return_decoy) 
				{
					chomp \$theoretical_formula;
					my \$carbon_match = 0;
					my \$nitrogen_match = 0;
					
					if(\$params->{'labeled_ID_method'} == 2)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					elsif(\$params->{'labeled_ID_method'} == 3)
					{
						\$carbon_match = (\$theoretical_formula=~\/\$C_num\\D+\/);
						\$nitrogen_match = 1;
					}
					elsif(\$params->{'labeled_ID_method'} == 4)
					{
						\$carbon_match = 1;
						\$nitrogen_match = (\$theoretical_formula=~\/\$N_num\\D+\/);
					}
					else
					{
						print "please set the labeled_ID_method parameter [1-4]\n";
						exit(1);
					}					
						
					if(\$carbon_match == 1 and \$nitrogen_match ==1)
					{				
						
						my (\$formula_db,\$mass_db) = split(/\:/,\$theoretical_formula);
=head					
						my \$query = new Spiders::DatabaseQuery();

						\$query->setDatabase(\$params->{structure_database});
						\$query->setFormula(\$formula_db);
						\$query->setMass(\$mass_db);
						\$str_db_return = \$query->QueryStructureDatabase();				
=cut
						print MRESULT \$index,"\\t",\$mono_mass,"\\t",\$N_mz,"\\t",\$C_mz,"\\t",\$pair_score,"\\t",\$c_diff_score,"\\t",\$c_int_score,"\\t", \$c_similarity,"\\t", \$n_diff_score,"\\t",\$n_int_score,"\\t", \$n_similarity,"\\t", \$scan,"\\t",\$formula_comb->{\$mono_mass}->{'intensity'},"\\t","Decoy","\\t",\$formula_db,"\\n";					
						
						
						\$index++;		
					}	
				}
			}				
		}	
    }
	
}
	
	
	
EOF
	close(RUNSHELL);	
}


1;
