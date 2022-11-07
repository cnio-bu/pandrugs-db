##!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes;
use Switch;
use Cwd;
use File::Path;
use Net::FTP;
use POSIX ();
use LWP::Simple;
use File::Copy;
use DBM::Deep;
use Parallel::ForkManager;
use List::MoreUtils qw(uniq);

# Variable Initialization
my $outdir = "./output/";
my $dbdir;
my $outpath = getcwd;
my (%cosmic_list1, %cosmic_list2, %cosmic_list3, %cosmic_list4, %cosmic_list5, %cosmic_list6, %cosmic_list7, %cosmic_list8, %cosmic_list9, %cosmic_list10, %cosmic_list11, %cosmic_list12, %cosmic_list13, %cosmic_list14, %cosmic_list15, %cosmic_list16, %cosmic_list17);
my ($cosmic_list1, $cosmic_list2, $cosmic_list3, $cosmic_list4, $cosmic_list5, $cosmic_list6, $cosmic_list7, $cosmic_list8, $cosmic_list9, $cosmic_list10, $cosmic_list11, $cosmic_list12, $cosmic_list13, $cosmic_list14, $cosmic_list15, $cosmic_list16, $cosmic_list17);
my $vep_results_file;
my $logfile = "";
my %pathw_desc;
my $pathw_desc;
my %genes_ids;
my $genes_ids;
my $root_name = "";
my $jobid;
my $genes_affected;

my $conseq_file = 1;
my $MAX_PROCESSES = 8;

my %kegg_gene_pathway_DB;
my $kegg_gene_pathway_DB;

my ($start, $end, $time);
my (%pfam_a, %uniprot_b, %interpro_a);
my ($pfam_a, $uniprot_b, $interpro_a);

my %zygosity;
my %last_domain;
my %essential;
my %cancer_domain;
my %clinvar;
my %appris;

my $generole;
my $last_domain;
my $essential;
my $cancer_domain;
my $clinvar;
my $appris;

#Command line arguments handle
if (!@ARGV || grep (/^((\-\-help)|(\-h))$/,@ARGV)) {
	&help_info;
}

for my $a (0..$#ARGV){

	switch ($ARGV[$a]){

		#Input file
		case /^((\-\-vepfile=)|(\-f=))/ {
			$ARGV[$a] =~ /\-(\-vepfile|f)=(.+)/;
			$vep_results_file = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. -f=file.vcf\n\n";
            my @vep_results_file_tmp = glob ("$2");
            $vep_results_file = $vep_results_file_tmp[0];
		}

		# Output file
		case /^((\-\-output=)|(\-o=))/ {
			$ARGV[$a] =~ /\-(\-output|o)=(.+)/;
			$outdir = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. -o=/home/user/z13-222\n\n";
			$outpath = "";
		}

		# Root for the file
		case /^((\-\-root=)|(\-r=))/{
			$ARGV[$a] =~ /\-(\-root|r)=(.+)/;
			if ($2) {
				$root_name = $2 . "_";
			}
			else {
				die "\nEmpty argument. Please enter the parameter information.\n\neg. -r=analysis\n\n";
			}
		}

		# JobID
		case /^((\-\-int=)|(\-i=))/ {
			$ARGV[$a] =~ /\-(\-int|i)=(.+)/;
			$jobid = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. -i=20140213_000000\n\n";
		}

		# databases path
		case /^((\-\-databases=)|(\-d=))/ {
			$ARGV[$a] =~ /\-(\-databases|d)=(.+)/;
			$dbdir = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. -d=/home/epineiro/Programs/PCDA/databases\n\n";
			$dbdir = $dbdir . "/";
		}

		# Fork processes
		case /^(\-\-forkp=)/ {
			$ARGV[$a] =~ /\-(\-forkp)=(.+)/;
			$MAX_PROCESSES = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. --forkp=8\n\n";
		}
		else {
			die "\nArgument $ARGV[$a] not valid.\n\n";
		}

	}

}

if (!$vep_results_file) {
	die "\nVEP file not indicated. Please enter the absolute path to VEP sorted output file.\n\neg. -f=z13-222\n\n";
}
if (!$dbdir) {
	die "\nPath to databases not indicated. Please, enter the databases path.\n\neg. -d=databases\n\n";
}

# Create folders
mkpath($dbdir, 0);

# Start time counter
$start = Time::HiRes::gettimeofday();

# get new experiment id
if ($jobid) {
	$jobid = $root_name . $jobid . "_VEP";
}
else {
	$jobid = &get_runid;
}
my $outexp = $outdir . "/" . $jobid . "/";

# Create experiment folder in output folder
mkpath($outexp);

# Call VEP main subroutine
#print("$cosmic1\n");
&VEP_Parser_Csv;

$end = Time::HiRes::gettimeofday();
$time =  sprintf("%.2f", $end - $start);
printl ("\nTotal time: $time seconds\n");

# Log file creation
open LOGFILE, ">$outpath$jobid" . ".log" or die $!;
print LOGFILE $logfile;
close LOGFILE;

exit;

sub load_vars2 {
#Load files into variables

   	$cosmic_list1 = DBM::Deep->new("$dbdir/cosmic01.db");
   	$cosmic_list2 = DBM::Deep->new("$dbdir/cosmic02.db");
   	$cosmic_list3 = DBM::Deep->new("$dbdir/cosmic03.db");
   	$cosmic_list4 = DBM::Deep->new("$dbdir/cosmic04.db");
   	$cosmic_list5 = DBM::Deep->new("$dbdir/cosmic05.db");
   	$cosmic_list6 = DBM::Deep->new("$dbdir/cosmic06.db");
   	$cosmic_list7 = DBM::Deep->new("$dbdir/cosmic07.db");
   	$cosmic_list8 = DBM::Deep->new("$dbdir/cosmic08.db");
   	$cosmic_list9 = DBM::Deep->new("$dbdir/cosmic09.db");
   	$cosmic_list10 = DBM::Deep->new("$dbdir/cosmic10.db");
   	$cosmic_list11 = DBM::Deep->new("$dbdir/cosmic11.db");
   	$cosmic_list12 = DBM::Deep->new("$dbdir/cosmic12.db");
   	$cosmic_list13 = DBM::Deep->new("$dbdir/cosmic13.db");
   	$cosmic_list14 = DBM::Deep->new("$dbdir/cosmic14.db");
   	$cosmic_list15 = DBM::Deep->new("$dbdir/cosmic15.db");
   	$cosmic_list16 = DBM::Deep->new("$dbdir/cosmic16.db");
   	$cosmic_list17 = DBM::Deep->new("$dbdir/cosmic17.db");

	$genes_ids = DBM::Deep->new("$dbdir/genesids.db");

	$kegg_gene_pathway_DB = DBM::Deep->new("$dbdir/gene_pathway.db");

	$pathw_desc = DBM::Deep->new("$dbdir/pathways_desc.db");

	$pfam_a = DBM::Deep->new("$dbdir/pfam.db");

	$uniprot_b = DBM::Deep->new("$dbdir/uniprot_b.db");

	$interpro_a = DBM::Deep->new("$dbdir/interpro_a.db");
	$last_domain = DBM::Deep->new("$dbdir/last_domain.db");

	$depmap = DBM::Deep->new("$dbdir/depmap.db");

	$cancer_domain = DBM::Deep->new("$dbdir/cancer_domain.db");

	$clinvar = DBM::Deep->new("$dbdir/clinvar.db");

    $generole = DBM::Deep->new("$dbdir/generole.db");

}

sub help_info {

	print "\n\n--vepfile=filename or -f=filename \t\t Input file containing results of VEP from Ensembl analysis. Mandatory.\n\n";
	print "--output=directory or -o=directory \t\t Execution output dir. Default ./output.\n\n";
	print "--int=jobID or -i=jobID \t\t\t Job ID code (when executing from sequencingAP). Default: Generated during execution.\n\n";
	print "--databases=directory or -d=directory \t\t\t Absolute path to databases directory. Mandatory.\n\n";
	print "--forkp \t\t\t\t\t\t Number of forks to improve runtime. Default: 4.\n\n";

	print "\ni.e. VEP_parser.pl -f=file.vcf -o=/home/user/z13-222 -r=analysis -i=20140213_000000 -d=/home/epineiro/Programs/PCDA/databases --forkp=8\n\n";
	exit;

}

sub get_runid {

	my $timestamp = POSIX::strftime("%Y%m%d_%H%M%S", localtime);
	return $root_name . "$timestamp" . "_VEP";

}

sub VEP_Parser_Csv($$) {#Require a DB_conection and source data file

	# Main process
	# Variable inti
	my $data = '';

	# Open vcf file
	open (RFILE, "gunzip -c $vep_results_file |") || die "Could not find file $vep_results_file";
	printl ("\nProcessing file $vep_results_file...\n");
	my @rfile = <RFILE>;
	close RFILE;

	my %pos;

    my $count = 0;

	foreach my $i (0..$#rfile) {

		# Remove lines with ## & # to be able to handle file as csv
		$rfile[$i] =~ s/\r/\n/g;

		if ($rfile[$i] =~ /^#[^#]/) {
			$rfile[$i] =~ s/#//;
			$data = $data . $rfile[$i];
		}
		elsif ($rfile[$i] =~ /^[^#]/) {
			$data = $data . $rfile[$i];
            $count++;
		}
		elsif ($rfile[$i] =~ /^##INFO=<ID=CSQ.+Format: (.+)">/) {
			my @vep_fields = split ('\|', $1);
			foreach my $i (0..$#vep_fields) {
				$pos{Consequence} = $i if ($vep_fields[$i] eq "Consequence");
				$pos{Impact} = $i if ($vep_fields[$i] eq "IMPACT");
				$pos{Existing_variation} = $i if ($vep_fields[$i] eq "Existing_variation");
				$pos{Feature} = $i if ($vep_fields[$i] eq "Feature");
				$pos{PolyPhen} = $i if ($vep_fields[$i] eq "PolyPhen");
				$pos{SIFT} = $i if ($vep_fields[$i] eq "SIFT");
                $pos{CADD_PHRED} = $i if ($vep_fields[$i] eq "CADD_PHRED");
                $pos{CADD_RAW} = $i if ($vep_fields[$i] eq "CADD_RAW");
				$pos{SYMBOL} = $i if ($vep_fields[$i] eq "SYMBOL");
				$pos{Protein_position} = $i if ($vep_fields[$i] eq "Protein_position");
				$pos{Amino_acids} = $i if ($vep_fields[$i] eq "Amino_acids");
				$pos{HGVSc} = $i if ($vep_fields[$i] eq "HGVSc");
				$pos{HGVSp} = $i if ($vep_fields[$i] eq "HGVSp");
				$pos{GMAF} = $i if ($vep_fields[$i] eq "AF");
				$pos{CDS_position} = $i if ($vep_fields[$i] eq "CDS_position");
				$pos{Allele} = $i if ($vep_fields[$i] eq "Allele");
				$pos{Gene} = $i if ($vep_fields[$i] eq "Gene");
				$pos{Feature_type} = $i if ($vep_fields[$i] eq "Feature_type");
				$pos{cDNA_position} = $i if ($vep_fields[$i] eq "cDNA_position");
				$pos{Codons} = $i if ($vep_fields[$i] eq "Codons");
				$pos{VARIANT_CLASS} = $i if ($vep_fields[$i] eq "VARIANT_CLASS");
                $pos{gnomAD} = $i if ($vep_fields[$i] eq "gnomAD_AF");
                $pos{gnomAD_NFE} = $i if ($vep_fields[$i] eq "gnomAD_NFE_AF");
				$pos{EXON} = $i if ($vep_fields[$i] eq "EXON");
				$pos{APPRIS} = $i if ($vep_fields[$i] eq "APPRIS");

			}

		} else {
            
        }

	}

	# Save modifications to ensembl_vep.csv
    my $outexp_path = $outexp;
	$outexp_path =~ s/^\.//g;

	$outpath .= $outexp_path;
	my $outpathfile = $outpath . "ensembl_vep.csv";

	open (SFILE, ">$outpathfile") || die "Could not save temp file\n";
	print SFILE $data;
	close SFILE;

	printl ("\nensembl_vep.csv file created!\n");

	my $sth;

	my @last_gene = ("","","","");

    print "\n\nCreating annotations for $count variants...\n";

    open (INPUT, "$outexp/ensembl_vep.csv");
    my $ensemblhead = "";
    my $chr = "";
	my $filechr = "";
	my @chr_names = ();
	while (<INPUT>) {
		if ($_ =~ /^CHROM/) {
			$ensemblhead = $_;
		} else {
			my @line = split ("\t", $_);
            if ($line[0] ne $chr) {
				$chr = $line[0];
				push (@chr_names, $chr);
				if ($filechr ne "") {close $filechr}
				$filechr = "$outexp/ensembl_vep_$chr.csv";
				open (OUTPUTCHR, ">$filechr");
				print OUTPUTCHR $ensemblhead;
				print OUTPUTCHR $_;
			} else {
				print OUTPUTCHR $_;
			}
		}
	}
	close INPUT;

	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	DATA_LOOP:
	foreach $chr(@chr_names) {
        print "Processing chromosome $chr...\n";
		&load_vars2;
		my $pid = $pm->start and next DATA_LOOP;
        open (INPUT, "$outexp/ensembl_vep_$chr.csv");
       	open (OUT, ">$outexp/vep_data_$chr.csv");
   	   	open (OUTSORT, ">$outexp/vep_data_sorted_$chr.csv");

        print OUT "Chr\tLoc\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\tmut\tlocation\tallele\tgene\tfeature\tfeature_type\tconsequence\timpact\tcdna_position\tcds_position\tprotein_position\tamino_acids\tcodons\texisting_variation\textra\tprincipal\tpoly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\tCADD_raw\tgene_hgnc\tgene_role\tcosmic_id\tKEGG_data\tKEGG_path_id\tclinvar_acc\tclinvar_disease\tclinvar_clinical_significance\tvariation_type\tHGVS_cDNA\tHGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\tinterpro\tgene_cosmic_freq\tmut_cosmic_freq\tvscore\tbranch\tzygosity\n";

        print OUTSORT "chr\tloc\tmut\tgene\tfeature\tfeature_type\tconsequence\timpact\tprincipal\tpoly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\tCADD_raw\tgene_hgnc\tgene_role\tcosmic_id\tkegg_data\tkegg_path_id\tprotein_position\tamino_acids\tclinvar_acc\tclinvar_disease\tclinvar_clinical_significance\tvariation_type\tHGVS_cDNA\tHGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\tinterpro\tgene_cosmic_freq\tmut_cosmic_freq\tvscore\tbranch\tzygosity\n";

        $count = 1;
        while (<INPUT>) {
            unless ($_ =~ /^CHROM/) {
                chomp $_;
                my @line = split ("\t", $_);

    			$line[0] =~ s/chr//;
    			$line[2] = "" if ($line[2] eq ".");

                my $VCF_pos = $line[1];
                my $VCF_ref = $line[3];
                my $VCF_alt = $line[4];

			    if ((length($line[3]) ne length($line[4])) && ($line[4] !~ /,/)) {
				    if (length($line[3]) > length($line[4])) {
					    #Puede haber una deleción interna que no se procesaría correctamente
					    if ($line[3] =~ /^$line[4]/) {
						    $line[3] =~ s/^$line[4]//;
						    $line[1] = $line[1] + length($line[4]);
						    $line[4] = "-";
					    }
				    } else {
					    #Puede haber una inserción interna que no se procesaría correctamente (ej. TGCTCTACC/TATAGATCGGAAGCTCTACC)
					    if ($line[4] =~ /^$line[3]/) {
						    $line[4] =~ s/^$line[3]//;
						    $line[1] = $line[1] + length($line[3]);
						    $line[3] = "-";
					    }
				    }
			    }

    			$line[10] = $line[3] . "/" . $line[4];

    			$line[11] = $line[0] . ":" . $line[1];

				if ($line[9]) {
					my @gentype = split(":", $line[9]);
					if ($gentype[0] eq "1/1") {
	       				$zygosity{"$line[0]_$line[1]_$line[10]"} = "Homozygous";
				    }
				    else {
	       				$zygosity{"$line[0]_$line[1]_$line[10]"} = "Heterozygous";
				    }
				} else {
					$zygosity{"$line[0]_$line[1]_$line[10]"} = "";
				}

                my @q_data;

       			if ($line[7] =~ /CSQ=(.+)/) {
    				my @trans = split (",", $1);
    				foreach my $t (0..$#trans) {
	    				my $transcript = "";
	    				my $pol_cons = "";
	    				my $pol_sco = "";
	    				my $sift_cons = "";
	    				my $sift_sco = "";
                        my $CADD_phred = "";
                        my $CADD_raw = "";
	    				my $cosmic_id = "";
	    				my $cosmic_fathmm = "";
                        my $gene_freq = "";
                        my $mut_freq = "";
                        my $cosmic_total = "";
	    				my $kegg_data = "";
	    				my $kegg_ids = "";
	    				my $clinvar_acc = "";
	    				my $clinvar_dis = "";
	    				my $clinvar_pat = "";
    					my $HGVSc = "";
    					my $HGVSp = "";
    					my $GMAF = "";
    					my $uniprot = "";
    					my $pfam = "";
    					my $interpro = "";
						my $var_type = "";
	                    my $gnomAD = "";
						my $gnomAD_NFE = "";

						my @fields = split ('\|', $trans[$t]);

						$fields[$pos{Consequence}] =~ s/&/,/g;
						$fields[$pos{Existing_variation}] =~ s/&/,/g if ($fields[$pos{Existing_variation}]);

    					if ($fields[$pos{HGVSp}] && $fields[$pos{HGVSp}] =~ /:(.+)/) {
    						$HGVSp = $1;
    					}

    					if ($fields[$pos{HGVSc}] && $fields[$pos{HGVSc}] =~ /:(.+)/) {
    						$HGVSc = $1;
    					}

						if ($fields[$pos{PolyPhen}] && $fields[$pos{PolyPhen}] =~ /(\w+)\((\d+\.*\d*)\)/) {
							$pol_cons = $1;
							$pol_sco = $2;
						}

#						elsif ($fields[$pos{Consequence}] =~ /stop_gained/ || $fields[$pos{Consequence}] =~ /frameshift_variant/) {
#							$pol_cons = "inferred";
#							$pol_sco = 1;
#						}

    					if ($fields[$pos{SIFT}] && $fields[$pos{SIFT}] =~ /(\w+)\((\d+\.*\d*)\)/) {
    						$sift_cons = $1;
    						$sift_sco = $2;
    					}

#    					elsif ($fields[$pos{Consequence}] =~ /stop_gained/ || $fields[$pos{Consequence}] =~ /frameshift_variant/) {
#    						$sift_cons = "inferred";
#    						$sift_sco = 0;
#    					}

                        $CADD_phred = $fields[$pos{CADD_PHRED}];
                        $CADD_raw = $fields[$pos{CADD_RAW}];

  						if ($fields[$pos{HGVSc}]) {
							($cosmic_id, $cosmic_fathmm, $gene_freq, $mut_freq, $cosmic_total) = &chkmut_cosmic($fields[$pos{SYMBOL}], $fields[$pos{Feature}], $HGVSc);
							$cosmic_id .= ":$cosmic_fathmm" if ($cosmic_fathmm);
                            $gene_freq = $gene_freq if ($gene_freq);
                            $mut_freq = $mut_freq if ($mut_freq);
  						}

   						if ($fields[$pos{SYMBOL}]) {
                            if ($fields[$pos{SYMBOL}] ne $last_gene[0]) {
							    #Get information about gene symbol: pathway_description, pathway_ids and entrez_gene_id
       							@last_gene = get_kegg_id_sym(uc($fields[$pos{SYMBOL}]));
   	    					}
   		    				$kegg_data = $last_gene[1];
   			    			$kegg_ids = $last_gene[2];

                        }

						$var_type = $fields[$pos{"VARIANT_CLASS"}];

    					if ($fields[$pos{GMAF}]) {
    						my $num = $fields[$pos{GMAF}];
							my @GMAF_a = split ("&", $fields[$pos{GMAF}]);
							foreach my $gf (@GMAF_a) {
								if ($gf ne "") {
									$GMAF = $gf * 100;
								}
							}
    					}

                        if (exists($clinvar->{"$line[0]:$VCF_pos:$VCF_ref:$VCF_alt"})) {
     						$clinvar_acc = @{$clinvar->{"$line[0]:$VCF_pos:$VCF_ref:$VCF_alt"}}[1];
  	    					$clinvar_dis = @{$clinvar->{"$line[0]:$VCF_pos:$VCF_ref:$VCF_alt"}}[0];
   		    				$clinvar_pat = @{$clinvar->{"$line[0]:$VCF_pos:$VCF_ref:$VCF_alt"}}[2];
                        }

	    				my $prot_pos = $fields[$pos{Protein_position}];
	    				my $prot_end;
	    				my $ident = "";

                        if (exists($uniprot_b->{$fields[$pos{SYMBOL}]})) {
   							$ident = $uniprot_b->{$fields[$pos{SYMBOL}]};

   							my $prot_end = 0;

   							if ($fields[$pos{Protein_position}] =~ /(\d+)\-(\d+)/) {
   								$prot_pos = $1;
   								$prot_end = $2;
   							}

   							if ($fields[$pos{Protein_position}] =~ /(\d+)\-(\?)/) {
   								$prot_pos = $1;
   								$prot_end = 0;
   							}


   							if ($prot_pos ne "") {
                                if (exists($pfam_a->{$ident})) {
                                    foreach my $ia (0..scalar(@{$pfam_a->{$ident}})-1) {
    									if (($prot_pos >= ${@{$pfam_a->{$ident}}[$ia]}[2] && $prot_pos <= ${@{$pfam_a->{$ident}}[$ia]}[3]) || ($prot_end >= ${@{$pfam_a->{$ident}}[$ia]}[2] && $prot_end <= ${@{$pfam_a->{$ident}}[$ia]}[3])) {
   											$pfam = "${@{$pfam_a->{$ident}}[$ia]}[0]: ${@{$pfam_a->{$ident}}[$ia]}[1]";
                                        }
   									}
   								}
                                if (exists($interpro_a->{$ident})) {
                                    foreach my $ii (0..scalar(@{$interpro_a->{$ident}})-1) {
   										if (($prot_pos >= ${@{$interpro_a->{$ident}}[$ii]}[2] && $prot_pos <= ${@{$interpro_a->{$ident}}[$ii]}[3]) || ($prot_end >= ${@{$interpro_a->{$ident}}[$ii]}[2] && $prot_end <= ${@{$interpro_a->{$ident}}[$ii]}[3])) {
   											$interpro = "${@{$interpro_a->{$ident}}[$ii]}[0]: ${@{$interpro_a->{$ident}}[$ii]}[1]";
   										}
   									}
   								}
   								if ($fields[$pos{Consequence}] =~ /(stop_gained|frameshift_variant)/ && $interpro eq "" && exists($last_domain->{$ident}) && $prot_pos <= $last_domain->{$ident}) {
   									$interpro = "Mutation previous last protein domain";
   								}
   							}

	    				}

	    				my $gene_role = '';
	    				$gene_role = $generole->{$fields[$pos{SYMBOL}]} if ($generole->{$fields[$pos{SYMBOL}]});

						if ($fields[$pos{gnomAD}]) {
							my @gnomAD = split ("&", $fields[$pos{gnomAD}]);
							foreach my $ex (@gnomAD) {
								if ($ex ne "") {
									$gnomAD = $ex * 100;
								}
							}
						}
						if ($fields[$pos{gnomAD_NFE}]) {
							my @gnomAD = split ("&", $fields[$pos{gnomAD_NFE}]);
							foreach my $ex (@gnomAD) {
								if ($ex ne "") {
									$gnomAD_NFE = $ex * 100;
								}
							}
						}
		               	my @q_data_line = [$line[0], $line[1], $line[2], $line[3], $line[4], $line[5], $line[6], $line[7], $line[8], $line[9], $line[10], $line[11], $fields[$pos{Allele}], $fields[$pos{Gene}], $fields[$pos{Feature}], $fields[$pos{Feature_type}], $fields[$pos{Consequence}], $fields[$pos{Impact}], $fields[$pos{cDNA_position}], $fields[$pos{CDS_position}], $fields[$pos{Protein_position}], $fields[$pos{Amino_acids}], $fields[$pos{Codons}], $fields[$pos{Existing_variation}], "", $fields[$pos{APPRIS}], $pol_cons, $pol_sco, $sift_cons, $sift_sco, $CADD_phred, $CADD_raw, $fields[$pos{SYMBOL}], $gene_role, $cosmic_id, $kegg_data, $kegg_ids, $clinvar_acc, $clinvar_dis, $clinvar_pat, $var_type, $HGVSc, $HGVSp, $GMAF, $gnomAD, $gnomAD_NFE, $pfam, $interpro, $gene_freq, $mut_freq];

						@q_data_line = &create_vscore (@q_data_line);

                		no warnings 'uninitialized';
        				push @q_data, @q_data_line;
                        print OUT join("\t", @{$q_data[scalar(@q_data)-1]}), "\n";
                    }
                    $count++;
				}

                @q_data = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[10] cmp $b->[10] || $a->[13] cmp $b->[13] || $a->[42] cmp $b->[42] || $a->[14] cmp $b->[14]} @q_data;

                foreach my $qi (0..$#q_data) {
            		no warnings 'uninitialized';
                    if ($conseq_file == 1) {
                        if ($q_data[$qi][17] =~ /HIGH|MODERATE/) {
                            print OUTSORT "$q_data[$qi][0]\t$q_data[$qi][1]\t$q_data[$qi][10]\t$q_data[$qi][13]\t$q_data[$qi][14]\t$q_data[$qi][15]\t$q_data[$qi][16]\t$q_data[$qi][17]\t$q_data[$qi][25]\t$q_data[$qi][26]\t$q_data[$qi][27]\t$q_data[$qi][28]\t$q_data[$qi][29]\t$q_data[$qi][30]\t$q_data[$qi][31]\t$q_data[$qi][32]\t$q_data[$qi][33]\t$q_data[$qi][34]\t$q_data[$qi][35]\t$q_data[$qi][36]\t$q_data[$qi][20]\t$q_data[$qi][21]\t$q_data[$qi][37]\t$q_data[$qi][38]\t$q_data[$qi][39]\t$q_data[$qi][40]\t$q_data[$qi][41]\t$q_data[$qi][42]\t$q_data[$qi][43]\t$q_data[$qi][44]\t$q_data[$qi][45]\t$q_data[$qi][46]\t$q_data[$qi][47]\t$q_data[$qi][48]\t$q_data[$qi][49]\t$q_data[$qi][50]\t$q_data[$qi][51]\t$q_data[$qi][52]\t$q_data[$qi][53]\t$q_data[$qi][54]\t$q_data[$qi][55]\t$q_data[$qi][56]\t$q_data[$qi][57]\t$q_data[$qi][58]\t$q_data[$qi][59]\t$q_data[$qi][60]\t$q_data[$qi][61]\n";
                        }
                    } else {
                        print OUTSORT "$q_data[$qi][0]\t$q_data[$qi][1]\t$q_data[$qi][10]\t$q_data[$qi][13]\t$q_data[$qi][14]\t$q_data[$qi][15]\t$q_data[$qi][16]\t$q_data[$qi][17]\t$q_data[$qi][25]\t$q_data[$qi][26]\t$q_data[$qi][27]\t$q_data[$qi][28]\t$q_data[$qi][29]\t$q_data[$qi][30]\t$q_data[$qi][31]\t$q_data[$qi][32]\t$q_data[$qi][33]\t$q_data[$qi][34]\t$q_data[$qi][35]\t$q_data[$qi][36]\t$q_data[$qi][37]\t$q_data[$qi][20]\t$q_data[$qi][21]\t$q_data[$qi][38]\t$q_data[$qi][39]\t$q_data[$qi][40]\t$q_data[$qi][41]\t$q_data[$qi][42]\t$q_data[$qi][43]\t$q_data[$qi][44]\t$q_data[$qi][45]\t$q_data[$qi][46]\t$q_data[$qi][47]\t$q_data[$qi][48]\t$q_data[$qi][49]\t$q_data[$qi][50]\t$q_data[$qi][51]\t$q_data[$qi][52]\t$q_data[$qi][53]\t$q_data[$qi][54]\t$q_data[$qi][55]\t$q_data[$qi][56]\t$q_data[$qi][57]\t$q_data[$qi][58]\t$q_data[$qi][59]\t$q_data[$qi][60]\t$q_data[$qi][61]\n";
                    }
                }
            }
        }
        close INPUT;
        close OUT;
        close OUTSORT;
        print "Chromosome $chr processed!\n";
		$pm->finish; # Terminates the child process
	}
	$pm->wait_all_children;
	open (GREATOUT, ">$outexp/vep_data.csv");
    print GREATOUT "Chr\tLoc\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\tmut\tlocation\tallele\tgene\tfeature\tfeature_type\tconsequence\timpact\tcdna_position\tcds_position\tprotein_position\tamino_acids\tcodons\texisting_variation\textra\tprincipal\tpoly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\tCADD_raw\tgene_hgnc\tgene_role\tcosmic_id\tKEGG_data\tKEGG_path_id\tclinvar_acc\tclinvar_disease\tclinvar_clinical_significance\tvariation_type\tHGVS_cDNA\tHGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\tinterpro\tgene_cosmic_freq\tmut_cosmic_freq\tvscore\tbranch\n";

    foreach $chr (@chr_names) {
		open (INPUT, "$outexp/vep_data_$chr.csv");
		while (<INPUT>) {
			unless ($_ =~ /^Chr/) {
				print GREATOUT $_;
			}
		}
		close INPUT;
		unlink ("$outexp/vep_data_$chr.csv");
		unlink ("$outexp/ensembl_vep_$chr.csv");
	}
	close GREATOUT;

	open (GREATOUT, ">$outexp/vep_data_sorted.csv");
    print GREATOUT "chr\tloc\tmut\tgene\tfeature\tfeature_type\tconsequence\timpact\tprincipal\tpoly_effect\tpoly_score\tsift_effect\tsift_score\tCADD_phred\tCADD_raw\tgene_hgnc\tgene_role\tcosmic_id\tkegg_data\tkegg_path_id\tprotein_position\tamino_acids\tclinvar_acc\tclinvar_disease\tclinvar_clinical_significance\tvariation_type\tHGVS_cDNA\tHGVS_protein\tGMAF\tgnomAD\tgnomAD_NFE\tpfam\tinterpro\tgene_cosmic_freq\tmut_cosmic_freq\tvscore\tbranch\n";

    foreach $chr (@chr_names) {
		open (INPUT, "$outexp/vep_data_sorted_$chr.csv");
		while (<INPUT>) {
			unless ($_ =~ /^chr/) {
				print GREATOUT $_;
			}
		}
		close INPUT;
		unlink ("$outexp/vep_data_sorted_$chr.csv")
	}
	close GREATOUT;

	# Save query to file
	my $veptable_sort = $outpath . "vep_data_sorted.csv";

	# Select the highest score for the principal isoform
	my @vep_file;
	open FILE, "<$veptable_sort" or die $!;
	while (<FILE>) {
		chomp $_;
		my @line = split ("\t", $_);
		push @vep_file, [@line];
	}
	close FILE;

	my %isoform;
	foreach my $is (1..$#vep_file) {
		if ($vep_file[$is][7] =~ /HIGH|MODERATE/) {
			if (exists($isoform{"$vep_file[$is][15]"})) {
				$isoform{"$vep_file[$is][15]"}[0] = 1 if ($vep_file[$is][8] eq "P1");
				$isoform{"$vep_file[$is][15]"}[1] = 1 if ($vep_file[$is][8] eq "P2");
				$isoform{"$vep_file[$is][15]"}[2] = 1 if ($vep_file[$is][8] eq "P3");
				$isoform{"$vep_file[$is][15]"}[3] = 1 if ($vep_file[$is][8] eq "P4");
				$isoform{"$vep_file[$is][15]"}[4] = 1 if ($vep_file[$is][8] eq "P5");
				$isoform{"$vep_file[$is][15]"}[5] = 1 if ($vep_file[$is][8] eq "A1" || $vep_file[$is][8] eq "A2" || $vep_file[$is][8] eq "");
			} else {
				$isoform{"$vep_file[$is][15]"} = [1, 0, 0, 0, 0, 0] if ($vep_file[$is][8] eq "P1");
				$isoform{"$vep_file[$is][15]"} = [0, 1, 0, 0, 0, 0] if ($vep_file[$is][8] eq "P2");
				$isoform{"$vep_file[$is][15]"} = [0, 0, 1, 0, 0, 0] if ($vep_file[$is][8] eq "P3");
				$isoform{"$vep_file[$is][15]"} = [0, 0, 0, 1, 0, 0] if ($vep_file[$is][8] eq "P4");
				$isoform{"$vep_file[$is][15]"} = [0, 0, 0, 0, 1, 0] if ($vep_file[$is][8] eq "P5");
				$isoform{"$vep_file[$is][15]"} = [0, 0, 0, 0, 0, 1] if ($vep_file[$is][8] eq "A1" || $vep_file[$is][8] eq "A2" || $vep_file[$is][8] eq "");
			}
		}
	}

	foreach my $key (keys (%isoform)) {
		if ($isoform{$key}[0] == 1) {
			$isoform{$key} = "P1";
		} elsif ($isoform{$key}[1] == 1) {
			$isoform{$key} = "P2";
		} elsif ($isoform{$key}[2] == 1) {
			$isoform{$key} = "P3";
		} elsif ($isoform{$key}[3] == 1) {
			$isoform{$key} = "P4";
		} elsif ($isoform{$key}[4] == 1) {
			$isoform{$key} = "P5";
		} else {
			$isoform{$key} = "|A1|A2";
		}
	}

	my %genes_affected;
	my $head = "";
	foreach my $is (0..$#vep_file) {
		if ($is == 0) {
			$head = join("\t", @{$vep_file[$is]});
		} else {
			if ($vep_file[$is][7] =~ /HIGH|MODERATE/) {
				if ($vep_file[$is][8] =~ /^$isoform{$vep_file[$is][15]}$/) {
					if (exists($genes_affected{$vep_file[$is][15]})) {
						if ($vep_file[$is][35] >= $genes_affected{$vep_file[$is][15]}[0]) {
							$genes_affected{$vep_file[$is][15]} = [$vep_file[$is][35], $vep_file[$is][36]];
							push @{$genes_affected{$vep_file[$is][15]}}, @{$vep_file[$is]};
						}
					} else {
						$genes_affected{$vep_file[$is][15]} = [$vep_file[$is][35], $vep_file[$is][36]];
						push @{$genes_affected{$vep_file[$is][15]}}, @{$vep_file[$is]};
					}
				}
			}
		}
	}

	$genes_affected = $outpath . "genes_affected.csv";

	open (FILE, ">$genes_affected");
	print FILE "gene_hgnc\tmax(vscore)\tbranch\t$head\n";
	foreach my $key (keys(%genes_affected)) {
		print FILE $key, "\t", join ("\t", @{$genes_affected{$key}}), "\n";
	}
	close FILE;

	printl ("\ngenes_affected.csv file created!\n\n");

}

sub chkmut_cosmic() {

	my ($gene, $transcript, $HGVSc) = @_;

    if (exists($cosmic_list1->{"$gene:$transcript:$HGVSc"})) {
        return (@{$cosmic_list1->{"$gene:$transcript:$HGVSc"}});
	} elsif (exists($cosmic_list2->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list2->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list3->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list3->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list4->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list4->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list5->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list5->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list6->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list6->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list7->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list7->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list8->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list8->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list9->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list9->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list10->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list10->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list11->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list11->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list12->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list12->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list13->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list13->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list14->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list14->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list15->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list15->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list16->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list16->{"$gene:$transcript:$HGVSc"}});
    } elsif (exists($cosmic_list17->{"$gene:$transcript:$HGVSc"})){
        return (@{$cosmic_list17->{"$gene:$transcript:$HGVSc"}});
    } 

}

sub get_kegg_id_sym {

    my ($hgnc_symbol) = @_;

    if (exists($genes_ids->{$hgnc_symbol})) {
		my @kegg_pathways = &get_kegg_path("$genes_ids->{$hgnc_symbol}");
		return ($hgnc_symbol,@kegg_pathways,$1);
	}
	else {
		return ("","","","");
	}

}

sub get_kegg_path {
	my ($kegg_id) = @_;

    if (exists($kegg_gene_pathway_DB->{$kegg_id})) {
		my @kegg_path_desc = split /\|/, $kegg_gene_pathway_DB->{$kegg_id};
		return get_kegg_path_desc(\@kegg_path_desc);
	}
	else{
		return ("","");
	}
}

sub get_kegg_path_desc {
    my @kegg_paths = @{$_[0]};
	my $results = "";

	foreach my $path_id (@kegg_paths) {
		if (exists($pathw_desc->{$path_id})) {
			$results .= "$pathw_desc->{$path_id}|";
		}
	}

	my $keggpaths = join ("|",@kegg_paths);

	$results =~ s/\|$//g;

	return ($results, $keggpaths);
}

sub create_vscore() {

	# vscore creation
	my @vep_file = @_;
	my $score = 0;
	my $gene_freq = 0;
	my $cosmic_freq = 0;
	my $lastgene = "";
	my $skip = 0;
	my %scored_columns = ();

	# Get column position to apply the score (poly_score, sift_score, cosmic_id, kegg_path_id )
	# Get Gene_hgnc column for getting sample mutation frequency
	# Add column CCLE genes matches column
	my @linedata = ("Chr", "Loc", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample", "mut", "location", "allele", "gene", "feature", "feature_type", "consequence", "impact", "cdna_position", "cds_position", "protein_position", "amino_acids", "codons", "existing_variation", "extra", "principal", "poly_effect", "poly_score", "sift_effect", "sift_score", "CADD_phred", "CADD_raw", "gene_hgnc", "gene_role", "cosmic_id", "KEGG_data", "KEGG_path_id", "clinvar_acc", "clinvar_disease", "clinvar_clinical_significance", "variation_type", "HGVS_cDNA", "HGVS_protein", "GMAF", "gnomAD", "gnomAD_NFE", "pfam", "interpro", "gene_freq", "mut_freq");

	foreach my $i (0..$#linedata) {

		if ($linedata[$i] eq "Chr") { $scored_columns{Chr} = $i; }
		if ($linedata[$i] eq "Loc") { $scored_columns{Loc} = $i; }
		if ($linedata[$i] eq "mut") { $scored_columns{mut} = $i; }
		if ($linedata[$i] eq "poly_score") { $scored_columns{poly_score} = $i; }
		if ($linedata[$i] eq "sift_score") { $scored_columns{sift_score} = $i; }
		if ($linedata[$i] eq "CADD_phred") { $scored_columns{CADD_phred} = $i; }
		if ($linedata[$i] eq "gene_role") { $scored_columns{gene_role} = $i; }
		if ($linedata[$i] eq "cosmic_id") { $scored_columns{cosmic_id} = $i; }
		if ($linedata[$i] eq "protein_position") { $scored_columns{protein_position} = $i; }
		if ($linedata[$i] eq "mut_freq") { $scored_columns{mut_cosmic_freq} = $i; }
		if ($linedata[$i] eq "gene_freq") { $scored_columns{gene_cosmic_freq} = $i; }

		if ($linedata[$i] eq "consequence") { $scored_columns{consequence} = $i; }
		if ($linedata[$i] eq "impact") { $scored_columns{impact} = $i; }
		if ($linedata[$i] eq "GMAF") { $scored_columns{GMAF} = $i; }
		if ($linedata[$i] eq "gnomAD") { $scored_columns{gnomAD} = $i; }
		if ($linedata[$i] eq "pfam") { $scored_columns{pfam} = $i; }
		if ($linedata[$i] eq "interpro") { $scored_columns{interpro} = $i; }
		if ($linedata[$i] eq "clinvar_clinical_significance") { $scored_columns{clinvar} = $i; }
		if ($linedata[$i] eq "INFO") {$scored_columns{INFO} = $i; }
		if ($linedata[$i] eq "gene_hgnc") { $scored_columns{gene_hgnc} = $i; }

	}

	foreach my $i (0..$#vep_file) {

		@linedata = @{$vep_file[$i]};

		if ($lastgene ne $linedata[$scored_columns{gene_hgnc}]) {

			$lastgene = $linedata[$scored_columns{gene_hgnc}];

		}

		if ($skip == 0) {

			#Decide the branch for the calculation according to the variation or gene definition
			my $branch;

			my @components = split ("; ", $linedata[$scored_columns{gene_role}]);
			my @genetype = ();
			foreach my $ic (0..$#components) {
				if ($components[$ic] =~ /:([\w ]+)/) {
                    if ($1 eq "ONC" || $1 eq "TSG") {
                        push(@genetype, $1);
                    }
                }
			}
            @genetype = uniq @genetype;

            my $genetype_string = join(":", @genetype);

			if ($genetype_string) {
				if ($genetype_string eq "ONC") {
					$branch = "ONC";
				} elsif ($genetype_string eq "TSG") {
					$branch = "TSG";
				} else {
					$branch = "UNCLASSIFIED";
				}
			} else {
				$branch = "UNCLASSIFIED";
			}

			if ($branch eq "UNCLASSIFIED") {
				if ($linedata[$scored_columns{consequence}] =~ /stop_gain|stop_lost|frameshift_variant|splice_donor_variant|splice_acceptor_variant|splice_region_variant/) {
					$branch = "TSG";
				}
			}

			# Add scores
			my $prediction_damaging = 0;
			
			#Cosmic ID
			if ($linedata[$scored_columns{cosmic_id}] =~ /(COSV\d+):*/) {
				if ($linedata[$scored_columns{cosmic_id}] =~ /PATHOGENIC/) {
					$prediction_damaging += 1;
				}
				if ($linedata[$scored_columns{mut_cosmic_freq}] && $linedata[$scored_columns{mut_cosmic_freq}] =~ /(\d+) \/ (\d+)/) {
					if ($branch eq "ONC" || $branch eq "UNCLASSIFIED") {
						my $ss;
						if ($1 >= 100) {
							$ss = 0.125 / 2;
						} else {
							$ss = (0.125 / 2) * ((log($1) - 0) / (log($2) - 0));
						}
						$score += $ss;
					}
				}
				if ($linedata[$scored_columns{gene_cosmic_freq}] && $linedata[$scored_columns{gene_cosmic_freq}] =~ /(\d+) \/ (\d+)/) {
					my $ss;
					if ($1 >= 100) {
						$ss = 0.125 / 2;
					} else {
						$ss = (0.125 / 2) * ((log($1) - 0) / (log($2) - 0));
					}
					$score += $ss;
				}
			}
			else {
				$cosmic_freq = 0;
			}
#print "cosm:$score\n";
			#Prediction score
			if ($linedata[$scored_columns{poly_score}] && $linedata[$scored_columns{poly_score}] > 0.435) {
				$prediction_damaging += 1;
			}

			if ($linedata[$scored_columns{sift_score}] ne "" && $linedata[$scored_columns{sift_score}] <= 0.05) {
				$prediction_damaging += 1;
			}

			if ($linedata[$scored_columns{CADD_phred}] && $linedata[$scored_columns{CADD_phred}] > 20) {
				$prediction_damaging += 1;
			}

			if ($prediction_damaging >= 3) {
				$score += 0.125;
			} elsif ($prediction_damaging == 2) {
				$score += 0.08;
			} elsif ($prediction_damaging == 1) {
				$score += 0.04
			}
#print "pred:$score\n";
			#Mutation type
			if ($linedata[$scored_columns{impact}] =~ /(HIGH)/) {
				$score += 0.125;
			}
#print "mut:$score\n";
			#Frequencies
			if ($linedata[$scored_columns{GMAF}] ne "") {
				if ($linedata[$scored_columns{GMAF}] < 1) {
					$score += 0.125 / 2;
				}
			} else {
				$score = ($score + 0.125 / 2);
			}
			if ($linedata[$scored_columns{gnomAD}] ne "") {
				if ($linedata[$scored_columns{gnomAD}] < 1) {
					$score += 0.125 / 2;
				}
			} else {
				$score = ($score + 0.125 / 2);
			}
#print "freq:$score\n";
			#Domains
			if ($linedata[$scored_columns{pfam}] =~ /^(\w+)\./ && exists($cancer_domain->{$1})) {
				$score += 0.125;
			} elsif ($linedata[$scored_columns{interpro}] eq "Mutation previous last protein domain") {
				$score += 0.125;
			} else {
				if ($linedata[$scored_columns{pfam}] || $linedata[$scored_columns{interpro}]) {
					$score += 0.125 / 2;
				}
			}
#print "dom:$score\n";
			#Clinvar
			if ($linedata[$scored_columns{clinvar}] =~ /Pathogenic/) {
				if ($zygosity{"$linedata[$scored_columns{Chr}]_$linedata[$scored_columns{Loc}]_$linedata[$scored_columns{mut}]"} ne "") {
					$score += 0.125;
				} else {
					if ($branch eq "ONC" || $branch eq "UNCLASSIFIED") {
						$score += 0.250;
					} else {
						$score += 0.3125;
					}
				}
			}
#print "clin:$score\n";
			#Homozigous
			if ($zygosity{"$linedata[$scored_columns{Chr}]_$linedata[$scored_columns{Loc}]_$linedata[$scored_columns{mut}]"} eq "Homozygous") {
				if ($branch eq "ONC" || $branch eq "UNCLASSIFIED") {
					$score += 0.125;
				} else {
					$score += 0.1875;
				}
			}
#print "hom:$score\n";
			#DepMap
			

			push @{$vep_file[$i]}, (sprintf("%.4f", $score), $branch);

		}

		$score = 0;

	}

    return @vep_file;
}

sub printl {
	$logfile = $logfile . $_[0];
	print $_[0];
}
