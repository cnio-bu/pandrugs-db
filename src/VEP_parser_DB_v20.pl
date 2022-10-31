##!/usr/bin/perl
use lib "modules/";
use strict;
use warnings;
use Time::HiRes;
use Switch;
use Cwd;
use File::Path;
use Net::FTP;
use Archive::Extract ();
use POSIX ();
use LWP::Simple;
use File::Copy;
use DBM::Deep;
use Parallel::ForkManager;

# Variable Initialization
my $dbdir;
my $genesids = "custom";
my $pathways = "pathway_desc.tsv";
my $logfile = "";

my ($start, $end, $time);

#Command line arguments handle
if (!@ARGV || grep (/^((\-\-help)|(\-h))$/,@ARGV)) {
	&help_info;
}

for my $a (0..$#ARGV){

	switch ($ARGV[$a]){

		# databases path
		case /^((\-\-databases=)|(\-d=))/ {
			$ARGV[$a] =~ /\-(\-databases|d)=(.+)/;
			$dbdir = $2 ? $2 : die "\nEmpty argument. Please enter the parameter information.\n\neg. -d=/home/epineiro/Programs/PCDA/databases\n\n";
			$dbdir = $dbdir . "/";
		}

		else {
			die "\nArgument $ARGV[$a] not valid.\n\n";
		}

	}

}

if (!$dbdir) {
	die "\nPath to databases not indicated. Please, enter the databases path.\n\neg. -d=databases\n\n";
}

# Create folders
mkpath($dbdir, 0);

&create_dbs;

# Start time counter
$start = Time::HiRes::gettimeofday();

$end = Time::HiRes::gettimeofday();
$time =  sprintf("%.2f", $end - $start);
printl ("\nTotal time: $time seconds\n");

exit;

sub create_dbs {
	
#	Load files into variables
	print "\n\nLoading database files...\n\n";

	my (%pfam_a, %interpro_a, %last_domain, %cancer_domain);

	my @cosmic_files = glob("$dbdir/cosmic*.tsv");

	foreach (@cosmic_files) {
		my $file = $_;
		$_ =~ s/.tsv/.db/;
		my $cosmic_list = DBM::Deep->new($_);
	print("$file\n");
		open (FILE, "<$file") or die "Couldn't open file: $!"; 
		while (<FILE>){
			chomp $_;
			if ($. % 100000 == 0) {print("$.\n")};
			my @line = split ("\t", $_);
			$cosmic_list->{$line[0]} = [$line[1], $line[2], "$line[3] / $line[5]", "$line[4] / $line[5]"];
		}
		close FILE;
	}
    
	print "COSMIC loaded!\n";

	my $genes_ids = DBM::Deep->new("$dbdir/genesids.db");
	my $genes_ids_file = $dbdir . $genesids;
	open (FILE, "<$genes_ids_file") or die "Couldn't open file: $!"; 
	while (<FILE>){
		chomp $_;
		my @line = split ("\t", $_);
		$genes_ids -> {$line[0]} = $line[1] if ($line[1]);
	}
	close FILE;

	print "genes IDs loaded!\n";

	my $kegg_gene_pathway_DB = DBM::Deep->new("$dbdir/gene_pathway.db");
	my $gene_pathway_file = $dbdir . "gene_pathway.tsv";
	open FILE, "<$gene_pathway_file" or die "Couldn't open file: $!"; 
	while (<FILE>){
		chomp $_;
		my @line = split("\t", $_);
		$kegg_gene_pathway_DB -> {$line[0]} = $line[1];
	}
	close FILE;

	print "gene-pathway loaded!\n";

	my $pathw_desc = DBM::Deep->new("$dbdir/pathways_desc.db");
	my $pathway_desc = $dbdir . $pathways;
	open FILE, "<$pathway_desc" or die "Couldn't open file: $pathway_desc $!";
	while (<FILE>){
		chomp $_;
		my @line = split("\t", $_);
		$pathw_desc -> {$line[0]} = $line[1];
	}
	close FILE;

	print "pathway description loaded!\n";

	my $pfam_a = DBM::Deep->new("$dbdir/pfam.db");
    my @pfam_file = glob("$dbdir/Pfam-A.full_*.tsv");
	open FILE, "<$pfam_file[0]" or die "Couldn't open file: $!";
	while (<FILE>){
		chomp ($_);
		my @line = split("\t", $_);
		if (exists($pfam_a{$line[4]})) {
			push @{$pfam_a{$line[4]}}, [$line[1], $line[2], $line[5], $line[6]];
		} else {
			@{$pfam_a->{$line[4]}} = [$line[1], $line[2], $line[5], $line[6]];
		}
	}
	close FILE;

	print "pfam loaded!\n";

	my $uniprot_b = DBM::Deep->new("$dbdir/uniprot_b.db");
    my @uniprot_file = "$dbdir/Uniprot.tsv";
	open FILE, "<$uniprot_file[0]" or die "Couldn't open file: $!";
	while (<FILE>) {
		chomp ($_);
		my @line = split("\t", $_);
		my $name = $1 if ($line[0] =~ /^([A-Z0-9]+)/);
		$uniprot_b->{$line[1]} = $name if ($line[1] ne "");
	}
	close FILE;

	print "uniprot loaded!\n";

	my $interpro_a = DBM::Deep->new("$dbdir/interpro_a.db");
    my @interpro_file = "$dbdir/Interpro.tsv";
	my $last_domain = DBM::Deep->new("$dbdir/last_domain.db");
	open FILE, "<$interpro_file[0]" or die "Couldn't open file: $!";
	while (<FILE>) {
		chomp ($_);
		my @line = split ("\t",$_);
		if (exists($interpro_a{$line[3]})) {
			push @{$interpro_a{$line[3]}}, [$line[0], $line[1], $line[4], $line[5]];
		} else {
			@{$interpro_a->{$line[3]}} = [$line[0], $line[1], $line[4], $line[5]];
		}
		if (exists($last_domain{$line[3]})) {
			if ($last_domain{$line[3]} < $line[4]) {
				$last_domain{$line[3]} = $line[4];
			}
		} else {
			$last_domain->{$line[3]} = $line[4];
		}
	}
	close FILE;

	print "interpro loaded!\n";

	my $oncorole = DBM::Deep->new("$dbdir/generole.db");
	open FILE, "<$dbdir/processed/generole.tsv" or die "Couldn't open file: $!";
    my %pos;
	while (<FILE>) {
		chomp ($_);
		my @line = split ("\t",$_);
		if ($_ =~ /^gene/) {
			for my $i (0..$#line) {
				$pos{$i} = $line[$i];
			}
		} else {
			my @roles = @line[1..$#line];
			my @role_list;
			for my $i (0 .. $#roles) {
				my $role = $roles[$i];
				if ($role ne "") {
					push(@role_list,"$pos{$i+1}:$role");
				}
			}
			$oncorole->{$line[0]} = join ("; ", @role_list);
		}
	}
	close FILE;
	print "Gene Role loaded!\n";

	my $essential = DBM::Deep->new("$dbdir/essential.db");
	open (ESSEN, "$dbdir/gene_essentiality_score.tsv");
	while (<ESSEN>) {
		chomp $_;
		my @line = split ("\t", $_);
		unless ($line[0] eq "symbol") {
			$essential->{$line[0]} = $line[1] ;
		}
	}
	close ESSEN;
	print "essential loaded!\n";

	my $cancer_domain = DBM::Deep->new("$dbdir/cancer_domain.db");
	open DOM, "<$dbdir/domains.tsv" or die "Couldn't open file: $!";
	while (<DOM>){
		chomp ($_);
		my @line = split ("\t",$_);
		unless (exists($cancer_domain{$line[4]})) {

			$cancer_domain->{$line[4]} = "";
		}
	}
	close DOM;

	print "cancer domains loaded!\n";

	my $clinvar = DBM::Deep->new("$dbdir/clinvar.db");
	my @clinvar_file = "$dbdir/Clinvar.tsv";
	open CLINVAR, "<$clinvar_file[0]" or die "Couldn't open file: $!";
	while (<CLINVAR>) {
		chomp $_;
		my @line = split ("\t", $_);
		if ($line[1] eq "GRCh38") {
			if (exists($clinvar->{"$line[2]:$line[3]:$line[4]:$line[5]"})) {
				@{$clinvar->{"$line[2]:$line[3]:$line[4]:$line[5]"}}[0] .= "; $line[7]";
				@{$clinvar->{"$line[2]:$line[3]:$line[4]:$line[5]"}}[1] .= "; $line[0]";
				@{$clinvar->{"$line[2]:$line[3]:$line[4]:$line[5]"}}[2] .= "; $line[8]";
			} else {
				$clinvar->{"$line[2]:$line[3]:$line[4]:$line[5]"} = [$line[7], $line[0], $line[8]];
			}
		}
	}
	close CLINVAR;

	print "clinvar loaded!\n";
}

sub printl {
	$logfile = $logfile . $_[0];
	print $_[0];
}

sub help_info {

	print "--databases=directory or -d=directory \t\t\t Absolute path to databases directory. Mandatory.\n\n";

	print "\ni.e. VEP_parser.pl -d=/home/epineiro/Programs/PCDA/databases\n\n";
	exit;

}
