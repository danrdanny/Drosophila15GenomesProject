#!/usr/bin/perl

# Description
#
# Please use -h to see a full description.
# Stocks to align/analyze are listed in $pwd/sampleSheet.tsv.
# Lines beginning with '#' are ignored.

use strict;
use Getopt::Std;

## Command-line options
my %opts;
getopts('t:hw:', \%opts); # Values in %opts

## Usage Statement
if ($opts{'h'} || !$opts{'t'} || !$opts{'w'}) {
	print "
	Required:
	t	threads
	w	read type: must be either allReads or pass

	Optional:
	-none
	\n";
	exit 0;
}

die "wrong read type" unless $opts{'w'} eq "pass" || $opts{'w'} eq "allReads";

## Get our current working directory
my $pwd = `pwd`; chomp($pwd);
## Threads to run multithreadded programs with
my $threads = $opts{'t'};
## Which reads
my $readType = "allReads";
   $readType = "pass" if $opts{'w'} eq "pass";

## Two subroutines are used, one to run commands, the second to log all activities
sub executeCommand {
	open LOGF, ">>$_[1]";
	print LOGF "[".localtime()."] CMD: $_[0]\n";
	close LOGF;
	my $output = `$_[0]`;
	return($output);
}
## End Subroutines

## Decide what stocks we're working on
my(%stocksToAlign,$numStocksToAlign);
if ($opts{'f'} =~ /[a-z]/i) {
	$stocksToAlign{$opts{'n'}} = $opts{'f'};
} else {
	# Open sampleSheet.tsv and work off of that
	my $sampleSheet = "sampleSheet.tsv";
	open INF,"$pwd/$sampleSheet" or die "Can't open $sampleSheet: $!";
	while (<INF>) {
        	chomp($_);
        	next if $_ =~ /^#/; # skip any line that starts with #
        	next unless $_ =~ /[a-z0-9]/i;
        	my(@F) = split /\t/, $_;
#print "$_\n";
        	my $sampleName  = $F[0];
        	my $barcode     = $F[1];
        	my $LIMS_Order  = $F[2];
		my $flowcell	= $F[3];

		my $file1 = "$LIMS_Order/$flowcell/n_1_1_$barcode.fastq.gz";
		my $file2 = "$LIMS_Order/$flowcell/n_1_2_$barcode.fastq.gz";

        	print "[".localtime(time)."] $sampleName | Checking for files in $LIMS_Order ... ";
               	print "$file1 missing - exiting\n" if !-e "$file1";
               	die if !-e "$file1";
               	print "$file2 missing - exiting\n" if !-e "$file1";
               	die if !-e "$file2";

        	print " OK\n";

        	$stocksToAlign{$sampleName} = $_;
        	++$numStocksToAlign;
	}
	close INF;
}

## Loop over the stocks we should be aligning.
foreach my $newStockName (sort keys %stocksToAlign) {
	print "[".localtime(time)."] $newStockName | ----- Start -----\n";
	my $targetDir = "$pwd/polish_all_racon";
	   $targetDir = "$pwd/polish_pass_racon" if $readType eq "pass";
	my $logFile = "$targetDir/$newStockName.log";
	#executeCommand("mkdir -p $targetDir",$logFile);

       	my(@F) = split /\t/, $stocksToAlign{$newStockName};

       	my $barcode     = $F[1];
       	my $LIMS_Order  = $F[2];
	my $flowcell	= $F[3];
	my $refGenome	= $F[4];

	my($f_reads,$r_reads);

	$f_reads  = "$pwd/$LIMS_Order/$flowcell/n_1_1_$barcode.fastq.gz ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_2_1_$barcode.fastq.gz ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_3_1_$barcode.fastq.gz ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_4_1_$barcode.fastq.gz ";
	$r_reads  = "$pwd/$LIMS_Order/$flowcell/n_1_2_$barcode.fastq.gz ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_2_2_$barcode.fastq.gz ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_3_2_$barcode.fastq.gz ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_4_2_$barcode.fastq.gz ";

	foreach my $raconCount (1..4) {
		my $lastRaconCount = $raconCount - 1;
		my $lastFasta = "$newStockName.$readType.minimap2.racon.x$lastRaconCount.fasta";
		   $lastFasta = "$newStockName.2.1.0.$readType.minimap2.assembly.fasta" if $lastRaconCount == 0;
		my $newFasta = "$newStockName.$readType.minimap2.racon.x$raconCount.fasta";

		executeCommand("minimap2 -t $threads $targetDir/$lastFasta /n/projects/dem/nanopore/rawData/$newStockName/$newStockName.2.1.0.$readType.fastq > $targetDir/$newStockName.mappings.x$raconCount.paf",$logFile);
		executeCommand("/home/mec/t/racon/bin/racon -t $threads /n/projects/dem/nanopore/rawData/$newStockName/$newStockName.2.1.0.$readType.fastq $targetDir/$newStockName.mappings.x$raconCount.paf $targetDir/$lastFasta $targetDir/$newFasta > $targetDir/racon.$raconCount.out",$logFile);
		executeCommand("script_assemblyStats.pl $targetDir/$newFasta > $targetDir/$newFasta.stats",$logFile);
	}

	executeCommand("rm -f $targetDir/*.paf",$logFile);

	my $newFasta = "$newStockName.$readType.minimap2.racon.x4.fasta";
	my $newBam = "$newStockName.$readType.minimap2.racon.x4.bam";
	# Align with BWA
	print "[".localtime(time)."] $newStockName | RUN: Aligning reads to $newFasta\n";
	executeCommand("bwa index $targetDir/$newFasta",$logFile);
	executeCommand("bwa mem -t $threads $targetDir/$newFasta '<zcat $f_reads' '<zcat $r_reads' > $targetDir/$newStockName.sam 2>/dev/null",$logFile);
	executeCommand("samtools view -bS $targetDir/$newStockName.sam > $targetDir/$newStockName.bam",$logFile);
	executeCommand("samtools sort -\@ $threads $targetDir/$newStockName.bam -o $targetDir/$newStockName.sorted",$logFile); #creates x.sorted.bam
	executeCommand("mv $targetDir/$newStockName.sorted $targetDir/$newBam",$logFile);
	executeCommand("rm -f $targetDir/$newStockName.bam $targetDir/*.sam $targetDir/*.amb $targetDir/*.ann $targetDir/*.bwt $targetDir/*.fai $targetDir/*.pac $targetDir/*.sa",$logFile);
	executeCommand("samtools index $targetDir/$newBam",$logFile);
	print "[".localtime(time)."] $newStockName | RUN: $newBam alignment complete\n";

	# Call snps using samtools
	my $pid = fork();
	die "unable to fork: $!" unless defined($pid);
	print "[".localtime(time)."] $newStockName | RUN: Calling SNPs for $newBam\n";
	if (!$pid) {
		executeCommand("samtools mpileup -ugf $targetDir/$newFasta $targetDir/$newBam | bcftools call -vmO z -o $targetDir/$newStockName.$readType.minimap2.racon.x4.vcf.gz",$logFile);
		executeCommand("rm -f $targetDir/$newBam $targetDir/$newBam.bai",$logFile);
		print "[".localtime(time)."] $newStockName | RUN: Called SNPs for $newBam\n";
               	exit 1;
	}

	print "[".localtime(time)."] $newStockName | ----- Done -----\n";
	print "[".localtime(time)."]\n";
}

print "[".localtime(time)."] Goodbye.\n";
