#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options
my %opts;
getopts('t:hw:', \%opts); # Values in %opts

## Usage Statement
if ($opts{'h'} || !$opts{'t'}) {
	print "
	\n";
	exit 0;
}

## Get our current working directory
my $pwd = `pwd`; chomp($pwd);
## Threads to run multithreadded programs with
my $threads = $opts{'t'};
my $readType = "allReads";
   $readType = "pass" if $opts{'w'} eq "pass";
my $maxRuns = 6;

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

		my $file1 = "$LIMS_Order/$flowcell/n_1_1_$barcode.fastq";
		my $file2 = "$LIMS_Order/$flowcell/n_1_2_$barcode.fastq";

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
	my $targetDir = "$pwd/polish_all_pilon";
	   $targetDir = "$pwd/polish_pass_pilon" if $readType eq "pass";
	my $logFile = "$targetDir/$newStockName.log";

	my $finalFasta = "$newStockName.$readType.minimap2.pilon.x$maxRuns.fasta";
	my $finalBam   = "$newStockName.$readType.minimap2.pilon.x$maxRuns.bam";
	my $finalVCF   = "$newStockName.$readType.minimap2.pilon.x$maxRuns.vcf.gz";

	print "[".localtime(time)."] $newStockName | Skipping because $finalVCF exists\n" if -e "$targetDir/$finalVCF";
	next if -e "$targetDir/$finalVCF";

       	my(@F) = split /\t/, $stocksToAlign{$newStockName};

       	my $barcode     = $F[1];
       	my $LIMS_Order  = $F[2];
	my $flowcell	= $F[3];

	my($f_reads,$r_reads);

	$f_reads  = "$pwd/$LIMS_Order/$flowcell/n_1_1_$barcode.fastq ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_2_1_$barcode.fastq ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_3_1_$barcode.fastq ";
	$f_reads .= "$pwd/$LIMS_Order/$flowcell/n_4_1_$barcode.fastq ";
	$r_reads  = "$pwd/$LIMS_Order/$flowcell/n_1_2_$barcode.fastq ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_2_2_$barcode.fastq ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_3_2_$barcode.fastq ";
	$r_reads .= "$pwd/$LIMS_Order/$flowcell/n_4_2_$barcode.fastq ";

	foreach my $pilonRunNum (1..$maxRuns) {
		my $lastPilonRunNum = $pilonRunNum - 1;

		my $lastBam       = "$newStockName.$readType.minimap2.pilon.x$lastPilonRunNum.bam";
		my $lastFasta     = "$newStockName.$readType.minimap2.pilon.x$lastPilonRunNum.fasta";
	     	   $lastFasta     = "$newStockName.2.1.0.$readType.minimap2.assembly.fasta" if $pilonRunNum == 1;
		my $newBamPrefix  = "$newStockName.$readType.minimap2.pilon.x$pilonRunNum";
		my $newBam        = "$newStockName.$readType.minimap2.pilon.x$pilonRunNum.bam";
		my $newFasta	  = "$newStockName.$readType.minimap2.pilon.x$pilonRunNum.fasta";
		my $VCF 	  = "$newStockName.$readType.minimap2.pilon.x$lastPilonRunNum.vcf.gz";

		die "can't find $lastFasta" if !-e "$targetDir/$lastFasta";

		print "[".localtime(time)."] $newStockName | $newFasta exists, skipping iteration $pilonRunNum\n" if -e "$targetDir/$newFasta";
		next if -e "$targetDir/$newFasta";
		print "[".localtime(time)."] $newStockName | Iteration $pilonRunNum\n";

		if (!-e "$targetDir/$lastBam") {
			# Align with BWA
			print "[".localtime(time)."] $newStockName | Indexing $lastFasta.\n";
			executeCommand("bwa index $targetDir/$lastFasta 2>/dev/null",$logFile);
			print "[".localtime(time)."] $newStockName | Aligning reads to $lastFasta.\n";
			executeCommand("bwa mem -t $threads $targetDir/$lastFasta '<cat $f_reads' '<cat $r_reads' > $targetDir/$newStockName.sam 2>/dev/null",$logFile);
			print "[".localtime(time)."] $newStockName | Samtools view.\n";
			executeCommand("samtools view -bS $targetDir/$newStockName.sam > $targetDir/$newStockName.bam 2>/dev/null",$logFile);
			print "[".localtime(time)."] $newStockName | Sorting $newStockName.bam.\n";
			executeCommand("samtools sort -\@ $threads $targetDir/$newStockName.bam -o $targetDir/$newStockName.sorted 2>/dev/null",$logFile); #creates x.sorted.bam
			executeCommand("mv $targetDir/$newStockName.sorted $targetDir/$lastBam",$logFile);
			executeCommand("rm -f $targetDir/$newStockName.bai $targetDir/$newStockName.bam $targetDir/$newStockName.sam",$logFile);
			executeCommand("samtools index $targetDir/$lastBam 2>/dev/null",$logFile);
			print "[".localtime(time)."] $newStockName | Alignment complete.\n";
		} else {
			print "[".localtime(time)."] $newStockName | $lastBam exists, skipping alignment.\n";
		}

		if (!-e "$targetDir/$newFasta") {
			print "[".localtime(time)."] $newStockName | Running pilon on $lastBam.\n";
			executeCommand("java -Xmx70G -jar /home/dem/bin/pilon-1.22.jar --threads $threads --genome $targetDir/$lastFasta --frags $targetDir/$lastBam --output $targetDir/$newBamPrefix --changes --diploid",$logFile);
			executeCommand("script_assemblyStats.pl $targetDir/$newFasta > $targetDir/$newFasta.stats",$logFile);
			#executeCommand("rm -f $targetDir/*.amb $targetDir/*.ann $targetDir/*.bwt $targetDir/*.fai $targetDir/*.pac $targetDir/*.sa",$logFile);
		} else {
			print "[".localtime(time)."] $newStockName | $newFasta exists, skipping pilon step.\n";
		}

		if (!-e "$targetDir/$VCF") {
			print "[".localtime(time)."] $newStockName | Forking to make $VCF from $lastBam.\n";
			my $pid = fork();
			die "unable to fork: $!" unless defined($pid);
			if (!$pid) {
				print "[".localtime(time)."] $newStockName | Forked to $pid and calling SNPs on $lastBam\n";
				executeCommand("samtools mpileup -ugf $targetDir/$lastFasta $targetDir/$lastBam | bcftools call -vmO z -o $targetDir/$VCF 2>/dev/null",$logFile);
				#executeCommand("rm -f $targetDir/$lastBam $targetDir/$lastBam.bai",$logFile);
				print "[".localtime(time)."] $newStockName | $pid finished calling SNPs for $lastBam\n";
               			exit 1;
			}
		} else {
			print "[".localtime(time)."] $newStockName | $VCF exists, skipping SNP calling.\n";
		}
	}

	if (!-e "$targetDir/$finalBam") {
		print "[".localtime(time)."] $newStockName | Final alignment, indexing $finalFasta.\n";
		executeCommand("bwa index $targetDir/$finalFasta 2>/dev/null",$logFile);
		print "[".localtime(time)."] $newStockName | Final alignment, aligning reads to $finalFasta.\n";
		executeCommand("bwa mem -t $threads $targetDir/$finalFasta '<cat $f_reads' '<cat $r_reads' > $targetDir/$newStockName.sam 2>/dev/null",$logFile);
		print "[".localtime(time)."] $newStockName | Final alignment, samtools view.\n";
		executeCommand("samtools view -bS $targetDir/$newStockName.sam > $targetDir/$newStockName.bam 2>/dev/null",$logFile);
		print "[".localtime(time)."] $newStockName | Final alignment, Sorting reads.\n";
		executeCommand("samtools sort -\@ $threads $targetDir/$newStockName.bam -o $targetDir/$newStockName.sorted 2>/dev/null",$logFile);
		executeCommand("mv $targetDir/$newStockName.sorted $targetDir/$finalBam",$logFile);
		executeCommand("rm -f $targetDir/$newStockName.bai $targetDir/$newStockName.bam $targetDir/$newStockName.sam",$logFile);
		executeCommand("samtools index $targetDir/$finalBam 2>/dev/null",$logFile);
		#executeCommand("rm -f $targetDir/$newStockName*.amb $targetDir/$newStockName*.ann $targetDir/$newStockName*.bwt $targetDir/$newStockName*.fai $targetDir/$newStockName*.pac $targetDir/$newStockName*.sa",$logFile);
	} else {
		print "[".localtime(time)."] $newStockName | Final alignment, $finalBam exists, skipping alignment.\n";
	}

	if (!-e "$targetDir/$finalVCF") {
		my $pid = fork();
		die "unable to fork: $!" unless defined($pid);
		if (!$pid) {
			print "[".localtime(time)."] $newStockName | Forked to $pid and Calling SNPs on $finalBam\n";
			executeCommand("samtools mpileup -ugf $targetDir/$finalFasta $targetDir/$finalBam | bcftools call -vmO z -o $targetDir/$finalVCF 2>/dev/null",$logFile);
			print "[".localtime(time)."] $newStockName | PID: $pid finished calling SNPs for $finalBam\n";
			exit 1;
		}
	} else {
		print "[".localtime(time)."] $newStockName | Final alignment, $finalVCF exists, skipping SNP calling.\n";
	}

	print "[".localtime(time)."] $newStockName | ----- Done -----\n";
	print "[".localtime(time)."]\n";
}

print "[".localtime(time)."] Goodbye.\n";
