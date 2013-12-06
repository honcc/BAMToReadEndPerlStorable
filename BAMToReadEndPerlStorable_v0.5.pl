#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle);
use threads;
use threads::shared;
use Statistics::Descriptive;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to read a BAM file and generate the read-end pileup in format of perl storables. The method is largely based on BAMPolyAFilterer_v0.1
#
#	Input
#		--BAMPath=				file path [compulsory]; the path of the sorted bam file;
#		--fastaPath=			file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--IGVGenomePath=		file path [compulsory]; the path of IGV .genome file;
#		--baseComposition=		'yes' or 'no' ['no']; will generate base composition breakdown of reads at each position instead of total coverage; if no, the total coverage will be output;
#		--countMode=			full or 5 or 3 or midPt [5]; to record the 5'end or 3'end, full read or mid point of the read [mid point defined at the middle of 5’end and 3’end]; midPt is invalid if baseComposition is 'yes', will be reset to 5
#		--offset=				integer [0]; to record the position of $offset nt downstream (if +ve) or upstream (if -ve) of the position, works in all countModes;
#		--maxThread=			integer [4]; max number of threads to be used;
#		--outDir=				directory path ['./BAMToReadEndPerlStorable/']; output directory;
#
#	v0.2
#		[Wed 17 Jul 2013 13:57:57 CEST] added the offset option, to record a certain offset distance downstream/upstream of the read end;
#
#	v0.3
#		[Sat 20 Jul 2013 11:57:47 CEST] multithread abilities added, using the use Thread module, maxThread option added; reference: http://perldoc.perl.org/perlthrtut.html
#
#	v0.4
#		[Thu  1 Aug 2013 15:39:47 CEST] added --baseComposition= option to generate base composition breakdown of reads at each position instead of total coverage;
#		[Fri  2 Aug 2013 10:30:56 CEST] --5Or3= option is changed to --countMode=
#
#	v0.5
#		[Sat  7 Sep 2013 14:33:34 CEST] --countMode= can be 'midPt'
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-27 15:12]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/BAMToReadEndPerlStorable/v0.5/BAMToReadEndPerlStorable_v0.5.pl --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta --IGVGenomePath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927_v4.2.genome --baseComposition=no --countMode=5 --offset=0 --maxThread=10 --BAMPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/bowtie/pooled.bt2.local.sorted.bam --outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/BAMToReadEndPerlStorable/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/BAMToReadEndPerlStorable/v0.5/BAMToReadEndPerlStorable_v0.5.pl
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/plasmoDB/PF3D7/Pf3D7_01_v3_JAN2012_withMitoPlstd.fa
#	--IGVGenomePath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/plasmoDB/PF3D7/Pf3D7_01_v3_JAN2012_withMitoPlstd.genome
#	--baseComposition=no
#	--countMode=midPt
#	--offset=0
#	--maxThread=4
#	--BAMPath=/Volumes/C_Analysis/NGS/results/plasmodium/gDNARemap/gSimulation/bowtie/mapped.sorted.bam
#	--outDir=/Volumes/C_Analysis/NGS/results/plasmodium/gDNARemap/gSimulation/BAMToReadEndPerlStorable/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: checkIGVtoolsVersion|206, checkSamtoolsVersion|254, printCMDLogOrFinishMessage|362, readParameters|840
#	secondaryDependOnSub: currentTime|310, reportStatus|890
#
#<section ID="startingTasks" num="0">
########################################################################## 
&printCMDLogOrFinishMessage("CMDLog");#->362
my ($BAMPath, $fastaPath, $IGVGenomePath, $countMode, $offset, $maxThread, $baseComposition, $outDir) = &readParameters();#->840

&checkSamtoolsVersion();#->254
&checkIGVtoolsVersion();#->206
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $paramTag = "countMode.$countMode.offset.$offset.baseComp.$baseComposition";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
push @mkDirAry, $outDir;
my $readEndDir = "$outDir/$paramTag/"; push @mkDirAry, $readEndDir;
my $cntgCovPlsDir = "$readEndDir/cntgCovPls/"; push @mkDirAry, $cntgCovPlsDir;
my $wigDir = "$readEndDir/wiggle/"; push @mkDirAry, $wigDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
my $wigPrefix = "$wigDir/$paramTag";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: createEmptyGenomeCovPerlStorable|274, getIndivCntgCovPlsPath|328, readBAMAllCntgMultiThread|524, readMultiFasta|786
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|226, readBAMIndividualCntgBaseComposition|579, readBAMIndividualCntgCov|689, reportStatus|890
#
#<section ID="processInputData" num="4">
########################################################################## 
#----------Read Fasta
my $fastaHsh_ref = &readMultiFasta($fastaPath);#->786

#----------Create empty stroable
my $cntgCovIdxHshPath = &createEmptyGenomeCovPerlStorable($cntgCovPlsDir, $fastaHsh_ref);#->274
my $covPlsPathHsh_ref = &getIndivCntgCovPlsPath($cntgCovIdxHshPath);#->328

$fastaHsh_ref = {};
&readBAMAllCntgMultiThread($BAMPath, $countMode, $covPlsPathHsh_ref, $offset, $maxThread, $baseComposition);#->524
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_outputData
#	primaryDependOnSub: printWigFromBaseComPerlStorable|395, printWigFromCovPerlStorable|464
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|226, reportStatus|890
#
#<section ID="outputData" num="5">
if ($baseComposition eq 'no') {
	&printWigFromCovPerlStorable($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath);#->464
} elsif ($baseComposition eq 'yes') {
	&printWigFromBaseComPerlStorable($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath);#->395
}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|362
#	secondaryDependOnSub: currentTime|310
#
#<section ID="finishingTasks" num="6">
&printCMDLogOrFinishMessage("finishMessage");#->362
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	checkTools [n=2]:
#		checkIGVtoolsVersion, checkSamtoolsVersion
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=7]:
#		checkIGVtoolsVersion, checkSamtoolsVersion, currentTime
#		printCMDLogOrFinishMessage, readMultiFasta, readParameters
#		reportStatus
#
#	multithread [n=1]:
#		checkRunningThreadAndWaitToJoin
#
#	reporting [n=1]:
#		currentTime
#
#	storable [n=2]:
#		createEmptyGenomeCovPerlStorable, getIndivCntgCovPlsPath
#
#	unassigned [n=5]:
#		printWigFromBaseComPerlStorable, printWigFromCovPerlStorable, readBAMAllCntgMultiThread
#		readBAMIndividualCntgBaseComposition, readBAMIndividualCntgCov
#
#====================================================================================================================================================#

sub checkIGVtoolsVersion {
#....................................................................................................................................................#
#	subroutineCategory: general, checkTools
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80
#	secondaryAppearInSection: >none
#	input: none
#	output: none
#	toCall: &checkIGVtoolsVersion();
#	calledInLine: 90
#....................................................................................................................................................#
	
	my $IGVtoolsStdout = `igvtools 2>&1`;
	if ($IGVtoolsStdout =~ m/\s+(Version \S+)\s+/) {
		&reportStatus("Checking: IGVtools: $1", 0, "\n");#->890
	} else {
		die "IGVtools not installed properly. Quitting.\n";
	}
}
sub checkRunningThreadAndWaitToJoin {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: reportStatus|890
#	appearInSub: printWigFromCovPerlStorable|464, readBAMAllCntgMultiThread|524
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|130, 5_outputData|149
#	input: $sleepTime, $verbose
#	output: none
#	toCall: &checkRunningThreadAndWaitToJoin($verbose, $sleepTime);
#	calledInLine: 520, 575
#....................................................................................................................................................#
	
	my ($verbose, $sleepTime) = @_;
	
	my @runningThrAry = threads->list(threads::running);
	my @joinableThrAry = threads->list(threads::joinable);
	while (@runningThrAry or @joinableThrAry) {
		@runningThrAry = threads->list(threads::running);
		@joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
		my $numThreadRunning = scalar @runningThrAry;
		&reportStatus("The last $numThreadRunning threads are still running", 20, "\r") if $verbose eq 'yes';#->890
		sleep $sleepTime;
	}
}
sub checkSamtoolsVersion {
#....................................................................................................................................................#
#	subroutineCategory: general, checkTools
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80
#	secondaryAppearInSection: >none
#	input: none
#	output: none
#	toCall: &checkSamtoolsVersion();
#	calledInLine: 89
#....................................................................................................................................................#

	my $samtoolsStdout = `samtools 2>&1`;
	if ($samtoolsStdout =~ m/\s+(Version: \S+)\s+/) {
		&reportStatus("Checking: samtools: $1", 0, "\n");#->890
	} else {
		die "samtools not installed properly. Quitting.\n";
	}
}
sub createEmptyGenomeCovPerlStorable {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|130
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsDir, $fastaHsh_ref
#	output: $cntgCovIdxHshPath
#	toCall: my ($cntgCovIdxHshPath) = &createEmptyGenomeCovPerlStorable($cntgCovPlsDir, $fastaHsh_ref);
#	calledInLine: 140
#....................................................................................................................................................#

	my ($cntgCovPlsDir, $fastaHsh_ref) = @_;
	
	my $cntgCovPlsIdxHsh_ref = {};
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		
		&reportStatus("Creating empty storabe for $cntg", 20,"\r");#->890
		
		my $cntgLen = length($fastaHsh_ref->{$cntg});
		my $cntgCovAry_ref = ();
		foreach (1..$cntgLen) {
			push @{$cntgCovAry_ref}, undef;
		}
		my $cntgCovPlsName = "$cntg.ary.pls";
		my $cntgCovPlsPath = "$cntgCovPlsDir/$cntgCovPlsName";
		$cntgCovPlsIdxHsh_ref->{$cntg} = $cntgCovPlsName;
		store($cntgCovAry_ref, "$cntgCovPlsPath");
	}

	my $cntgCovIdxHshPath = "$cntgCovPlsDir/index.hsh.pls";
	store($cntgCovPlsIdxHsh_ref, "$cntgCovIdxHshPath");
	
	return $cntgCovIdxHshPath;
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|362, reportStatus|890
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|80, 6_finishingTasks|163
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 382, 385, 390, 906
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|130
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 141, 340
#....................................................................................................................................................#
	
	#my $cntgCovInPlsPathHsh_ref = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);#->328
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	&reportStatus("pls path of $numCntg contig stored", 0, "\n");#->890
	
	return $cntgCovInPlsPathHsh_ref;
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|310
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80, 6_finishingTasks|163
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 86, 168
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->310
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->310
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->310
		print "=========================================================================\n\n";
	}
}
sub printWigFromBaseComPerlStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 5_outputData|149
#	secondaryAppearInSection: >none
#	input: $IGVGenomePath, $covPlsPathHsh_ref, $wigPrefix
#	output: none
#	toCall: &printWigFromBaseComPerlStorable($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath);
#	calledInLine: 157
#....................................................................................................................................................#

	my ($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath) = @_;

	my $wigInfoHsh_ref = {};
	
	foreach my $plusOrMinus (qw/plus minus/) {
		foreach my $ATGC (qw/A T G C/) {
			my $path = "$wigPrefix.$plusOrMinus.$ATGC.wig.gz";
			$wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'path'} = $path;
			open $wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'FH'}, "| gzip -fc >$path";
		}
	}
	
	foreach my $cntg (sort {$a cmp $b} keys %{$covPlsPathHsh_ref}) {
		my $cntgCovPlsPath = $covPlsPathHsh_ref->{$cntg};
		&reportStatus("Writing $cntg", 20, "\r");#->890

		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);

		foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
			foreach my $ATGC (keys %{$wigInfoHsh_ref->{$plusOrMinus}}) {
				print {$wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'FH'}} "variableStep chrom=$cntg span=1\n";
			}
		}

		foreach my $index (0..$#{$cntgCovAry_ref}) {
			if (defined $cntgCovAry_ref->[$index]) {
				my $pos = $index+1;
				my %strandCountStringHsh = {};
				($strandCountStringHsh{'plus'}, $strandCountStringHsh{'minus'}) = split /,/, $cntgCovAry_ref->[$index];
				foreach my $plusOrMinus (keys %strandCountStringHsh) {
					my %baseCountHsh = {};
					($baseCountHsh{'A'}, $baseCountHsh{'T'}, $baseCountHsh{'G'}, $baseCountHsh{'C'}) = split /:/, $strandCountStringHsh{$plusOrMinus};
					foreach my $ATGC (keys %baseCountHsh) {
						my $count = $baseCountHsh{$ATGC};
						print {$wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'FH'}} join '', ((join "\t", ($pos, $count)), "\n") if ($count > 0);
					}
				}
			}
		}
	}
	
	foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
		foreach my $ATGC (keys %{$wigInfoHsh_ref->{$plusOrMinus}}) {
			close $wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'FH'};
		}
	}
	
	foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
		foreach my $ATGC (keys %{$wigInfoHsh_ref->{$plusOrMinus}}) {
			my $path = $wigInfoHsh_ref->{$plusOrMinus}{$ATGC}{'path'};
			my $cmd = "igvtools toTDF $path $wigPrefix.$plusOrMinus.$ATGC.tdf $IGVGenomePath 2&>/dev/null;";
			&reportStatus("Converting $path", 0, "\n");#->890
			system $cmd;
		}
	}
}
sub printWigFromCovPerlStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|226, reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 5_outputData|149
#	secondaryAppearInSection: >none
#	input: $IGVGenomePath, $covPlsPathHsh_ref, $wigPrefix
#	output: none
#	toCall: &printWigFromCovPerlStorable($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath);
#	calledInLine: 155
#....................................................................................................................................................#

	my ($covPlsPathHsh_ref, $wigPrefix, $IGVGenomePath) = @_;

	my $wigInfoHsh_ref = {};
	
	foreach my $plusOrMinus (qw/plus minus/) {
		my $path = "$wigPrefix.$plusOrMinus.wig.gz";
		$wigInfoHsh_ref->{$plusOrMinus}{'path'} = $path;
		open $wigInfoHsh_ref->{$plusOrMinus}{'FH'}, "| gzip -fc >$path";
	}
	
	foreach my $cntg (sort {$a cmp $b} keys %{$covPlsPathHsh_ref}) {
		my $cntgCovPlsPath = $covPlsPathHsh_ref->{$cntg};
		&reportStatus("Writing $cntg", 20,"\r");#->890

		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);

		foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
			print {$wigInfoHsh_ref->{$plusOrMinus}{'FH'}} "variableStep chrom=$cntg span=1\n";
		}

		foreach my $index (0..$#{$cntgCovAry_ref}) {
			if (defined $cntgCovAry_ref->[$index]) {
				my $pos = $index+1;
				my %strandCountHsh = ();
				($strandCountHsh{'plus'}, $strandCountHsh{'minus'}) = split /,/, $cntgCovAry_ref->[$index];
				foreach my $plusOrMinus (keys %strandCountHsh) {
					my $count = $strandCountHsh{$plusOrMinus};
					print {$wigInfoHsh_ref->{$plusOrMinus}{'FH'}} join '', ((join "\t", ($pos, $count)), "\n") if ($count > 0);
				}
			}
		}
	}
	
	foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
		close $wigInfoHsh_ref->{$plusOrMinus}{'FH'};
	}
	
	foreach my $plusOrMinus (keys %{$wigInfoHsh_ref}) {
		my $path = $wigInfoHsh_ref->{$plusOrMinus}{'path'};
		my $cmd = "igvtools toTDF $path $wigPrefix.$plusOrMinus.tdf $IGVGenomePath 2&>/dev/null;";
		&reportStatus("Issue a thread to convert $plusOrMinus wig", 0, "\n");#->890
		threads->create(sub{system $cmd;});
	}
	
	&checkRunningThreadAndWaitToJoin('yes', 1);#->226

}
sub readBAMAllCntgMultiThread {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|226, readBAMIndividualCntgBaseComposition|579, readBAMIndividualCntgCov|689, reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|130
#	secondaryAppearInSection: >none
#	input: $BAMPath, $baseComposition, $countMode, $covPlsPathHsh_ref, $maxThread, $offset
#	output: none
#	toCall: &readBAMAllCntgMultiThread($BAMPath, $countMode, $covPlsPathHsh_ref, $offset, $maxThread, $baseComposition);
#	calledInLine: 144
#....................................................................................................................................................#
	
	my ($BAMPath, $countMode, $covPlsPathHsh_ref, $offset, $maxThread, $baseComposition) = @_;
	
	my $cntgProc = 0;
	my $cntgNum = keys %{$covPlsPathHsh_ref};
	my $sleepTime = 1/$maxThread;
	
	&reportStatus("Start multi-threads counting of BAM file", 20, "\n");#->890

	foreach my $cntg (keys %{$covPlsPathHsh_ref}) {
		my $treadIssued = 'no';
		$cntgProc++;
		my $cntgCovPlsPath = $covPlsPathHsh_ref->{$cntg};

		while ($treadIssued eq 'no') {
			#---check num of running threads
			my @runningThrAry = threads->list(threads::running);
			#----if running threads lower than the limit
			if (@runningThrAry < $maxThread) {
				#---spawn a new thread
				if ($baseComposition eq 'yes') {
					threads->create(\&readBAMIndividualCntgBaseComposition, ($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg));#->579
				} else {
					threads->create(\&readBAMIndividualCntgCov, ($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg));#->689
				}
				$treadIssued = 'yes';
				@runningThrAry = threads->list(threads::running);
				my $threadNum = scalar @runningThrAry;
				&reportStatus("$cntgProc cntg counted with $threadNum threads running", 20, "\r");#->890
			}
			sleep $sleepTime;
		}

		#----detach the finished threads
		my @joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
	}
	
	&checkRunningThreadAndWaitToJoin('yes', 1);#->226

}
sub readBAMIndividualCntgBaseComposition {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readBAMAllCntgMultiThread|524
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|130
#	input: $BAMPath, $cntg, $cntgCovPlsPath, $countMode, $offset
#	output: none
#	toCall: &readBAMIndividualCntgBaseComposition($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg);
#	calledInLine: 556
#....................................................................................................................................................#

	my ($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg) = @_;
	
	open BAMIN, "samtools view $BAMPath \'$cntg\' |"; #---no header

	my $cntgPosHsh_ref = {};
	my $readProc = 0;
	while (my $theLine = <BAMIN>) { 

		$readProc++;
		my (undef, $flag, undef, $readStart, undef, $cigarStr, undef, undef, undef, $readSeq, undef) = split /\t/, $theLine;

		#--- skip the I reads ad hoc for simpler code
		next if $cigarStr =~ m/(\d+)I/;
		
		my $genomicLength = 0;
		while ($cigarStr =~ /(\d+)([M|N|D])/g) {$genomicLength += $1;}
	
		my $readEndPosHsh_ref = {};
		my $readStrand;
		my $readPosBaseHsh_ref = ();
		my $offsetOnStrnd = $offset;
		
		if ($flag & 16) {
			$readStrand = "-";
			$readEndPosHsh_ref->{3} = $readStart;
			$readEndPosHsh_ref->{5} = $readStart + $genomicLength - 1;
			$readSeq = reverse $readSeq;
			$readSeq =~ tr/ACGTacgt/TGCAtgca/;
			$offsetOnStrnd = -1*$offset;
			
		} else {
			$readStrand = "+";
			$readEndPosHsh_ref->{5} = $readStart;
			$readEndPosHsh_ref->{3} = $readStart + $genomicLength - 1;
		}
		
		if ($countMode eq '5') {
			
			my $base = substr $readSeq, $offset, 1;
			my $offSetPos = $readEndPosHsh_ref->{5} + $offsetOnStrnd;
			$readPosBaseHsh_ref->{$offSetPos} = $base;
			
		} elsif ($countMode eq '3') {

			my $base = substr $readSeq, -1*$offset, 1;
			my $offSetPos = $readEndPosHsh_ref->{3} + $offsetOnStrnd;
			$readPosBaseHsh_ref->{$offSetPos} = $base;

		} elsif ($countMode eq 'full') {
			
			#---get all aligned position
			my @alignedPosAry = ();
			my @cigarSegmentAry = ($cigarStr =~ /(\d+M|\d+N|\d+D)/g); #--refer to the definition of "The Sequence Alignment/Map format and SAMtools" BIOINFORMATICS APPLICATIONS NOTE
			my $curntStartPos = $readStart;
			foreach my $cigarSegment (@cigarSegmentAry) {
				my ($numPos, $CIGARType) = ($cigarSegment =~ /(\d+)(\w)/);
				my $curntEndPos = $curntStartPos + $numPos - 1;
				push @alignedPosAry, ($curntStartPos..$curntEndPos) if $CIGARType eq 'M';
				$curntStartPos = $curntEndPos + 1;
			}

			@alignedPosAry = reverse @alignedPosAry if $readStrand eq '-';
			foreach my $i (0..$#alignedPosAry) {
				my $offSetPos = $alignedPosAry[$i] + $offsetOnStrnd;
				my $base = substr $readSeq, $i, 1;
				$readPosBaseHsh_ref->{$offSetPos} = $base;
			}
		}
		
		#---store the position and bases
		foreach my $offSetPos (keys %{$readPosBaseHsh_ref}) {
			my $base = $readPosBaseHsh_ref->{$offSetPos};
			$cntgPosHsh_ref->{$offSetPos}{$readStrand}{$base}++;
		}
	}
	close BAMIN;

	#---if there are reads
	if ($readProc > 0) {
		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);

		foreach my $offSetPos (keys %{$cntgPosHsh_ref}) {
			my $i = $offSetPos - 1;
			next if ($i > $#{$cntgCovAry_ref} or $i < 0);
			foreach my $base ("A", "T", "G", "C") {
				foreach my $readStrand ("+", "-") {
					$cntgPosHsh_ref->{$offSetPos}{$readStrand}{$base} = 0 if not $cntgPosHsh_ref->{$offSetPos}{$readStrand}{$base};
				}
			}
			my $plusCovStr = join ":", ($cntgPosHsh_ref->{$offSetPos}{'+'}{'A'}, $cntgPosHsh_ref->{$offSetPos}{'+'}{'T'}, $cntgPosHsh_ref->{$offSetPos}{'+'}{'G'}, $cntgPosHsh_ref->{$offSetPos}{'+'}{'C'});
			my $minusCovStr = join ":", ($cntgPosHsh_ref->{$offSetPos}{'-'}{'A'}, $cntgPosHsh_ref->{$offSetPos}{'-'}{'T'}, $cntgPosHsh_ref->{$offSetPos}{'-'}{'G'}, $cntgPosHsh_ref->{$offSetPos}{'-'}{'C'});
			$cntgCovAry_ref->[$i] = join ',', ($plusCovStr, $minusCovStr);
			my $pos = $i+1;
		}
		store($cntgCovAry_ref, "$cntgCovPlsPath");
	}
}
sub readBAMIndividualCntgCov {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readBAMAllCntgMultiThread|524
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|130
#	input: $BAMPath, $cntg, $cntgCovPlsPath, $countMode, $offset
#	output: none
#	toCall: &readBAMIndividualCntgCov($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg);
#	calledInLine: 558
#....................................................................................................................................................#

	my ($BAMPath, $countMode, $cntgCovPlsPath, $offset, $cntg) = @_;
	
	open BAMIN, "samtools view $BAMPath \'$cntg\' |"; #---no header

	my $cntgPosHsh_ref = {};
	my $readProc = 0;
	while (my $theLine = <BAMIN>) { 

		$readProc++;
		my (undef, $flag, undef, $readStart, undef, $cigarStr, undef, undef, undef, undef) = split /\t/, $theLine;

		#--- skip the I reads ad hoc for simpler code
		next if $cigarStr =~ m/(\d+)I/;

		my $genomicLength = 0;
		while ($cigarStr =~ /(\d+)([M|N|D])/g) {$genomicLength += $1;}
	
		my $readEndPosHsh_ref = {};
		my $readStrand;
		my $offsetOnStrnd = $offset;
		
		if ($flag & 16) {
			$readStrand = "-";
			$readEndPosHsh_ref->{3} = $readStart;
			$readEndPosHsh_ref->{5} = $readStart + $genomicLength - 1;
			$readEndPosHsh_ref->{midPt} = $readEndPosHsh_ref->{3}+ int(($readEndPosHsh_ref->{5}-$readEndPosHsh_ref->{3})/2);
			$offsetOnStrnd = -1*$offset;
			
		} else {
			$readStrand = "+";
			$readEndPosHsh_ref->{5} = $readStart;
			$readEndPosHsh_ref->{3} = $readStart + $genomicLength - 1;
			$readEndPosHsh_ref->{midPt} = $readEndPosHsh_ref->{5}+ int(($readEndPosHsh_ref->{3}-$readEndPosHsh_ref->{5})/2);
		}
		
		my @posAry = ();
		
		if ($countMode eq '5') {
			@posAry = ($readEndPosHsh_ref->{5});

		} elsif ($countMode eq '3') {
			@posAry = ($readEndPosHsh_ref->{3});

		} elsif ($countMode eq 'midPt') {
			@posAry = ($readEndPosHsh_ref->{midPt});

		} elsif ($countMode eq 'full') {
			
			#---get all aligned position
			#$cigarStr =~ s/\d+S//g;
			my @cigarSegmentAry = ($cigarStr =~ /(\d+M|\d+N|\d+D)/g); #--refer to the definition of "The Sequence Alignment/Map format and SAMtools" BIOINFORMATICS APPLICATIONS NOTE
			my $curntStartPos = $readStart;
			foreach my $cigarSegment (@cigarSegmentAry) {
				my ($numPos, $CIGARType) = ($cigarSegment =~ /(\d+)(\w)/);
				my $curntEndPos = $curntStartPos + $numPos - 1;
				push @posAry, ($curntStartPos..$curntEndPos) if $CIGARType eq 'M';
				$curntStartPos = $curntEndPos + 1;
			}
		}

		#---store the position and bases
		foreach my $pos (@posAry) {
			my $offSetPos = $pos + $offsetOnStrnd;
			$cntgPosHsh_ref->{$offSetPos}{$readStrand}++;
		}
	}
	close BAMIN;

	#---if there are reads
	if ($readProc > 0) {
		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);

		foreach my $offSetPos (keys %{$cntgPosHsh_ref}) {
			my $i = $offSetPos - 1;
			next if ($i > $#{$cntgCovAry_ref} or $i < 0);
			my $plusCov = my $minusCov = 0;
			#----at least one strand has data
			$plusCov = $cntgPosHsh_ref->{$offSetPos}{'+'} if $cntgPosHsh_ref->{$offSetPos}{'+'};
			$minusCov = $cntgPosHsh_ref->{$offSetPos}{'-'} if $cntgPosHsh_ref->{$offSetPos}{'-'};
			$cntgCovAry_ref->[$i] = join ',', ($plusCov, $minusCov);
		}
		store($cntgCovAry_ref, "$cntgCovPlsPath");
	}
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|130
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 137
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->890
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/ +/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: reportStatus|890
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80
#	secondaryAppearInSection: >none
#	input: none
#	output: $BAMPath, $IGVGenomePath, $baseComposition, $countMode, $fastaPath, $maxThread, $offset, $outDir
#	toCall: my ($BAMPath, $fastaPath, $IGVGenomePath, $countMode, $offset, $maxThread, $baseComposition, $outDir) = &readParameters();
#	calledInLine: 87
#....................................................................................................................................................#
	
	my ($BAMPath, $fastaPath, $IGVGenomePath, $countMode, $offset, $maxThread, $baseComposition, $outDir);
	
	$countMode = 5;
	$offset = 0;
	$maxThread = 4;
	my $dirPath = dirname(rel2abs($0));
	$outDir ="$dirPath/BAMToReadEndPerlStorable/";
	$baseComposition = 'no';
	
	GetOptions 	("BAMPath=s" => \$BAMPath,
				 "fastaPath=s"  => \$fastaPath,
				 "IGVGenomePath=s"  => \$IGVGenomePath,
				 "countMode:s"  => \$countMode,
				 "maxThread:i"  => \$maxThread,
				 "offset:i"  => \$offset,
				 "baseComposition:s"  => \$baseComposition,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($BAMPath, $fastaPath, $IGVGenomePath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	#---for the offset and countMode values
	$offset = 0 if $countMode eq 'full';
	if ($countMode eq 'midPt' and $baseComposition eq 'yes') {
		$countMode = 5;
		&reportStatus("countMode reset to 5 since $baseComposition is 'yes'", 10, "\n");#->890
	}
	die 'countMode has to be 3 or 5 or full or midPt' if ($countMode ne '3' and $countMode ne '5' and $countMode ne 'full' and $countMode ne 'midPt');
	
	return($BAMPath, $fastaPath, $IGVGenomePath, $countMode, $offset, $maxThread, $baseComposition, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|310
#	appearInSub: checkIGVtoolsVersion|206, checkRunningThreadAndWaitToJoin|226, checkSamtoolsVersion|254, createEmptyGenomeCovPerlStorable|274, getIndivCntgCovPlsPath|328, printWigFromBaseComPerlStorable|395, printWigFromCovPerlStorable|464, readBAMAllCntgMultiThread|524, readMultiFasta|786, readParameters|840
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|80, 4_processInputData|130, 5_outputData|149
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 220, 249, 268, 291, 357, 421, 458, 488, 516, 542, 563, 804, 883
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->310

	return ();
}

exit;
