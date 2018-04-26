#!/usr/bin/perl

use strict;

use Getopt::Std; my%opt = ();
getopts("B:s:O:E:G:q:C:p:m:P:F:M:l:c:n:t:L:Q:vS:T:fe:a:b:d:D:k:K:Rx:X:AHr:Izj:u:U:J:N:o:V:wg:", \%opt);

my$die1 = "

================================================================= fragScaff =================================================================
                                       Version 140324.1 (180112 10X support update - beta -G X option)
                                                         Andrew Adey (adey\@ohsu.edu)

Attempts to scaffold contigs using fragment pool data. Performs best using CPT-Seq (Contiguity Preserving Transposase) libraries, but will
work using fosmid (decent) or LFR (Long Fragment Read) (poor) fragment data, though these methods are far less tested in fragScaff. It is
required that the input bam has the proper \@SQ lines, as well as \@RG lines and are aligned to the reference assembly provided as the -F
option. The bam (-B) must be reads aligned to the input assembled fasta (-F).

    fragScaff.pl -B <BAM/bamParse> -F <Input fasta> <Other options>
	
NOTE: It is recommended to run as 'fragScaff.pl -B <BAM> -b 1 -E <node size>' to first generate the bamParse file. This file can then
      be used for all subsequent stages to eliminate the need to re-parse the bam. Similarly it is recommended to run the node-node
      calculations once to generate a '.links.txt' file which can then be used as option '-K' for all subsequent runs thus eliminating the
      need to rerun the time-consuming node-node calculations.
	  
UPDATE: Date: Jan. 12, 2018. Updated to allow 10X genomics bam files by detecting the CR:Z field and building the RG info in pre-processing.
        Note: This has not been fully tested.
	  
RECOMMENDED   (1) 'fragScaff.pl -B <myBAM> -b 1 -E <node_size> (optional: -J <repeats.bed> -N <Nbases.bed> -o <max_node_size>)'
      RUNS:   (2) 'fragScaff.pl -B <myBAM.E#.bamParse> -A -O <myOUT> -t <myTHREAD#/Q>'
              (3) 'fragScaff.pl -B <myBAM.E#.bamParse> -K <myOUT.fragScaff.links.txt> -F <myIN.fa> -O <myOUT2> <other options>'

              This will allow you to run all pre-processing in steps (1) and (2) and then run the last stage (3) multiple times to compare
              outputs and find the optimum setting for your input. Though defaults should give an optimum output.

    -H   Display detailed options description.

";

my$die2 = "

================================================================= fragScaff =================================================================
                                       Version 140324.1 (180112 10X support update - beta -G X option)
                                                         Andrew Adey (adey\@ohsu.edu)

Attempts to scaffold contigs using fragment pool data. Performs best using CPT-Seq (Contiguity Preserving Transposase) libraries, but will
work using fosmid (decent) or LFR (Long Fragment Read) (poor) fragment data, though these methods are far less tested in fragScaff. It is
required that the input bam has the proper \@SQ lines, as well as \@RG lines and are aligned to the reference assembly provided as the -F
option. The bam (-B) must be reads aligned to the input assembled fasta (-F).

    fragScaff.pl -B <BAM/bamParse> -F <Input fasta> <Other options>
	
NOTE: It is recommended to run as 'fragScaff.pl -B <BAM> -b 1 -E <node size>' to first generate the bamParse file. This file can then
      be used for all subsequent stages to eliminate the need to re-parse the bam. Similarly it is recommended to run the node-node
      calculations once to generate a '.links.txt' file which can then be used as option '-K' for all subsequent runs thus eliminating the
      need to rerun the time-consuming node-node calculations.

UPDATE: Date: Jan. 12, 2018. Updated to allow 10X genomics bam files by detecting the CR:Z field and building the RG info in pre-processing.
        Note: This has not been fully tested.
	  
RECOMMENDED   (1) 'fragScaff.pl -B <myBAM> -b 1 -E <node_size> (optional: -J <repeats.bed> -N <Nbases.bed> -o <max_node_size>)'
      RUNS:   (2) 'fragScaff.pl -B <myBAM.E#.bamParse> -A -O <myOUT> -t <myTHREAD#/Q>'
              (3) 'fragScaff.pl -B <myBAM.E#.bamParse> -K <myOUT.fragScaff.links.txt> -F <myIN.fa> -O <myOUT2> <other options>'

              This will allow you to run all pre-processing in steps (1) and (2) and then run the last stage (3) multiple times to compare
              outputs and find the optimum setting for your input. Though defaults should give an optimum output.

              Also note that step (3) can be run multiple times without the -F option and with -I so that it will only perform edge filtering
              and graph manipulations which require the most optimization. Once the final otimized parameters are decided based on the stats
              printed in the log file, it can then be run with -F to specify the input assembly fasta and without -I to produce the final
              output assembly.

    -H   Display the following detailed options description.

---------------------- INPUT OPTIONS --------------------------------------------------------------------------------------------------------
    -B   [FILE]   bam file (must have SQ (w/LN) & RG tags) (REQ)        -O   [STR]    output prefix (def = -B)
                      (can be replaced with bamParse file)              -F   [FILE]   fasta file from assembler output
    -P   [T/F/L]  platform:                                             -G   [R/N/H]  read group identifier:
                      T = CPT-Seq / 10X linked read (def)                                 N = read name (group:read_number) (def)
                      F = Fosmid Pool                                                     R = RG:Z:group bam tag
                      L = Long Fragment Read                                              H = after hash (name#group)
    -K   [FILE]   previous link file.                                                     X = CR:Z:group bam tag (10X bam file - untested)
    -N   [FILE]   Nbase bed file for input scaffolds                    -J   [FILE]   Repeatmasker bed file to exclude reads in windows

---------------------- PRUNING OPTIONS ------------------------------------------------------------------------------------------------------
    -q   [INT]    min mapping qual (def = 10)                           -m   [INT]    min contig size to include (def = 1)
    -d   [FLOAT]  min fraction group hit cutoff (def = 0.05)            -l   [INT]    max number of links to use (def = 5)
    -D   [FLOAT]  max fraction group hit cutoff (def = 0.95)            -a   [INT]    max number of links to allow (def = 20)
    -U   [INT]    min num group hits per node (def = 0)                                   (will remove nodes with more)
    -p   [INT/A]  -log10 score use minimum (def = A)                    -j   [FLOAT]  if -p A: mean links per p-bin (def = 1.25)
                      (set to A for auto determine)                     -c   [INT/D]  size of brach to check for split (def = D)
    -r   [INT]    -log10 score report minimum (def = 1)                                   (set to D to disable, experimental)
    -M   [INT]    max score (def = 200)                                 -u   [FLOAT]  combined reciprocation score multiplier (def = 2)
    -g   [INT/X]  max size to make join (ie. one scaff must be
                      less than this, def = X (null))

---------------------- PLATFORM-DEPENDENT OPTIONS -------------------------------------------------------------------------------------------
    -E   [INT]    contig end node size (T/L=5000, F=1000)               -C   [INT]    read count threshold (T=1, F/L=20)
                      (Important for bamParse creation)                 -o   [INT]    max E (if dynamic) (T=10000, F/L=5000)

---------------------- OUTPUT OPTIONS -------------------------------------------------------------------------------------------------------
    -n   [n/N]    case of N to add (def = N)                            -b   [0/1]    generate bamParse file to remove need to
    -L   [INT]    wrap length in output fasta (def = 100)                                 parse bam each run (...bam.bamParse)
    -v            print/keep intermediate files (def = no)                                0 = make bamParse & continue
    -f            print assembly qualities in fasta file                                  1 = make bamParse then exit
                      (def = separate file)                                               (can be run with just -b 1 -B <bam>)
    -R            print out link fraction overlaps (def = no)           -z            gzip bamParse (gzip must be command line callable)
    -x   [INT]    min N spacer size (rec = max MP size, def = 3000)     -X   [INT]    max N spacer size (def = 8000)
    -V   [INT/A]  print [INT] graph files (def = 0, A for all)          -w            exit after printing cytoscape graph (def = no, req -V)

---------------------- PERFORMANCE OPTIONS --------------------------------------------------------------------------------------------------
    -t   [INT/Q]  threads (def = 1)                                     -S   [INT]    if -t Q/>1: number of nodes per job (def = 100)
                      runs with 1 during bam parsing then               -T   [INT]    if -t Q: number of qsub jobs at a time (def = 100)
                      expands to multiple. The mutithread                                 set to 'N' for no limit
                      portion may require a large amount of             -e   [STR]    if -t Q: memory per job (def = 2G)
                      memory. Prepare for 2G each.                      -Q   [STR]    for qsub - do not manually toggle
                  assign to Q to qsub jobs. (best)                      -A            exit after node all-by-all claculations
    -I            exit after graph manipulation
    
---------------------- DEPENDENCY OPTIONS ---------------------------------------------------------------------------------------------------
    -k   [STR]    fragScaff call (def = fragScaff.pl)                   -s   [STR]    samtools call (def = samtools)
                      (only necessary if -t >1 or Q)

";

if (defined $opt{'H'}) {
	die $die2;
}

$_ = 0 for my($pi,$gauss_x,$gauss_m,$gauss_sd,$gauss_p,$i,$ts);
$pi = 3.1415926535;

# QUAL SCORE PRUNING OPTIONS #
if (!defined $opt{'p'}) {$opt{'p'} = "A"}
if (!defined $opt{'j'}) {$opt{'j'} = 1.25}
if (!defined $opt{'r'}) {$opt{'r'} = 1}
my$max_qual = 0;
if (!defined $opt{'M'}) {$opt{'M'} = 200; $max_qual = 1E-200} else {$max_qual = (10**(-1*$opt{'M'}))}

# BEGIN QSUB CHECK #
if (defined $opt{'Q'}) {
	run_qsub();
	exit;
}
## END QSUB CHECK ##

# BEGIN OPTION PARSING AND INITIALIZATION #

# INPUT OPTIONS #
if (!defined $opt{'B'}) {die "\n\tPROVIDE A BAM OR BAMPARSE FILE (OPTION -B)\n$die1"}
if (!defined $opt{'F'}&&(!defined $opt{'b'} || $opt{'b'} < 0.5)&&!defined $opt{'I'}&&!defined $opt{'A'}&&(!defined $opt{'w'}&&defined $opt{'V'})) {die "\n\tPROVIDE A FASTA FILE FOR THE INPUT ASSEMBLY (OPTION -F)\n\tTHIS IS REQUIRED WHEN RUNNING W/O -b 1, -A or -I\n$die1"}
if (!defined $opt{'O'}) {$opt{'O'} = $opt{'B'}}
if (!defined $opt{'P'}) {$opt{'P'} = "T"}
if (!defined $opt{'q'}) {$opt{'q'} = 10}
if (!defined $opt{'G'}) {$opt{'G'} = "N"} elsif ($opt{'G'} !~ /[RNHX]/) {
	die "
    Can not determine group option (-G), please select from:
           N = read name (group:read_number) (def)
           R = RG:Z:group bam tag
           H = after hash (name#group)
		   X = CR:Z:group bam tag (10X bam file)
"}
if ($opt{'P'} =~ /T/i) {
	if (!defined $opt{'E'}) {$opt{'E'} = 5000}
	if (!defined $opt{'o'}) {$opt{'o'} = 10000}
	if (!defined $opt{'C'}) {$opt{'C'} = 1}
} elsif ($opt{'P'} =~ /L/i) {
	if (!defined $opt{'E'}) {$opt{'E'} = 1000}
	if (!defined $opt{'o'}) {$opt{'o'} = 5000}
	if (!defined $opt{'C'}) {$opt{'C'} = 20}
} elsif ($opt{'P'} =~ /F/i) {
	if (!defined $opt{'E'}) {$opt{'E'} = 1000}
	if (!defined $opt{'o'}) {$opt{'o'} = 5000}
	if (!defined $opt{'C'}) {$opt{'C'} = 20}
} else {
	die "
    Can not determine platform option (-P), please select from:
           T = TC-Seq / 10X Linked Read (def)
           F = Fosmid Pool
           L = Long Fragment Read
"}

# CONTIG PARSING OPTIONS #
if (!defined $opt{'m'}) {$opt{'m'} = 1}
if (!defined $opt{'d'}) {$opt{'d'} = 0.05}
if (!defined $opt{'D'}) {$opt{'D'} = 0.95}
if (!defined $opt{'U'}) {$opt{'U'} = 0}

# GRAPH OPTIONS #
if (!defined $opt{'l'}) {$opt{'l'} = 5}
if (!defined $opt{'a'}) {$opt{'a'} = 20}
if (!defined $opt{'c'}) {$opt{'c'} = "D"}
if (!defined $opt{'V'}) {$opt{'V'} = 0}

# OUTPUT OPTIONS #
if (!defined $opt{'L'}) {$opt{'L'} = 100}
if (!defined $opt{'n'}) {$opt{'n'} = "N"}
if (!defined $opt{'x'}) {$opt{'x'} = 3000}
if (!defined $opt{'X'}) {$opt{'X'} = 8000}
if (!defined $opt{'g'}) {$opt{'g'} = "X"}

# PERFORMANCE OPTIONS #
if (!defined $opt{'S'}) {$opt{'S'} = 100}
if (($opt{'S'}/2) =~ /\./) {$opt{'S'}++}
if (!defined $opt{'T'}) {$opt{'T'} = 100}
if (!defined $opt{'e'}) {$opt{'e'} = "2G"}

# DEPENDENCY OPTIONS #
if (!defined $opt{'s'}) {$opt{'s'} = 'samtools'}
if (!defined $opt{'k'}) {$opt{'k'} = 'fragScaff.pl'}

# COMMON VARIABLE DEFINITIONS #
my$PFX = "";
my$dPFX = "";
my$bPFX = "";
my@P = ();
my@H = ();
my$l = "";
my$i = 0;
my$j = 0;

# LOG / WORKING DIR OPTIONS #
if (!defined $opt{'b'} || $opt{'b'} < 0.5) {

	$PFX = $opt{'O'};

	if (-e "$PFX.fragScaff.log") {
		die "\n$PFX.fragScaff.log logfile already exists! Exiting!\n";
	} else {
		open LOG, ">$PFX.fragScaff.log";
		select((select(LOG), $|=1)[0]); # make LOG hot
		$ts = localtime(time);
		print LOG "$ts\tProgram Called.\n";
		
		if ($opt{'B'} =~ /bamParse/) {
			my$bamParseE = 0;
			my@H;
			@H = split(/\./, $opt{'B'});
			foreach (@H) {
				if ($_ =~ /^E\d+$/) {$opt{'E'} = $_; $opt{'E'} =~ s/^E//; $bamParseE = 1}
			}
			if ($bamParseE>0) {
				print LOG "\tbamParse file provided. E detected as $opt{'E'}\n";
			} else {
				print LOG "\tbamParse file provided. E not detected, proceeding anyway.\n";
			}
		}
		
		my$linkFile = "$PFX.fragScaff.links.txt";
		if (defined $opt{'K'}) {
			$linkFile = $opt{'K'};
			open IN, "$linkFile";
			$l = <IN>; chomp $l;
			if ($l =~ /^#LINKFILE/) {
				@P = split(/\t/, $l);
				$opt{'d'} = $P[1];
				$opt{'D'} = $P[2];
				$opt{'E'} = $P[3];
				$opt{'r'} = $P[4];
				print LOG "\tK file provided, forcing options: -d $opt{'d'}, -D $opt{'D'}, -E $opt{'E'}, -r $opt{'r'}. Expected nodeCT = $P[5]\n";
			} else {
				print LOG "\tK file provided, but no #LINKFILE header line! Proceeding anyway.\n";
			}
			close IN;
		}
		
		$dPFX = "$PFX.r$opt{'r'}.fragScaff";
		
		$ts = localtime(time);
		print LOG "$ts\tOptions: ";
		
		foreach my$option (sort keys %opt) {
			if ($option !~ /[vfRzwIAH]/) {
				print LOG " -$option $opt{$option}";
			}
		}
		foreach my$option (sort keys %opt) {
			if ($option =~ /[vfRzwIAH]/) {
				print LOG " -$option";
			}
		}
		
		print LOG "\n";
	}
	
	if (!defined $opt{'K'}) {
		$ts = localtime(time);
		print LOG "$ts\tCreating temporary directory ($dPFX).\n";
		if (-d "$dPFX") {
			print LOG "\n$dPFX directory already exists! Exiting!\n";
			die "\n$dPFX directory already exists! Exiting!\n";
		} else {
			system("mkdir $dPFX");
		}
	}
	
	if (defined $opt{'b'}) {
		my$BPoutName = "$opt{'O'}.E$opt{'E'}.o$opt{'o'}";
		if (defined $opt{'J'}) {$BPoutName .= ".J"}
		if (defined $opt{'N'}) {$BPoutName .= ".N"}
		$bPFX = "$BPoutName.bamParse";
	}
	
} else {
	my$BPoutName = "$opt{'O'}.E$opt{'E'}.o$opt{'o'}";
	if (defined $opt{'J'}) {$BPoutName .= ".J"}
	if (defined $opt{'N'}) {$BPoutName .= ".N"}
	$bPFX = "$BPoutName.bamParse";
	open LOG, ">$bPFX.log";
	select((select(LOG), $|=1)[0]);
}

## END OPTION PARSING AND INITIALIZATION ##

# BEGIN READ IN CONTIGS AND READ GROUPS #

$ts = localtime(time);
print LOG "$ts\tReading in contigs and groups ... \n";

$_ = 0 for my($groupID,$contigID,$nodeID,$length,$groupCT,$contigCT,$nodeCT);
$_ = "" for my($group_name,$contig_name);
my(@GROUP_NAME,@CONTIG_NAME,@NODE_NAME,@NODE_PARTNER_NAME,@CONTIG_LENGTH,@NODE_HITS);
my(%GROUP_ID,%CONTIG_ID,%NODE_ID);

# BEGIN RG IDENTIFICATION SUB #
my($group_found,$fieldID,$null) = (0,0,0);
sub get_rg {
	$group_found = 0;
	if ($opt{'G'} =~ /R/i) {
		for ($fieldID = 11; $fieldID < @P; $fieldID++) {
			if ($P[$fieldID] =~ /^RG:Z:/) {
				$group_name = $P[$fieldID]; $group_name =~ s/^RG:Z://;
				if (defined $GROUP_ID{$group_name}) {
					$groupID = $GROUP_ID{$group_name};
					$group_found = 1;
				}
				$fieldID += 999;
			}
		}
	} elsif ($opt{'G'} =~ /N/i) {
		($group_name,$null) = split(/:/, $P[0]);
		if (defined $GROUP_ID{$group_name}) {
			$groupID = $GROUP_ID{$group_name};
			$group_found = 1;
		}
	} elsif ($opt{'G'} =~ /H/i) {
		($null,$group_name) = split(/#/, $P[0]);
		if (defined $GROUP_ID{$group_name}) {
			$groupID = $GROUP_ID{$group_name};
			$group_found = 1;
		}
	} elsif ($opt{'G'} =~ /X/i) {
		for ($fieldID = 11; $fieldID < @P; $fieldID++) {
			if ($P[$fieldID] =~ /^CR:Z:/) {
				$group_name = $P[$fieldID]; $group_name =~ s/^CR:Z://;
				if (defined $GROUP_ID{$group_name}) {
					$groupID = $GROUP_ID{$group_name};
					$group_found = 1;
				}
				$fieldID += 999;
			}
		}
	}
}
## END RG IDENTIFICATION SUB ##


if ($opt{'B'} =~ /bamParse/) {
	my($groupHits,$readHits) = ("",0);
	if ($opt{'B'} =~ /gz$/) {
		open BP, "zcat $opt{'B'} |";
	} else {
		open BP, "$opt{'B'}";
	}
	while ($l = <BP>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[0] eq "GROUP") {
			$GROUP_NAME[$P[1]] = $P[2];
			$GROUP_ID{$P[2]} = $P[1];
			$groupCT++;
		} elsif ($P[0] eq "CONTIG") {
			$CONTIG_NAME[$P[1]] = $P[2];
			$CONTIG_ID{$P[2]} = $P[1];
			$CONTIG_LENGTH[$P[1]] = $P[3];
			$contigCT++;
		} elsif ($P[0] eq "NODE") {
			$NODE_NAME[$P[1]] = $P[2];
			$NODE_ID{$P[2]} = $P[1];
			$NODE_PARTNER_NAME[$P[1]] = $P[3];
			$nodeCT++;
		} elsif ($P[0] eq "NODE_HITS") {
			@H = split(/,/, $P[2]);
			if (@H>=2) {
				foreach $groupHits (@H) {
					($groupID,$readHits) = split(/:/, $groupHits);
					$NODE_HITS[$P[1]]{$groupID} = $readHits;
				}
			}
		} else {
			die "\nUNRECOGNIZED FIELD IN BAMPARSE! ($P[0])\n";
		}
	} close BP;
	
	$ts = localtime(time);
	print LOG "$ts\tbamParse read: $contigCT contigs ($nodeCT nodes), and $groupCT groups.\n";
	
} else {
	if (defined $opt{'b'}) {
		if (defined $opt{'z'}) {
			open BP, "| gzip > $bPFX";
		} else {
			open BP, ">$bPFX";
		}
	}
	if ($opt{'G'} =~ /X/i) { # ADD IN 10X read group identification
		$ts = localtime(time);
		print LOG "$ts\t10X Genomics bam file specified. Parsing bam to collect all read groups. (beta function)\n";
		open BAM, "$opt{'s'} view $opt{'B'} |";
		while ($l = <BAM>) {
			for ($fieldID = 11; $fieldID < @P; $fieldID++) {
				if ($P[$fieldID] =~ /^CR:Z:/) {
					$group_name = $P[$fieldID]; $group_name =~ s/^CR:Z://;
					$fieldID += 999;
					if (!defined $GROUP_ID{$group_name}) {
						$GROUP_NAME[$groupID] = $group_name;
						$GROUP_ID{$group_name} = $groupID;
						if (defined $opt{'b'}) {print BP "GROUP\t$groupID\t$group_name\n"}
						$groupID++;
					}
				}
			}
		} close BAM;
	}
	open HEAD, "$opt{'s'} view -H $opt{'B'} |";
	while ($l = <HEAD>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[0] eq "\@RG" && $opt{'G'} !~ /X/i) {
			$group_name = $P[1]; $group_name =~ s/^ID://;
			if ($group_name !~ /NO_MATCH/) {
				if (!defined $GROUP_ID{$group_name}) {
					$GROUP_NAME[$groupID] = $group_name;
					$GROUP_ID{$group_name} = $groupID;
					if (defined $opt{'b'}) {print BP "GROUP\t$groupID\t$group_name\n"}
					$groupID++;
				}
			}
		} elsif ($P[0] eq "\@SQ") {
			$length = $P[2];
			$length =~ s/\D//g;
			if ($length>=$opt{'m'}) {
				$contig_name = $P[1]; $contig_name =~ s/^SN://;
				if (!defined $CONTIG_ID{$contig_name}) {
					$CONTIG_NAME[$contigID] = $contig_name;
					$CONTIG_ID{$contig_name} = $contigID;
					$CONTIG_LENGTH[$contigID] = $length;
					if (defined $opt{'b'}) {print BP "CONTIG\t$contigID\t$contig_name\t$length\n"}
					$NODE_ID{"$contigID.L"} = $nodeID;
					$NODE_NAME[$nodeID] = "$contigID.L";
					$NODE_PARTNER_NAME[$nodeID] = "$contigID.R";
					if (defined $opt{'b'}) {print BP "NODE\t$nodeID\t$NODE_NAME[$nodeID]\t$NODE_PARTNER_NAME[$nodeID]\n"}
					$nodeID++;
					$NODE_ID{"$contigID.R"} = $nodeID;
					$NODE_NAME[$nodeID] = "$contigID.R";
					$NODE_PARTNER_NAME[$nodeID] = "$contigID.L";
					if (defined $opt{'b'}) {print BP "NODE\t$nodeID\t$NODE_NAME[$nodeID]\t$NODE_PARTNER_NAME[$nodeID]\n"}
					$nodeID++;
					$contigID++;
				}
			}
		}
	} close HEAD;
	$groupCT = $groupID; $groupID = 0;
	$contigCT = $contigID; $contigID = 0;
	$nodeCT = $nodeID; $nodeID = 0;
	$ts = localtime(time);
	print LOG "$ts\t$contigCT contigs ($nodeCT nodes), and $groupCT groups.\n";
	
	# BEGIN PARSE BAM #
	
	# BEGIN EXCLUSION BED #
	my%EXCLUSIONS; my$excludeWin; my%RIGHT_EXCLUSIONS;
	if (defined $opt{'J'}) {
		$ts = localtime(time);
		if ($opt{'B'} =~ /bamParse/) {
			print LOG "$ts\tExclusion windows defined (-J), but bamParse provided. Ignoring -J.\n";
		} else {
			print LOG "$ts\tReading in exclusion windows (-J) from $opt{'J'}\n";
			$_ = 0 for my($totalExclusionWindows,$totalExclusionBases,$totalExclusionWindowsE,$totalExclusionBasesE,$wrongWindowCT);
			open IN, "$opt{'J'}";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				if (!defined $CONTIG_LENGTH[$CONTIG_ID{$P[0]}]) {
					$wrongWindowCT++;
				} else {
					if ($P[1]<$opt{'E'}||$P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'E'})) {
						$totalExclusionWindowsE++;
						$totalExclusionBasesE+=($P[2]-$P[1]);
					}
					if ($P[1]<$opt{'o'}||$P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'o'})) {
						$EXCLUSIONS{$P[0]}{$P[1]} = $P[2];
						if ($P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'o'})) {
							$RIGHT_EXCLUSIONS{$P[0]}{$P[2]} = $P[1];
						}
					}
					$totalExclusionWindows++;
					$totalExclusionBases+=($P[2]-$P[1]);
				}
			} close IN;
			print LOG "\t$totalExclusionWindows total windows ($totalExclusionBases bp)
\t$totalExclusionWindowsE windows in default nodes ($totalExclusionBasesE bp)
\t$wrongWindowCT windows do not match contigs.\n";
		}
	}
	## END EXCLUSION BED ##
	
	# BEGIN N-BASE BED #
	my%NBASE_GAPS; my%RIGHT_NBASE_GAPS;
	if (defined $opt{'N'}) {
		$ts = localtime(time);
		if ($opt{'B'} =~ /bamParse/) {
			print LOG "$ts\tNbase windows defined (-N), but bamParse provided. Ignoring -N.\n";
		} else {
			print LOG "$ts\tReading in Nbase windows (-N) from $opt{'N'}\n";
			$_ = 0 for my($NbaseWinCount,$NbaseWinSize,$NbaseWinCountE,$NbaseWinSizeE,$NbaseNonWindowCount);
			open IN, "$opt{'N'}";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				if (!defined $CONTIG_LENGTH[$CONTIG_ID{$P[0]}]) {
					$NbaseNonWindowCount++;
				} else {
					if ($P[1]<$opt{'E'}||$P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'E'})) {
						$NbaseWinCountE++;
						$NbaseWinSizeE+=($P[2]-$P[1]);
					}
					if ($P[1]<$opt{'o'}||$P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'o'})) {
						$NBASE_GAPS{$P[0]}{$P[1]} = $P[2];
						if ($P[2]>($CONTIG_LENGTH[$CONTIG_ID{$P[0]}]-$opt{'o'})) {
							$RIGHT_NBASE_GAPS{$P[0]}{$P[2]} = $P[1];
						}
					}
					$NbaseWinCount++;
					$NbaseWinSize+=($P[2]-$P[1]);
				}
			} close IN;
			print LOG "\t$NbaseWinCount total windows ($NbaseWinSize bp)
\t$NbaseWinCountE windows in default nodes ($NbaseWinSizeE bp)
\t$NbaseNonWindowCount windows do not match contigs.\n";
		}
	}
	
	## END N-BASE BED ##
	
	# BEGIN DEFINE NODE BOUNDARIES #
	
	$ts = localtime(time);
	my@LEFT_END; my@RIGHT_START; my$clippedNodeCT = 0; my$clippedBpCT = 0;
	if (defined $opt{'J'} || defined $opt{'N'}) {
		open NODES, ">$bPFX.reDefinedNodes.bed";
		print LOG "$ts\tRedefining node boundaries to accomodate exclusions. ";
		$_ = 0 for my($nodeBaseCount,$nodeNullEnd,$nodeInNull,$nodePosition,$nodeNullBasesTotal);
		for ($contigID = 0; $contigID < $contigCT; $contigID++) {
			$nodeBaseCount = 0; $nodePosition = 0; $nodeNullEnd = 0; $nodeNullBasesTotal = 0;  $nodeInNull = 0;
			while ($nodeBaseCount<$opt{'E'}&&$nodePosition<$opt{'o'}&&$nodePosition<$CONTIG_LENGTH[$contigID]) {
				if (defined $NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition}) {
					if (($nodeInNull>0.5&&$NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition}>$nodeNullEnd) || $nodeInNull<0.5) {
						$nodeNullEnd = $NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition};
						$nodeInNull = 1;
					}
				}
				if (defined $EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition}) {
					if (($nodeInNull>0.5&&$EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition}>$nodeNullEnd) || $nodeInNull<0.5) {
						$nodeNullEnd = $EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition};
						$nodeInNull = 1;
					}
				}
				if ($nodeInNull<0.5) {
					$nodeBaseCount++;
				} elsif ($nodePosition>=$nodeNullEnd) {
					$nodeInNull = 0;
				} else {
					$nodeNullBasesTotal++;
				}
				$nodePosition++;
			}
			$LEFT_END[$contigID] = $nodePosition;
			if ($nodeBaseCount<$opt{'E'}) {
				$clippedNodeCT++;
				$clippedBpCT+=($opt{'E'}-$nodeBaseCount);
			}
			print NODES "$CONTIG_NAME[$contigID]\t1\t$LEFT_END[$contigID]\tSIDE=L;CONTIG_LENGTH=$CONTIG_LENGTH[$contigID];VALID_BP=$nodeBaseCount;NULL_BP=$nodeNullBasesTotal\n";
			
			$nodeBaseCount = 0; $nodePosition = $CONTIG_LENGTH[$contigID]; $nodeNullEnd = $CONTIG_LENGTH[$contigID]; $nodeNullBasesTotal = 0; $nodeInNull = 0;
			while ($nodeBaseCount<$opt{'E'}&&$nodePosition>($CONTIG_LENGTH[$contigID]-$opt{'o'})&&$nodePosition>0) {
				if (defined $RIGHT_NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition}) {
					if (($nodeInNull>0.5&&$RIGHT_NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition}<$nodeNullEnd) || $nodeInNull<0.5) {
						$nodeNullEnd = $RIGHT_NBASE_GAPS{$CONTIG_NAME[$contigID]}{$nodePosition};
						$nodeInNull = 1;
					}
				}
				if (defined $RIGHT_EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition}) {
					if (($nodeInNull>0.5&&$RIGHT_EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition}<$nodeNullEnd) || $nodeInNull<0.5) {
						$nodeNullEnd = $RIGHT_EXCLUSIONS{$CONTIG_NAME[$contigID]}{$nodePosition};
						$nodeInNull = 1;
					}
				}
				if ($nodeInNull<0.5) {
					$nodeBaseCount++;
				} elsif ($nodePosition<=$nodeNullEnd) {
					$nodeInNull = 0;
				} else {
					$nodeNullBasesTotal++;
				}
				$nodePosition--;
			}
			$RIGHT_START[$contigID] = $nodePosition;
			if ($nodeBaseCount<$opt{'E'}) {
				$clippedNodeCT++;
				$clippedBpCT+=($opt{'E'}-$nodeBaseCount);
			}
			print NODES "$CONTIG_NAME[$contigID]\t$RIGHT_START[$contigID]\t$CONTIG_LENGTH[$contigID]\tSIDE=R;CONTIG_LENGTH=$CONTIG_LENGTH[$contigID];VALID_BP=$nodeBaseCount;NULL_BP=$nodeNullBasesTotal\n";
		}
		close NODES;
		print LOG "$clippedNodeCT nodes with less than $opt{'E'} valid bp for a total of $clippedBpCT bp short.\n";
	} else {
		print LOG "$ts\tUsing default node boundaries. ";
		for ($contigID = 0; $contigID < $contigCT; $contigID++) {
			if ($CONTIG_LENGTH[$contigID]<$opt{'E'}) {
				$LEFT_END[$contigID] = $CONTIG_LENGTH[$contigID];
				$RIGHT_START[$contigID] = 1;
				$clippedNodeCT++;
				$clippedBpCT+=(($opt{'E'}-$CONTIG_LENGTH[$contigID])*2);
			} else {
				$LEFT_END[$contigID] = $opt{'E'};
				$RIGHT_START[$contigID] = ($CONTIG_LENGTH[$contigID]-$opt{'E'});
			}
		}
		print LOG "$clippedNodeCT nodes with less than $opt{'E'} bp for a total of $clippedBpCT bp short.\n";
	}
	%NBASE_GAPS = (); %RIGHT_EXCLUSIONS = (); %RIGHT_NBASE_GAPS = ();
	## END DEFINE NODE BOUNDARIES ##
	
	$ts = localtime(time);
	print LOG "$ts\tFinding covered contig ends ...\n";
	
	$_ = 0 for my($no_group_read_ct,$assigned_ct,$read_increment,$read_check,$read_count,
		   $dup_ct,$uniq_frac,$excludeFlag,$excludedReads,$readEnd,$exclude_frac,$no_group_read_frac);
	my@PREV_LOC;
	$read_increment = 10000000; $read_check = $read_increment;
	open IN, "$opt{'s'} view -q $opt{'q'} $opt{'B'} |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		
		get_rg();
		if ($group_found > 0.5) {
			$assigned_ct++;
			
			$excludeFlag = 0;
			if (defined $EXCLUSIONS{$P[2]}) {
				$readEnd = $P[3]+length($P[9]);
				foreach $excludeWin (keys %{$EXCLUSIONS{$P[2]}}) {
					if ($P[3]<=$excludeWin) {
						if ($readEnd>=$EXCLUSIONS{$P[2]}{$excludeWin}) {
							$excludeFlag = 1;
						} elsif ($readEnd>$excludeWin&&$readEnd<=$EXCLUSIONS{$P[2]}{$excludeWin}) {
							if (($readEnd-$excludeWin)/length($P[9])>=0.5) {
								$excludeFlag = 1;
							}
						}
					} elsif ($P[3]<$EXCLUSIONS{$P[2]}{$excludeWin}) {
						if ($readEnd<=$EXCLUSIONS{$P[2]}{$excludeWin}) {
							$excludeFlag = 1;
						} elsif (($readEnd-$EXCLUSIONS{$P[2]}{$excludeWin})/length($P[9])>=0.5) {
							$excludeFlag = 1;
						}
					}
				}
			}
		
			if ($excludeFlag < 0.5) {
				if ($PREV_LOC[$groupID] ne "$P[2]:$P[3]") {
					if (defined $CONTIG_ID{$P[2]}) {
						$contigID = $CONTIG_ID{$P[2]};
						if ($P[3]<=$LEFT_END[$contigID]) {
							$nodeID = $NODE_ID{"$contigID.L"};
							$NODE_HITS[$nodeID]{$groupID}++;
						}
						if (($P[3]+length($P[9]))>=$RIGHT_START[$contigID]) {
							$nodeID = $NODE_ID{"$contigID.R"};
							$NODE_HITS[$nodeID]{$groupID}++;
						}
					}
				} else {
					$dup_ct++;
				}
				$PREV_LOC[$groupID] = "$P[2]:$P[3]";
			} else {
				if ($PREV_LOC[$groupID] eq "$P[2]:$P[3]") {$dup_ct++}
				$PREV_LOC[$groupID] = "$P[2]:$P[3]";
				$excludedReads++;
			}
		} else {
			$no_group_read_ct++;
		}
		
		$read_count++;
		if ($read_count>=$read_check) {
			$no_group_read_frac = sprintf("%.3f", $no_group_read_ct/$read_count);
			$exclude_frac = sprintf("%.3f", $excludedReads/$assigned_ct);
			$uniq_frac = sprintf("%.3f", $dup_ct/$assigned_ct);
			$ts = localtime(time);
			print LOG "\t$ts\t$read_check reads processed. ($no_group_read_ct ($no_group_read_frac) have no group, $excludedReads ($exclude_frac) in repeats, $dup_ct ($uniq_frac) PCR duplicates)\n";
			$read_check+=$read_increment;
		}
	} close IN;
	@PREV_LOC = ();
	
	$no_group_read_frac = sprintf("%.3f", $no_group_read_ct/$read_count);
	$exclude_frac = sprintf("%.3f", $excludedReads/$assigned_ct);
	$uniq_frac = sprintf("%.3f", $dup_ct/$assigned_ct);
	$ts = localtime(time);
	print LOG "$ts\t$no_group_read_ct ($no_group_read_frac) reads had no group, $excludedReads ($exclude_frac) were in repeats, and $dup_ct ($uniq_frac) were PCR duplicates.\n";
	
	if (defined $opt{'b'}) {
		for ($nodeID = 0; $nodeID < $nodeCT; $nodeID++) {
			my$bamParse_out = "NODE_HITS\t$nodeID\t";
			foreach $groupID (keys %{$NODE_HITS[$nodeID]}) {
					$bamParse_out .= "$groupID:$NODE_HITS[$nodeID]{$groupID},";
			}
			$bamParse_out =~ s/,$//;
			print BP "$bamParse_out\n";
		}
		close BP;
		if ($opt{'b'} > 0.5) {
			$ts = localtime(time);
			print LOG "$ts\tbamParse created, exiting.\n";
			close LOG;
			exit;
		}
	}
	@LEFT_END = (); @RIGHT_START = ();
	%EXCLUSIONS = ();
	## END PARSE BAM ##
}
## END READ IN CONTIGS AND READ GROUPS ##

# BEGIN FILTER FOR REAL HITS AND MAKE MASTER FILE #

my($total_pass,$runningNodeCT) = (0,0);
my($minThresh,$maxThresh) = (-1,-1);
my$out_hits = "";
my(@NODE_INCLUSIONS,@TOTAL_PASS_CT);
my%GROUP_HIT_HIST;

$ts = localtime(time);
print LOG "$ts\tFiltering node hits ... ";

setup_node_hits();

sub setup_node_hits {
	for ($nodeID = 0; $nodeID < $nodeCT; $nodeID++) {
		$total_pass = 0;
		foreach $groupID (keys %{$NODE_HITS[$nodeID]}) {
			if ($NODE_HITS[$nodeID]{$groupID} >= $opt{'C'}) {
				$NODE_INCLUSIONS[$nodeID]{$groupID} = 1;
				$total_pass++;
			}
		}
		$TOTAL_PASS_CT[$nodeID] = $total_pass;
		$GROUP_HIT_HIST{$total_pass}++;
	}
	@NODE_HITS = ();
	
	if (defined $opt{'v'}) {open HIST, ">$PFX.fragScaff.group_hits.hist"}
	$runningNodeCT = 0;
	foreach $total_pass (sort {$a<=>$b} keys %GROUP_HIT_HIST) {
		$runningNodeCT+=$GROUP_HIT_HIST{$total_pass};
		if (($runningNodeCT/$nodeCT)>$opt{'d'}&&$minThresh<0) {
			print LOG "Min cut = $total_pass ($runningNodeCT excluded), ";
			$minThresh = $total_pass;
		}
		if (($runningNodeCT/$nodeCT)>$opt{'D'}&&$maxThresh<0) {
			print LOG "Max cut = $total_pass (".($nodeCT-$runningNodeCT)." excluded)\n";
			$maxThresh = $total_pass;
		}
		if (defined $opt{'v'}) {print HIST "$total_pass\t$GROUP_HIT_HIST{$total_pass}\n"}
	}
	if (defined $opt{'v'}) {close HIST}
	
	if (!defined $opt{'K'}) {
		open MST, ">$dPFX/hits.txt";
		for ($nodeID = 0; $nodeID < $nodeCT; $nodeID++) {
			if ($TOTAL_PASS_CT[$nodeID]>$minThresh&&$TOTAL_PASS_CT[$nodeID]<$maxThresh&&$TOTAL_PASS_CT[$nodeID]>=$opt{'U'}) {
				$out_hits = "";
				foreach $groupID (keys %{$NODE_INCLUSIONS[$nodeID]}) {
					$out_hits .= "$groupID,";
				}
				$out_hits =~ s/,$//;
				print MST "$nodeID\t$TOTAL_PASS_CT[$nodeID]\t$out_hits\n";
			}
		}
		close MST;
	}
}
## END FILTER FOR REAL HITS AND MAKE MASTER FILE ##

$ts = localtime(time);
print LOG "$ts\tCalculating link scores ... nodeCT = $nodeCT\n";

# BEGIN GAUSS MASS FUNCTION SUB #
sub gauss_p {
	$gauss_x = $_[0]; $gauss_m = $_[1]; $gauss_sd = $_[2];
	$gauss_p = ((1/sqrt(2*$pi*($gauss_sd**2)))*exp(-1*((($gauss_x-$gauss_m)**2)/(2*($gauss_sd**2)))));
	return $gauss_p;
}
## END GAUSS MASS FUNCTION SUB ##

# BEGIN EXECUTE PROBABILITY CALCULATIONS #
my($nodeID1,$nodeID2) = (0,0);
my%LOAD_FILES; my%NODE2;
if (!defined $opt{'K'}) {
	if (!defined $opt{'t'} || ($opt{'t'} =~ /\d/ && $opt{'t'} < 1.5)) {

		$ts = localtime(time);
		print LOG "$ts\tSingle thread mode. Will be slow.\n";
		
		pvalue_links(0,($nodeCT-1));
		$LOAD_FILES{"0"} = "$dPFX/node.0.".($nodeCT-1).".txt";
		
	} elsif ($opt{'t'} eq "Q") {

		$ts = localtime(time);
		print LOG "$ts\tQsubbing mode, submitting jobs. ";
		
		open IN, "pwd |"; my$PWD = <IN>; close IN; chomp $PWD;
		my$jobMax = 0;
		open JOBARRAY, ">$dPFX/job_array.txt";
		my(%QSUB_STATUS); my($complete_ct,$incomplete_ct,$nodeID1,$nodeID2) = (0,0,0,0);
		for ($nodeID1 = 0; $nodeID1 < $nodeCT; $nodeID1+=$opt{'S'}) {
			$incomplete_ct++;
			$nodeID2 = $nodeID1+($opt{'S'}-1);
			if ($nodeID2>=$nodeCT) {$nodeID2 = ($nodeCT-1)}
			$NODE2{$nodeID1} = $nodeID2;
			$QSUB_STATUS{$nodeID1} = 0;
			if (defined $opt{'R'}) {
				print JOBARRAY "$opt{'k'} -Q $PWD/$dPFX,$nodeCT,$nodeID1,$nodeID2 -r $opt{'r'} -M $opt{'M'} -R\n";
			} else {
				print JOBARRAY "$opt{'k'} -Q $PWD/$dPFX,$nodeCT,$nodeID1,$nodeID2 -r $opt{'r'} -M $opt{'M'}\n";
			}
			$LOAD_FILES{$jobMax} = "$dPFX/node.$nodeID1.$nodeID2.txt";
			$jobMax++;
		}
		close JOBARRAY; 
		
		my$qsubCommand = "$PWD/$dPFX/run_array.csh";
		open RUN_ARRAY, ">$qsubCommand";
		print RUN_ARRAY "#/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n#\$ -V\nCOMMAND=\$(head -n \$SGE_TASK_ID $PWD/$dPFX/job_array.txt | tail -n 1)\n\$COMMAND\n";
		close RUN_ARRAY; system("chmod +x $qsubCommand");
		
		my$jobName = "FSCF_".time;
		print LOG "Job array name = $jobName, Job count = $jobMax\n";
		
		if ($opt{'T'} eq "N") {
			system("qsub -t 1-$jobMax -N $jobName -b y -l h_vmem=$opt{'e'},virtual_free=$opt{'e'} $qsubCommand \"$PWD/$dPFX/job_array.txt\"");
		} else {
			system("qsub -t 1-$jobMax -N $jobName -b y -l h_vmem=$opt{'e'},virtual_free=$opt{'e'} -tc $opt{'T'} $qsubCommand \"$PWD/$dPFX/job_array.txt\"");
		}
		
		while ($complete_ct<$incomplete_ct) {
			$complete_ct = 0;
			foreach my$checkNodeID (keys %QSUB_STATUS) {
				$nodeID2 = $NODE2{$checkNodeID};
				if ($QSUB_STATUS{$checkNodeID} == 0 && -e "$dPFX/node.$checkNodeID.$nodeID2.complete") {
					$QSUB_STATUS{$checkNodeID} = 1;
					system("rm -f $dPFX/node.$checkNodeID.$nodeID2.complete");
				}
				if ($QSUB_STATUS{$checkNodeID} == 1) {
					$complete_ct++;
				}
			}
			sleep(1);
			if (-e "rm -f $dPFX/kill") {$complete_ct=$incomplete_ct}
		}
	} else {
	
		$ts = localtime(time);
		print LOG "$ts\tRunning with $opt{'t'} threads.\n";
		
		my(%THREAD_STATUS); my($incomplete_ct,$complete_ct,$runningThreads) = (0,0,0);
		($nodeID1,$nodeID2) = (0,0);
		for ($nodeID1 = 0; $nodeID1 < $nodeCT; $nodeID1+=$opt{'S'}) {
			$THREAD_STATUS{$nodeID1} = 0;
			$nodeID2 = $nodeID1+($opt{'S'}-1);
			if ($nodeID2>=$nodeCT) {$nodeID2 = ($nodeCT-1)}
			$NODE2{$nodeID1} = $nodeID2;
			$LOAD_FILES{$incomplete_ct} = "$dPFX/node.$nodeID1.$nodeID2.txt";
			$incomplete_ct++;
		}
		while ($complete_ct<$incomplete_ct) {
			for ($nodeID1 = 0; $nodeID1 < $nodeCT; $nodeID1+=$opt{'S'}) {
				$nodeID2 = $NODE2{$nodeID1};
				if (-e "$dPFX/node.$nodeID1.$nodeID2.complete") {
					if ($THREAD_STATUS{$nodeID1} == 1) {
						system("rm -f $dPFX/node.$nodeID1.$nodeID2.complete");
						$runningThreads--;
						$complete_ct++;
					}
					$THREAD_STATUS{$nodeID1} = 2;
				} elsif ($THREAD_STATUS{$nodeID1} == 0 && $runningThreads < $opt{'t'}) {
					$nodeID2 = $NODE2{$nodeID1};
					if (defined $opt{'R'}) {
						system("$opt{'k'} -Q $dPFX,$nodeCT,$nodeID1,$nodeID2 -r $opt{'r'} -M $opt{'M'} -R &");
					} else {
						system("$opt{'k'} -Q $dPFX,$nodeCT,$nodeID1,$nodeID2 -r $opt{'r'} -M $opt{'M'} &");
					}
					$THREAD_STATUS{$nodeID1} = 1;
					$runningThreads++;
				}
			}
			sleep(1);
			if (-e "rm -f $dPFX/kill") {$complete_ct=$incomplete_ct}
		}
	}
	$ts = localtime(time);
	print LOG "$ts\tAll link pvalues calculated.\n";
} else {
	$ts = localtime(time);
	print LOG "$ts\tK toggled, using $opt{'K'} for links.\n";
}
## END EXECUTE PROBABILITY CALCULATIONS ##

# BEGIN PVALUE LINKS SUBROUTINE #

my(@PASS_HITS,@GROUP_LIST,@TOTAL_PASS);
my($startNode,$endNode,$nodeID1,$nodeID2,$mean_num,$stdev_num,$tally,$shared,$total,$frac,$mean,
   $stdev,$pval,$side,$score,$total_links,$pair_frac,$calc_nodeID1);
my($nodeOut,$calcNodeOut) = ("","");
my%FRAC;
sub pvalue_links {

	($startNode,$endNode) = (0,0);
	$startNode = $_[0];
	$endNode = $_[1];
	
	(@PASS_HITS,@GROUP_LIST,@TOTAL_PASS) = ();
	
	open MST, "$dPFX/hits.txt";
	while ($l = <MST>) {
		chomp $l;
		@P = split(/\t/, $l);
		$TOTAL_PASS[$P[0]] = $P[1];
		@GROUP_LIST = split(/,/, $P[2]);
		foreach $groupID (@GROUP_LIST) {
			$PASS_HITS[$P[0]]{$groupID} = 1;
		}
	} close MST;
	
	$_ = 0 for ($nodeID1,$nodeID2,$mean_num,$stdev_num,$tally,$shared,$total,$frac,$mean,$stdev,
	            $pval,$side,$total_links,$pair_frac);
	
	$nodeOut = "";
	for ($nodeID1 = $startNode; $nodeID1 <= $endNode; $nodeID1++) {
		if ($TOTAL_PASS[$nodeID1]>0){$nodeOut .= calc_node($nodeID1,0)}
		$nodeID1++; 
		if ($TOTAL_PASS[$nodeID1]>0){$nodeOut .= calc_node($nodeID1,1)}
	}
	
	open OUT, ">$dPFX/node.$startNode.$endNode.txt";
	print OUT "$nodeOut"; close OUT;
}
## END PVALUE LINKS SUBROUTINE ##

# BEGIN CALC PVAL SUB #
sub calc_node {
	
	$calc_nodeID1 = $_[0];
	$side = $_[1];
	$_ = 0 for ($mean_num,$stdev_num,$tally,$shared,$total,$frac,$mean,$stdev,$pval,$score,
	            $total_links,$pair_frac);
	$calcNodeOut = "";
	
	%FRAC = ();
	for ($nodeID2 = 0; $nodeID2 < $nodeCT; $nodeID2++) {
		if ($side < 0.5) {
			if ($nodeID2 != $calc_nodeID1 && $nodeID2 != ($calc_nodeID1+1)) {
				$shared = 0;
				foreach $groupID (keys %{$PASS_HITS[$calc_nodeID1]}) {
					if (defined $PASS_HITS[$nodeID2]{$groupID}) {
						$shared++;
					}
				}
				if ($shared > 0) {
					$total = $TOTAL_PASS[$calc_nodeID1]+$TOTAL_PASS[$nodeID2]-$shared;
					if ($total>0) {
						$frac = $shared/$total;
					} else {
						$frac = 1;
					}
					$FRAC{$nodeID2} = $frac;
					$mean_num+=$frac;
					$tally++;
				}
			} elsif ($nodeID2 == ($calc_nodeID1+1)) {
				$shared = 0;
				foreach $groupID (keys %{$PASS_HITS[$calc_nodeID1]}) {
					if (defined $PASS_HITS[$nodeID2]{$groupID}) {
						$shared++;
					}
				}
				if ($shared > 0) {
					$total = $TOTAL_PASS[$calc_nodeID1]+$TOTAL_PASS[$nodeID2]-$shared;
					if ($total>0) {
						$pair_frac = $shared/$total;
					} else {
						$pair_frac = 1;
					}
				} else {
					$pair_frac = 0;
				}
			}
		} else {
			if ($nodeID2 != $calc_nodeID1 && $nodeID2 != ($calc_nodeID1-1)) {
				$shared = 0;
				foreach $groupID (keys %{$PASS_HITS[$calc_nodeID1]}) {
					if (defined $PASS_HITS[$nodeID2]{$groupID}) {
						$shared++;
					}
				}
				if ($shared > 0) {
					$total = $TOTAL_PASS[$calc_nodeID1]+$TOTAL_PASS[$nodeID2]-$shared;
					if ($total>0) {
						$frac = $shared/$total;
					} else {
						$frac = 1;
					}
					$FRAC{$nodeID2} = $frac;
					$mean_num+=$frac;
					$tally++;
				}
			}
		}
	}
	if ($tally>0) {
		$mean = $mean_num/$tally;
		foreach $nodeID2 (keys %FRAC) {
			$stdev_num += ($mean-$FRAC{$nodeID2})**2;
		}
		$stdev = sqrt($stdev_num/$tally);
		if (defined $opt{'R'}) {open FRACS, ">$dPFX/fragListNode$calc_nodeID1.txt"}
		$calcNodeOut .= "#NODE_ID\t$calc_nodeID1\n#HITS\t$TOTAL_PASS[$calc_nodeID1]\n#MEAN\t$mean\n#STDEV\t$stdev\n";
		if ($side < 0.5) {$calcNodeOut .= "#PAIR_FRAC\t$pair_frac\n"}
		foreach $nodeID2 (keys %FRAC) {
			if ($FRAC{$nodeID2}!=0&&$mean!=0&&$stdev!=0) {
				$pval = gauss_p($FRAC{$nodeID2},$mean,$stdev);
			} else {
				$pval = 1;
			}
			if ($pval > $max_qual) {
				$score = -1*(log($pval)/log(10));
			} else {
				$score = $opt{'M'};
			}
			if ($score >= $opt{'r'} && $FRAC{$nodeID2}>$mean) {
				$calcNodeOut .= "$nodeID2\t$TOTAL_PASS[$nodeID2]\t$pval\t$score\n";
				$total_links++;
			}
			if (defined $opt{'R'}) {print FRACS "$nodeID2\t$TOTAL_PASS[$nodeID2]\t$FRAC{$nodeID2}\t$pval\t$score\n"}
		}
		if (defined $opt{'R'}) {close FRACS}
	}
	return $calcNodeOut;
}
## END CALC PVAL SUB ##

# BEGIN READ IN LINK SCORES #
$ts = localtime(time);
print LOG "$ts\tLoading node information ... ";
my($includedCT,$trimmedCT,$excludedCT,$scoreBin,$binCT) = (0,0,0,0,0);
my$recommended_p = -1;
my$linkFile = "$dPFX.r$opt{'r'}.links.txt";
if (defined $opt{'A'} && $opt{'B'} =~ /bamParse/) {
	$linkFile = $opt{'B'};
	$linkFile =~ s/bamParse//;
	$linkFile .= "r$opt{'r'}.links.txt";
}
my$nodeFile = "";
my(%PRE_LINK,%LINK,%PAIR_FRAC,%INCLUSION_STATS);
my(@MEAN,@STDEV,@TOTAL_LINKS,@LINK_BIN_COUNTS);
for ($i = 0; $i <= $opt{'M'}; $i++) {$LINK_BIN_COUNTS[$i] = 0}
if (!defined $opt{'K'}) {
	open LINKFILE, ">$linkFile";
	print LINKFILE "#LINKFILE\t$opt{'d'}\t$opt{'D'}\t$opt{'E'}\t$opt{'r'}\t$nodeCT\n";
	foreach $nodeFile (sort {$a<=>$b} keys %LOAD_FILES) {
		if (-e "$LOAD_FILES{$nodeFile}") {
			open IN, "$LOAD_FILES{$nodeFile}";
			while ($l = <IN>) {
				chomp $l;
				print LINKFILE "$l\n";
				if ($l =~ /#NODE_ID/) {
					$binCT++;
				} elsif ($l !~ /^#/) {
					@P = split(/\t/, $l);
					$scoreBin = int($P[3]);
					for ($i = $scoreBin; $i > 0; $i--) {$LINK_BIN_COUNTS[$i]++}
				}
			}
			close IN;
		} else {
			print LOG "\n\tWARNING: CAN NOT OPEN $LOAD_FILES{$nodeFile}!\n\t";
		}
	}
	close LINKFILE;
	
	open P_BINS, ">$dPFX.links.mean.txt";
	for ($i = 1; $i <= $opt{'M'}; $i++) {
		$mean = sprintf("%.2f", $LINK_BIN_COUNTS[$i]/$binCT);
		print P_BINS "$i\t$mean\n";
		if ($recommended_p < 0 && $mean <= $opt{'j'}) {
			$recommended_p = $i;
		}
	}
	close P_BINS;
	
	if ($opt{'p'} =~ /A/i) {
		print LOG "Auto-p-score: $recommended_p\n";
		$opt{'p'} = $recommended_p;
	} else {
		print LOG "Rec-p-score: $recommended_p (set is $opt{'p'})\n";
	}
} else {
	$linkFile = $opt{'K'};
	if ($opt{'p'} =~ /A/i) {
		print LOG "-p set to A, auto-determining ";
		open IN, "$linkFile";
		while ($l = <IN>) {
			chomp $l;
			if ($l =~ /^#NODE_ID/) {
				$binCT++;
			} elsif ($l !~ /^#/) {
				@P = split(/\t/, $l);
				$scoreBin = int($P[3]);
				for ($i = $scoreBin; $i > 0; $i--) {$LINK_BIN_COUNTS[$i]++}
			}
		}
		close IN;
		
		open P_BINS, ">$dPFX.links.mean.txt";
		for ($i = 1; $i <= $opt{'M'}; $i++) {
			$mean = sprintf("%.2f", $LINK_BIN_COUNTS[$i]/$binCT);
			print P_BINS "$i\t$mean\n";
			if ($recommended_p < 0 && $mean <= $opt{'j'}) {
				$recommended_p = $i;
			}
		}
		close P_BINS;
		$opt{'p'} = $recommended_p;
		print LOG "(p = $opt{'p'})\n";
	} else {
		print LOG "p-score set to $opt{'p'}\n";
	}
}

if (!defined $opt{'v'}) {system("rm -f -R $dPFX/")}

if (defined $opt{'A'}) {
	$ts = localtime(time);
	print LOG "Option A toggled, exiting after node calculations.\n";
	close LOG;
	exit;
}

$nodeID1 = -1;
if ($linkFile =~ /\.gz$/) {
	open IN, "zcat $linkFile |";
} else {
	open IN, "$linkFile";
}
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($l =~ /^#/) {
		if ($P[0] =~ /#LINKFILE/) {
			if ($opt{'d'} != $P[1]) {print LOG "\n\tWARNING: -d option ($P[1]) in linkfile ($linkFile) does not equal provided -d ($opt{'d'}). Continuing anyway ...\n"};
			if ($opt{'D'} != $P[2]) {print LOG "\n\tWARNING: -D option ($P[2]) in linkfile ($linkFile) does not equal provided -D ($opt{'D'}). Continuing anyway ...\n"};
			if ($opt{'E'} != $P[3]) {print LOG "\n\tWARNING: -E option ($P[3]) in linkfile ($linkFile) does not equal provided -E ($opt{'E'}). Continuing anyway ...\n"};
			if ($opt{'p'} < $P[4]) {print LOG "\n\tWARNING: -r option ($P[4]) in linkfile ($linkFile) is greater than -p provided ($opt{'p'}). Continuing with -p effectively at $P[4] ...\n"};
			if ($nodeCT != $P[5]) {print LOG "\n\tFATAL ERROR: linkfile ($linkFile) node count ($P[5]) != node count ($nodeCT) determined from bamParse and -d/-D filtering.\n"; die};
		} elsif ($P[0] =~ /#NODE_ID/) {
			if ($nodeID1 > -0.5) {
				$TOTAL_LINKS[$nodeID1] = 0;
				foreach $nodeID2 (keys %PRE_LINK) {
					$TOTAL_LINKS[$nodeID1]++;
				}
				if ($TOTAL_LINKS[$nodeID1]<=$opt{'l'}) {
					foreach $nodeID2 (keys %PRE_LINK) {
						$LINK{$nodeID1}{$nodeID2} = $PRE_LINK{$nodeID2};
					}
					$INCLUSION_STATS{$nodeID1} = 0;
				} elsif ($TOTAL_LINKS[$nodeID1]<=$opt{'a'}&&$TOTAL_LINKS[$nodeID1]>$opt{'l'}) {
					$includedCT = 0;
					foreach $nodeID2 (sort {$PRE_LINK{$b}<=>$PRE_LINK{$a}} keys %PRE_LINK) {
						if ($includedCT<$opt{'l'}) {
							$LINK{$nodeID1}{$nodeID2} = $PRE_LINK{$nodeID2};
							$includedCT++;
						}
					}
					$INCLUSION_STATS{$nodeID1} = 1;
				} elsif ($TOTAL_LINKS[$nodeID1]>=$opt{'a'}) {
					$INCLUSION_STATS{$nodeID1} = 2;
				}
			}
			$nodeID1 = $P[1];
			%PRE_LINK = ();
		} elsif ($P[0] =~ /#MEAN/) {
			$MEAN[$nodeID1] = $P[1];
		} elsif ($P[0] =~ /#STDEV/) {
			$STDEV[$nodeID1] = $P[1];
		} elsif ($P[0] =~ /#PAIR_FRAC/) {
			$PAIR_FRAC{$NODE_NAME[$nodeID1]} = $P[1];
		}
	} else {
		if ($P[3] >= $opt{'p'}) {
			$PRE_LINK{$P[0]} = $P[3];
		}
	}
} close IN;

$includedCT = 0;
foreach $nodeID1 (keys %INCLUSION_STATS) {
	if ($INCLUSION_STATS{$nodeID1}==0) {$includedCT++}
	elsif ($INCLUSION_STATS{$nodeID1}==1) {$trimmedCT++; $includedCT++}
	elsif ($INCLUSION_STATS{$nodeID1}==2) {$excludedCT++}
}
$ts = localtime(time);
print LOG "$ts\t$includedCT included nodes of which $trimmedCT were trimmed, and $excludedCT were excluded. ";
%INCLUSION_STATS = ();
## END READ IN LINK SCORES ##

# BEGIN RECIPROCATE LINKS #
my$recip = 2*$opt{'p'};
if (defined $opt{'u'}) {$recip = $opt{'u'}*$opt{'p'}}
if ($recip > 1.5*$opt{'M'}) {$recip = 1.5*$opt{'M'}}
print LOG "recip. total set to $recip\n";

$ts = localtime(time);
print LOG "$ts\tBuilding graph ... ";

my($vertex_ct,$edge_ct,$contig_ct);
my(%EDGES,%VERTEX);

foreach $nodeID1 (keys %LINK) {
	foreach $nodeID2 (keys %{$LINK{$nodeID1}}) {
		if (defined $LINK{$nodeID2}{$nodeID1}) {
			
			if (($LINK{$nodeID1}{$nodeID2}+$LINK{$nodeID2}{$nodeID1}) >= $recip) {
			
				if (!defined $VERTEX{$NODE_NAME[$nodeID1]}) {
					$EDGES{$NODE_NAME[$nodeID1]}{$NODE_PARTNER_NAME[$nodeID1]} = 1000;
					$EDGES{$NODE_PARTNER_NAME[$nodeID1]}{$NODE_NAME[$nodeID1]} = 1000;
					$VERTEX{$NODE_NAME[$nodeID1]} = 1; $VERTEX{$NODE_PARTNER_NAME[$nodeID1]} = 1;
					$vertex_ct += 2; $edge_ct++; $contig_ct++;
				}
				
				if (!defined $VERTEX{$NODE_NAME[$nodeID2]}) {
					$EDGES{$NODE_NAME[$nodeID2]}{$NODE_PARTNER_NAME[$nodeID2]} = 1000;
					$EDGES{$NODE_PARTNER_NAME[$nodeID2]}{$NODE_NAME[$nodeID2]} = 1000;
					$VERTEX{$NODE_NAME[$nodeID2]} = 1; $VERTEX{$NODE_PARTNER_NAME[$nodeID2]} = 1;
					$vertex_ct += 2; $edge_ct++; $contig_ct++;
				}
				
				if (!defined $EDGES{$NODE_NAME[$nodeID2]}{$NODE_NAME[$nodeID1]}) {
					$EDGES{$NODE_NAME[$nodeID1]}{$NODE_NAME[$nodeID2]} = ($LINK{$nodeID1}{$nodeID2}+$LINK{$nodeID2}{$nodeID1});
					$EDGES{$NODE_NAME[$nodeID2]}{$NODE_NAME[$nodeID1]} = ($LINK{$nodeID1}{$nodeID2}+$LINK{$nodeID2}{$nodeID1});
					$edge_ct++;
				}
			}
		}
	}
}
%LINK = ();
$ts = localtime(time);
print LOG "$vertex_ct vertices, $edge_ct edges, $contig_ct contigs included\n";
## END RECIPROCATE LINKS ##

# BEGIN SUBGRAPH IDENTIFICATION #
$ts = localtime(time);
print LOG "$ts\tIdentifying MST subgraphs... ";

$_ = 0 for my($included_vertex_ct,$scaffoldID,$vertices_added,$max_score,$scaffold_ct,$added_in_round);
$_ = "" for my($vertex1,$vertex2,$vertex3,$max_vertex,$max_vertex2);
my(%SUBSET_VERTICES,%SUBSET_EDGES);

$added_in_round = 1;
while ($added_in_round>0) {
	$added_in_round = 0;
	foreach $vertex1 (keys %VERTEX) {
		if ($VERTEX{$vertex1} > 0) {
			$VERTEX{$vertex1} = -1;
			$SUBSET_VERTICES{$scaffoldID}{$vertex1} = 1;
			$vertex2 = $NODE_PARTNER_NAME[$NODE_ID{$vertex1}];
			$VERTEX{$vertex2} = -1;
			$SUBSET_VERTICES{$scaffoldID}{$vertex2} = 1;
			$SUBSET_EDGES{$scaffoldID}{$vertex1}{$vertex2} = 1000;
			$SUBSET_EDGES{$scaffoldID}{$vertex2}{$vertex1} = 1000;
			$included_vertex_ct++;
			$added_in_round++;
			$included_vertex_ct++;
			$added_in_round++;
			$vertices_added = 1;
			while ($vertices_added > 0) { # expansion portion (explores to find subgraphs)
				$vertices_added = 0;
				$max_score = 0; $max_vertex = ""; $max_vertex2 = "";
				foreach $vertex2 (keys %{$SUBSET_VERTICES{$scaffoldID}}) {
					foreach $vertex3 (keys %{$EDGES{$vertex2}}) {
						if (!defined $SUBSET_VERTICES{$scaffoldID}{$vertex3} && $VERTEX{$vertex3} > 0) {
							if ($EDGES{$vertex2}{$vertex3} > $max_score) {
								$max_score = $EDGES{$vertex2}{$vertex3};
								$max_vertex = $vertex3;
								$max_vertex2 = $vertex2;
							}
						}
					}
				}
				if ($max_score > 0) {
					$vertices_added++;
					$SUBSET_VERTICES{$scaffoldID}{$max_vertex} = 1;
					$SUBSET_EDGES{$scaffoldID}{$max_vertex}{$max_vertex2} = $EDGES{$max_vertex2}{$max_vertex};
					$SUBSET_EDGES{$scaffoldID}{$max_vertex2}{$max_vertex} = $EDGES{$max_vertex2}{$max_vertex};
					$VERTEX{$max_vertex} = -1;
					$added_in_round++;
					$included_vertex_ct++;
					$vertex3 = $NODE_PARTNER_NAME[$NODE_ID{$max_vertex}];
					$VERTEX{$vertex3} = -1;
					$SUBSET_VERTICES{$scaffoldID}{$vertex3} = 1;
					$SUBSET_EDGES{$scaffoldID}{$max_vertex}{$vertex3} = 1000;
					$SUBSET_EDGES{$scaffoldID}{$vertex3}{$max_vertex} = 1000;
					$included_vertex_ct++;
					$added_in_round++;
				}
			}
			$scaffoldID++;
		}
	}
} $scaffold_ct = $scaffoldID;
print LOG "included $included_vertex_ct of $vertex_ct vertices to $scaffold_ct scaffolds.\n";
## END SUBGRAPH IDENTIFICATION ##

# BEGIN MST STUFF #
$ts = localtime(time);
print LOG "$ts\tBuilding trunks via graph walking.\n";

$_ = 0 for my($degree_one_CT,$total_vertex_CT,$degree,$longest_trunk_length,$end_points,$trunk_length,
              $end_degree,$checked_vertices,$degree_scaffoldID,
			  $makeMST_scaffoldID,$split_scaffoldID,$trunkScaffoldID,$new_scaffoldID1,$new_scaffoldID2,
			  $path_check_max,$fake_degree_one_ct,$fake_degree_one_created);
$_ = "" for my($longest_trunk_path,$start_vertex,$trunk_path,$check_vertex,$next_vertex,$end_vertex,
               $current_vertex);
my(%DEGREE_ONE,%DEGREE,%VERTEX_CT,%DECISIONS,%PATH,%PATH_LENGTH,%TRUNK_VERTEX_CT,%TRUNK,%DEGREE_ONE_CT);

foreach $scaffoldID (keys %SUBSET_VERTICES) {calc_degrees($scaffoldID); TRUNK($scaffoldID)}

sub calc_degrees {
	$degree_scaffoldID = $_[0];
	%DEGREE_ONE = ();
	$degree_one_CT = 0;
	$total_vertex_CT = 0;
	foreach $vertex1 (keys %{$SUBSET_VERTICES{$degree_scaffoldID}}) {
		$total_vertex_CT++;
		$degree = 0;
		foreach $vertex2 (keys %{$SUBSET_EDGES{$degree_scaffoldID}{$vertex1}}) {$degree++};
		$DEGREE{$vertex1} = $degree;
		if ($degree == 1) {
			$DEGREE_ONE{$vertex1} = 1;
			$degree_one_CT++;
		}
	}
	$DEGREE_ONE_CT{$degree_scaffoldID} = $degree_one_CT;
	$VERTEX_CT{$degree_scaffoldID} = $total_vertex_CT;
}

sub TRUNK {
	$trunkScaffoldID = $_[0];
	$longest_trunk_length = 0;
	$longest_trunk_path = "";
	foreach $start_vertex (keys %DEGREE_ONE) {
		%PATH_LENGTH = ();
		%PATH = ();
		%DECISIONS = ();
		$checked_vertices = 1;
		$PATH_LENGTH{$start_vertex} = 1;
		$PATH{$start_vertex} = "$start_vertex";
		$DECISIONS{$start_vertex} = 1;
		while ($checked_vertices<$VERTEX_CT{$trunkScaffoldID}) {
			foreach $current_vertex (keys %DECISIONS) {
				if ($DECISIONS{$current_vertex} > 0) {
					foreach $check_vertex (keys %{$SUBSET_EDGES{$trunkScaffoldID}{$current_vertex}}) {
						if (!defined $PATH{$check_vertex}) {
							$PATH{$check_vertex} = "$PATH{$current_vertex},$check_vertex";
							$PATH_LENGTH{$check_vertex} = $PATH_LENGTH{$current_vertex}+1;
							$DECISIONS{$check_vertex} = $DEGREE{$check_vertex}-1;
							$checked_vertices++;
						}
					}
				}
			}
		}
		foreach $current_vertex (keys %DEGREE_ONE) {
			if ($PATH_LENGTH{$current_vertex}>$longest_trunk_length) {
				$longest_trunk_length = $PATH_LENGTH{$current_vertex};
				$longest_trunk_path = $PATH{$current_vertex};
			}
		}
		if ($DEGREE_ONE_CT{$trunkScaffoldID}<2) {
			foreach $current_vertex (keys %{$SUBSET_VERTICES{$trunkScaffoldID}}) {
				if ($PATH_LENGTH{$current_vertex}>$longest_trunk_length) {
					$longest_trunk_length = $PATH_LENGTH{$current_vertex};
					$longest_trunk_path = $PATH{$current_vertex};
				}
			}
		}
	}
	$TRUNK_VERTEX_CT{$trunkScaffoldID} = $longest_trunk_length;
	$TRUNK{$trunkScaffoldID} = $longest_trunk_path;
}
## END MST STUFF ##

# BEGIN VERBOSE MST PRINTING #
if (defined $opt{'v'}) {
	open TRUNKS, ">$PFX.fragScaff.trunks.txt";
	foreach $scaffoldID (keys %TRUNK) {
		print TRUNKS "$scaffoldID\t$TRUNK_VERTEX_CT{$scaffoldID}\t$TRUNK{$scaffoldID}\n";
	}
	close TRUNKS;
}
## END VERBOSE MST PRINTING ##

# BEGIN CYTOSCAPE PRINTOUT #
if ($opt{'V'} =~ /A/i) {$opt{'V'} = $scaffold_ct*2}
if ($opt{'V'}>0.5) {
	$ts = localtime(time);
	print LOG "$ts\t-V selected, printing out cytoscape file ... ";
	my(%CYTO_OUT,%CYTO_DONE); my$cyto_printed_scaff = 0; my@CYTO_VERTICES;
	foreach $vertex1 (keys %EDGES) {
		foreach $vertex2 (keys %{$EDGES{$vertex1}}) {
			$CYTO_OUT{$vertex1}{$vertex2} = "EDGE";
		}
	}
	foreach $scaffoldID (keys %SUBSET_EDGES) {
		foreach $vertex1 (keys %{$SUBSET_EDGES{$scaffoldID}}) {
			foreach $vertex2 (keys %{$SUBSET_EDGES{$scaffoldID}{$vertex1}}) {
				$CYTO_OUT{$vertex1}{$vertex2} = "MST";
			}
		}
	}
	foreach $scaffoldID (keys %TRUNK) {
		@CYTO_VERTICES = split(/,/, $TRUNK{$scaffoldID});
		foreach $vertex1 (@CYTO_VERTICES) {
			foreach $vertex2 (@CYTO_VERTICES) {
				if ($CYTO_OUT{$vertex1}{$vertex2} eq "MST") {
					$CYTO_OUT{$vertex1}{$vertex2} = "TRUNK";
				}
			}
		}
	}
	open CYTO, ">$PFX.fragScaff.cyto.csv";
	print CYTO "vertex_1,edge_type,vertex_2,edge_weight,group_ID\n";
	foreach $scaffoldID (keys %SUBSET_VERTICES) {
		if ($cyto_printed_scaff<=$opt{'V'}) {
			foreach $vertex1 (keys %{$SUBSET_VERTICES{$scaffoldID}}) {
				foreach $vertex2 (keys %{$CYTO_OUT{$vertex1}}) {
					if (!defined $CYTO_DONE{$vertex1}{$vertex2} && !defined $CYTO_DONE{$vertex2}{$vertex1}) {
						print CYTO "$vertex1,$CYTO_OUT{$vertex1}{$vertex2},$vertex2,$EDGES{$vertex1}{$vertex2},$scaffoldID\n";
						$CYTO_DONE{$vertex1}{$vertex2} = 1;
						$CYTO_DONE{$vertex2}{$vertex1} = 1;
					}
				}
			}
		}
		$cyto_printed_scaff++;
	} close CYTO;
	(%CYTO_OUT,%CYTO_DONE) = ();
	print LOG "done.\n";
	if (defined $opt{'w'}) {
		$ts = localtime(time);
		print LOG "$ts\t-w selected, exiting after cytoscape output.\n";
		close LOG;
		exit;
	}
}
## END CYTOSCAPE PRINTOUT ##

# BEGIN PLACE BRANCHES #
$ts = localtime(time);
print LOG "$ts\tPlacing branches.\n";

$_ = 0 for my($in_trunk,$added,$max_link_to_trunk,$linked_position,$added_scaffold,$sub_scaffoldID);
$_ = "" for my($vertex,$max_link_vertex,$previous_vertex,$subsequent_vertex,$final_path,$node_path,$contig);
my(@TRUNK_VERTICES,@INCLUDED_VERTICES);
my(%TRUNK_POSITION,%INCLUDED_CONTIG,%SCAFFOLD_CONTIGS,%CONTIG_ASSIGNMENT,%SCAFFOLD_ORIENT,%PLACED_WITH,
   %ANCHOR_CT,%FINAL_PATHS,%PLACED_TO,%PLACED_VERTICES,%CLEARED_ANCHOR);

if ($opt{'c'} eq "D") {$opt{'c'} = $nodeCT}

open ORDERED, ">$PFX.fragScaff.ordered.txt";
foreach $scaffoldID (keys %TRUNK) {

	@TRUNK_VERTICES = split(/,/, $TRUNK{$scaffoldID});
	$in_trunk = 0;
	%PLACED_WITH = (); %ANCHOR_CT = (); %PLACED_TO = (); %FINAL_PATHS = ();
	
	# PLACE BRANCHES INTO TRUNK
	# -c option in development
	while ($in_trunk < $VERTEX_CT{$scaffoldID}) {
		if ($in_trunk > 0) {
			@TRUNK_VERTICES = ();
			foreach $vertex (sort {$TRUNK_POSITION{$a}<=>$TRUNK_POSITION{$b}} keys %TRUNK_POSITION) {
				push @TRUNK_VERTICES, $vertex;
			}
		}
		$in_trunk = @TRUNK_VERTICES;
		%TRUNK_POSITION = ();
		for ($i = 0; $i < @TRUNK_VERTICES; $i++) {
			$TRUNK_POSITION{$TRUNK_VERTICES[$i]} = $i;
		}
		
		if ($in_trunk < $VERTEX_CT{$scaffoldID}) {
			$added = 0;
			foreach $vertex1 (keys %{$SUBSET_VERTICES{$scaffoldID}}) {
				if (!defined $TRUNK_POSITION{$vertex1} && $added == 0) {
					$max_link_to_trunk = 0;
					$max_link_vertex = "";
					foreach $vertex2 (keys %{$EDGES{$vertex1}}) {
						if (defined $TRUNK_POSITION{$vertex2}) {
							if ($EDGES{$vertex1}{$vertex2} > $max_link_to_trunk) {
								$max_link_to_trunk = $EDGES{$vertex1}{$vertex2};
								$max_link_vertex = $vertex2;
							}
						}
					}
					if ($max_link_to_trunk > 0) {
						$linked_position = $TRUNK_POSITION{$max_link_vertex};
						$previous_vertex = $TRUNK_VERTICES[$linked_position-1];
						$subsequent_vertex = $TRUNK_VERTICES[$linked_position+1];
						
						$vertex2 = $max_link_vertex;
						while (defined $PLACED_WITH{$vertex2}) {$vertex2 = $PLACED_WITH{$vertex2}}
						$PLACED_WITH{$vertex1} = $vertex2;
						$PLACED_TO{$vertex2}{$vertex1} = 1;
						$ANCHOR_CT{$vertex2}++;
						
						if ($EDGES{$max_link_vertex}{$previous_vertex} < $EDGES{$max_link_vertex}{$subsequent_vertex}) {
							$TRUNK_POSITION{$vertex1} = $linked_position-0.5;
						} else {
							$TRUNK_POSITION{$vertex1} = $linked_position+0.5;
						}
						$added++;
					}
				}
			}
		}
	}
	@TRUNK_VERTICES = ();
	foreach $vertex (sort {$TRUNK_POSITION{$a}<=>$TRUNK_POSITION{$b}} keys %TRUNK_POSITION) {
		push @TRUNK_VERTICES, $vertex;
	}
	%TRUNK_POSITION = ();
	for ($i = 0; $i < @TRUNK_VERTICES; $i++) {$TRUNK_POSITION{$TRUNK_VERTICES[$i]} = $i}
	
	@INCLUDED_VERTICES = ();
	%PLACED_VERTICES = ();
	%CLEARED_ANCHOR = ();
	$added_scaffold = 0;
	for ($i = 0; $i < @TRUNK_VERTICES; $i++) {
		$vertex = $TRUNK_VERTICES[$i];
		if (!defined $PLACED_VERTICES{$vertex}) {
			
			if (defined $PLACED_WITH{$vertex} && !defined $CLEARED_ANCHOR{$PLACED_WITH{$vertex}} && $ANCHOR_CT{$PLACED_WITH{$vertex}} >= $opt{'c'}) {
				if (@INCLUDED_VERTICES > 1) {
					if ($added_scaffold < 0.5) {
						@{$FINAL_PATHS{$scaffoldID}} = @INCLUDED_VERTICES;
						$added_scaffold = 1;
					} else {
						@{$FINAL_PATHS{$scaffold_ct}} = @INCLUDED_VERTICES;
						$scaffold_ct++;
					}
				}
				@INCLUDED_VERTICES = ();
				$CLEARED_ANCHOR{$PLACED_WITH{$vertex}} = 1;
			}
			
			push @INCLUDED_VERTICES, $vertex;
			$PLACED_VERTICES{$vertex} = 1;
			
			if (defined $ANCHOR_CT{$vertex} && !defined $CLEARED_ANCHOR{$vertex}) {
				if ($ANCHOR_CT{$vertex} >= $opt{'c'}) {
					
					$vertex2 = $TRUNK_VERTICES[$i+1];
					$vertex3 = $vertex2;
					if ($vertex2 =~ /L$/) {$vertex2 =~ s/L/R/} else {$vertex2 =~ s/R/L/}
					if ($vertex eq $vertex2) {
						push @INCLUDED_VERTICES, $vertex3;
						$PLACED_VERTICES{$vertex3} = 1;
					}
					
					
					for ($j = 0; $j < @TRUNK_VERTICES; $j++) {
						if (defined $PLACED_TO{$vertex}{$TRUNK_VERTICES[$j]}) {
							if (!defined $PLACED_VERTICES{$TRUNK_VERTICES[$j]}) {
								push @INCLUDED_VERTICES, $TRUNK_VERTICES[$j];
								$PLACED_VERTICES{$TRUNK_VERTICES[$j]} = 1;
							} else {
							}
						}
					}
					if ($added_scaffold < 0.5) {
						@{$FINAL_PATHS{$scaffoldID}} = @INCLUDED_VERTICES;
						$added_scaffold = 1;
					} else {
						@{$FINAL_PATHS{$scaffold_ct}} = @INCLUDED_VERTICES;
						$scaffold_ct++;
					}
					@INCLUDED_VERTICES = ();
				}
			}
		}
	}
	if ($added_scaffold < 0.5) {
		@{$FINAL_PATHS{$scaffoldID}} = @INCLUDED_VERTICES;
	} else {
		@{$FINAL_PATHS{$scaffold_ct}} = @INCLUDED_VERTICES;
		$scaffold_ct++;
	}
	
	# PRINT OUT TRUNK(S)
	%INCLUDED_CONTIG = ();
	foreach $sub_scaffoldID (sort {$a<=>$b} keys %FINAL_PATHS) {
		$final_path = "";
		$node_path = "";
		for ($i = 0; $i < @{$FINAL_PATHS{$sub_scaffoldID}}; $i++) {
			$node_path .= "$FINAL_PATHS{$sub_scaffoldID}[$i],";
			($contigID,$side) = split(/\./, $FINAL_PATHS{$sub_scaffoldID}[$i]);
			if (!defined $INCLUDED_CONTIG{$contigID}) {
				$INCLUDED_CONTIG{$contigID} = 1;
				push @{$SCAFFOLD_CONTIGS{$sub_scaffoldID}}, $CONTIG_NAME[$contigID];
				$CONTIG_ASSIGNMENT{$CONTIG_NAME[$contigID]} = $sub_scaffoldID;
				if ($side =~ /L/) {
					$final_path .= "$contigID.f,";
					push @{$SCAFFOLD_ORIENT{$sub_scaffoldID}}, "f";
				} else {
					$final_path .= "$contigID.r,";
					push @{$SCAFFOLD_ORIENT{$sub_scaffoldID}}, "r";
				}
			}
		} $final_path =~ s/,$//; $node_path =~ s/,$//;
		print ORDERED "$sub_scaffoldID\t$final_path\t$node_path\n";
	}
}
close ORDERED;
## END PLACE BRANCHES ##

# BEGIN PERFORM ESTIMATED N50 IMPROVEMENT #
$ts = localtime(time);
print LOG "$ts\tGenerating estimated contiguity improvement metrics.\n";
$_ = 0 for my($estN50_inputTotal,$estN50_runningTotal,$estN50_inputN10,$estN50_outputN10,$estN50_inputN50,
              $estN50_outputN50,$estN50_inputN90,$estN50_outputN90,$estN50_imp10,$estN50_imp50,$estN50_imp90,
			  $estN50_scaffoldTotal,$estN50_includedTotal,$estN50_includedFrac);
my(@ESTN50_SCAFF_SIZES);
my(%ESTN50_INCLUDED_CONTIGS);
foreach $estN50_scaffoldTotal (@CONTIG_LENGTH) {$estN50_inputTotal+=$estN50_scaffoldTotal}
foreach $contig_name (sort {$CONTIG_LENGTH[$CONTIG_ID{$b}]<=>$CONTIG_LENGTH[$CONTIG_ID{$a}]} keys %CONTIG_ID) {
	$contigID = $CONTIG_ID{$contig_name};
	$estN50_runningTotal+=$CONTIG_LENGTH[$contigID];
	if ($estN50_runningTotal>=($estN50_inputTotal*0.1)&&$estN50_inputN10<1) {
		$estN50_inputN10=$CONTIG_LENGTH[$contigID];
	}
	if ($estN50_runningTotal>=($estN50_inputTotal*0.5)&&$estN50_inputN50<1) {
		$estN50_inputN50=$CONTIG_LENGTH[$contigID];
	}
	if ($estN50_runningTotal>=($estN50_inputTotal*0.9)&&$estN50_inputN90<1) {
		$estN50_inputN90=$CONTIG_LENGTH[$contigID];
	}
}
foreach $scaffoldID (keys %SCAFFOLD_CONTIGS) {
	$estN50_scaffoldTotal = 0;
	for ($i = 0; $i < @{$SCAFFOLD_CONTIGS{$scaffoldID}}; $i++) {
		$estN50_scaffoldTotal+=$CONTIG_LENGTH[$CONTIG_ID{$SCAFFOLD_CONTIGS{$scaffoldID}[$i]}];
		$ESTN50_INCLUDED_CONTIGS{$CONTIG_ID{$SCAFFOLD_CONTIGS{$scaffoldID}[$i]}} = 1;
	}
	push @ESTN50_SCAFF_SIZES, $estN50_scaffoldTotal;
	$estN50_includedTotal+=$estN50_scaffoldTotal;
}
$estN50_includedFrac = sprintf("%.2f", ($estN50_includedTotal/$estN50_inputTotal)*100);
for ($i = 0; $i < $contigCT; $i++) {
	if (!defined $ESTN50_INCLUDED_CONTIGS{$i}) {
		push @ESTN50_SCAFF_SIZES, $CONTIG_LENGTH[$i];
	}
}
$estN50_runningTotal = 0;
foreach $estN50_scaffoldTotal (sort {$b<=>$a} @ESTN50_SCAFF_SIZES) {
	$estN50_runningTotal+=$estN50_scaffoldTotal;
	if ($estN50_runningTotal>=($estN50_inputTotal*0.1)&&$estN50_outputN10<1) {
		$estN50_outputN10=$estN50_scaffoldTotal;
	}
	if ($estN50_runningTotal>=($estN50_inputTotal*0.5)&&$estN50_outputN50<1) {
		$estN50_outputN50=$estN50_scaffoldTotal;
	}
	if ($estN50_runningTotal>=($estN50_inputTotal*0.9)&&$estN50_outputN90<1) {
		$estN50_outputN90=$estN50_scaffoldTotal;
	}
}

$estN50_imp10 = sprintf("%.2f", $estN50_outputN10/$estN50_inputN10);
$estN50_imp50 = sprintf("%.2f", $estN50_outputN50/$estN50_inputN50);
$estN50_imp90 = sprintf("%.2f", $estN50_outputN90/$estN50_inputN90);

print LOG "\tIncl:\t$estN50_includedTotal\tof\t$estN50_inputTotal\t\($estN50_includedFrac\)
\tN10:\t$estN50_inputN10\t\>\t$estN50_outputN10\t\($estN50_imp10\)
\tN50:\t$estN50_inputN50\t\>\t$estN50_outputN50\t\($estN50_imp50\)
\tN90:\t$estN50_inputN90\t\>\t$estN50_outputN90\t\($estN50_imp90\)\n";

@ESTN50_SCAFF_SIZES = ();
%ESTN50_INCLUDED_CONTIGS = ();
## END PERFORM ESTIMATED N50 IMPROVEMENT ##

if (defined $opt{'I'}) {
	$ts = localtime(time);
	print LOG "$ts\tOption I toggled, exiting after graph manipulations.\n";
	close LOG;
	exit;
}

## END PRUNE BY SIZE ##

# BEGIN REVCOMP SUB #
my($fwd_seq,$revcomp_seq) = ("","");
my$iii = 0;
my@FS;
my%REVCOMP_CONVERT;
sub rev_comp {
	%REVCOMP_CONVERT = ("A" => "T", "T" => "A", "C" => "G", "G" => "C", "N" => "N");
	$fwd_seq = $_[0];
	$revcomp_seq = "";
	@FS = split(//, $fwd_seq);
	for ($iii = (@FS-1); $iii >= 0; $iii--) {
		$revcomp_seq .= $REVCOMP_CONVERT{$FS[$iii]};
	}
	return $revcomp_seq;
}
## END REVCOMP SUB ##

# BEGIN READ IN FASTA #
$ts = localtime(time);
print LOG "$ts\tLoading input assembly fasta.\n";

my($total_size,$total_nonN,$scaff_nonN) = (0,0,0);
my($sequence,$base) = ("","");
my(%ORPHAN,%SCAFF_SIZES_nonN,%SCAFF_SIZES,%CONTIGS);
my@S;

$contig = "";
if ($opt{'F'} =~ /\.gz$/) {
	open IN, "zcat $opt{'F'} |";
} else {
	open IN, "$opt{'F'}";
}
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^>/) {
		if ($contig ne "") {
			$CONTIGS{$contig} = $sequence;
			if (!defined $CONTIG_ASSIGNMENT{$contig}) {
				$ORPHAN{$contig} = 1;
			}
			$SCAFF_SIZES{$contig} = length($sequence);
			$scaff_nonN = 0;
			@S = split(//, $sequence);
			foreach $base (@S) {if ($base !~ /N/i) {$scaff_nonN++}}
			$SCAFF_SIZES_nonN{$contig} = $scaff_nonN;
			$total_size += $SCAFF_SIZES{$contig};
			$total_nonN += $scaff_nonN;
		}
		$contig = $l;
		$contig =~ s/^>//;
		$sequence = "";
	} else {
		$sequence .= $l;
	}
} close IN;
## END READ IN FASTA ##

# BEGIN CALC CURRENT N50 STATS #
$_ = 0 for my($max_size,$max_size_nonN,$N10,$N10_nonN,$N50,$N50_nonN,$N90,$N90_nonN,$current_size);
my$scaffold = "";
sub fastaN50 {
	$_ = 0 for ($max_size,$max_size_nonN,$N10,$N10_nonN,$N50,$N50_nonN,$N90,$N90_nonN,$current_size);
	foreach $scaffold (sort {$SCAFF_SIZES{$b}<=>$SCAFF_SIZES{$a}} keys %SCAFF_SIZES) {
		$current_size += $SCAFF_SIZES{$scaffold};
		if ($max_size == 0) {
			$max_size = $SCAFF_SIZES{$scaffold};
		}
		if ($N10 == 0 && $current_size >= (0.1*$total_size)) {
			$N10 = $SCAFF_SIZES{$scaffold};
		}
		if ($N50 == 0 && $current_size >= (0.5*$total_size)) {
			$N50 = $SCAFF_SIZES{$scaffold};
		}
		if ($N90 == 0 && $current_size >= (0.9*$total_size)) {
			$N90 = $SCAFF_SIZES{$scaffold};
		}
	}
	$current_size = 0;
	foreach $scaffold (sort {$SCAFF_SIZES_nonN{$b}<=>$SCAFF_SIZES_nonN{$a}} keys %SCAFF_SIZES_nonN) {
		$current_size += $SCAFF_SIZES_nonN{$scaffold};
		if ($max_size_nonN == 0) {
			$max_size_nonN = $SCAFF_SIZES_nonN{$scaffold};
		}
		if ($N10_nonN == 0 && $current_size >= (0.1*$total_nonN)) {
			$N10_nonN = $SCAFF_SIZES_nonN{$scaffold};
		}
		if ($N50_nonN == 0 && $current_size >= (0.5*$total_nonN)) {
			$N50_nonN = $SCAFF_SIZES_nonN{$scaffold};
		}
		if ($N90_nonN == 0 && $current_size >= (0.9*$total_nonN)) {
			$N90_nonN = $SCAFF_SIZES_nonN{$scaffold};
		}
	}
}

fastaN50();

open N50OUT, ">$PFX.fragScaff.N50.txt";

print N50OUT "
PRE FRAGSCAFF ASSEMBLY:

	ASSEMBLY STATS (With N's):
	MAX SCAFFOLD = $max_size
	N10          = $N10
	N50          = $N50
	N90          = $N90

	ASSEMBLY STATS (no N's):
	MAX SCAFFOLD = $max_size_nonN
	N10          = $N10_nonN
	N50          = $N50_nonN
	N90          = $N90_nonN

";
## END CALC CURRENT N50 STATS ##

# BEGIN MAKE NEW ASSEMBLY #

$ts = localtime(time);
print LOG "$ts\tPrinting out fragScaff assembly.\n";

($total_size,$total_nonN) = (0,0);
($sequence,$base) = ("","");
(%SCAFF_SIZES_nonN,%SCAFF_SIZES) = ();
@S = ();

$_ = 0 for my($scaff_size,$scaff_size_nonN,$row_ct,$qual_scaffold,$qual_tally,$orient_quality);
$_ = "" for my($scaffold_out,$qualities_out,$qual_contig,$qual_return);

my$Nspacer = ""; my$Nbase = $opt{'n'}; my$score = $opt{'p'}*2; my$numNs = 1.5*$opt{'x'}; my$Nct = 0;

my$slope = ($opt{'x'}-$opt{'X'})/((2*$opt{'M'})-(2*$opt{'p'}));
my$intercept = (-1*($slope*(2*$opt{'p'}))) + $opt{'X'};

sub make_spacer {
	if (defined $EDGES{"$_[0].L"}{"$_[1].L"}) { #LL
		$score = $EDGES{"$_[0].L"}{"$_[1].L"};
		$numNs = ($slope*$score) + $intercept;
	} elsif (defined $EDGES{"$_[0].R"}{"$_[1].L"}) { #RL
		$score = $EDGES{"$_[0].R"}{"$_[1].L"};
		$numNs = ($slope*$score) + $intercept;
	} elsif (defined $EDGES{"$_[0].R"}{"$_[1].R"}) { #RR
		$score = $EDGES{"$_[0].R"}{"$_[1].R"};
		$numNs = ($slope*$score) + $intercept;
	} else {
		$numNs = 1.5*$opt{'x'};
	}
	$Nspacer = "";
	for ($Nct = 0; $Nct < $numNs; $Nct++) {
		$Nspacer .= $Nbase;
	}
	return $Nspacer;
}

open OUT, ">$PFX.fragScaff.assembly.fasta";
if (!defined $opt{'f'}) {
	open QUAL, ">$PFX.fragScaff.quals";
}


my$finalScaffoldID = 0; my$Nstretch = "";
foreach $scaffold (keys %SCAFFOLD_CONTIGS) {
	$scaffold_out = ""; $qualities_out = "";
	for ($i = 0; $i < @{$SCAFFOLD_CONTIGS{$scaffold}}; $i++) {
		
		if (!defined $CONTIGS{$SCAFFOLD_CONTIGS{$scaffold}[$i]}) {print STDERR "\nNO SEQUENCE FOR $SCAFFOLD_CONTIGS{$scaffold}[$i]!\n"};
		if ($CONTIGS{$SCAFFOLD_CONTIGS{$scaffold}[$i]} eq "") {print STDERR "\nEMPTY SEQUENCE FOR $SCAFFOLD_CONTIGS{$scaffold}[$i]!\n"};
		
		if ($i > 0) {
			$Nstretch = make_spacer($CONTIG_ID{$SCAFFOLD_CONTIGS{$scaffold}[$i-1]},$CONTIG_ID{$SCAFFOLD_CONTIGS{$scaffold}[$i]});
			$scaffold_out .= $Nstretch;
		}
		
		$qualities_out .= length($scaffold_out);
		$qualities_out .= score_quals($scaffold,$SCAFFOLD_CONTIGS{$scaffold}[$i]);
		
		if ($SCAFFOLD_ORIENT{$scaffold}[$i] eq "r") {
			$scaffold_out .= rev_comp($CONTIGS{$SCAFFOLD_CONTIGS{$scaffold}[$i]});
		} else {
			$scaffold_out .= $CONTIGS{$SCAFFOLD_CONTIGS{$scaffold}[$i]};
		}
	}
	if (defined $opt{'f'}) {
		print OUT ">fragScaff_scaffold_$finalScaffoldID\tFSCF_ID=$scaffold\t$qualities_out\n";
	} else {
		print OUT ">fragScaff_scaffold_$finalScaffoldID\n";
		print QUAL "fragScaff_scaffold_$finalScaffoldID\tFSCF_ID=$scaffold\t$qualities_out\n";
	}
	$scaff_size_nonN = 0;
	@S = split(//, $scaffold_out);
	$scaff_size = @S;
	$row_ct = 0;
	for ($i = 0; $i < @S; $i++) {
		print OUT "$S[$i]";
		if ($S[$i] !~ /N/i) {$scaff_size_nonN++}
		$row_ct++;
		if ($row_ct >= $opt{'L'}) {print OUT "\n"; $row_ct = 0}
	} print OUT "\n";
	$SCAFF_SIZES{$finalScaffoldID} = $scaff_size;
	$SCAFF_SIZES_nonN{$finalScaffoldID} = $scaff_size_nonN;
	$total_size += $scaff_size;
	$total_nonN += $scaff_size_nonN;
	$finalScaffoldID++;
}

sub score_quals {
	$qual_scaffold = $_[0];
	$qual_contig = $CONTIG_ID{$_[1]};
	$qual_return = "";
	if ($SCAFFOLD_ORIENT{$qual_scaffold}[$i] eq "f") {
		$vertex1 = "$qual_contig.L"; $vertex2 = "$qual_contig.R";
	} else {
		$vertex1 = "$qual_contig.R"; $vertex2 = "$qual_contig.L";
	}
	$qual_tally = 0;
	foreach $vertex3 (keys %{$SUBSET_EDGES{$qual_scaffold}{$vertex1}}) {
		if ($vertex3 ne $vertex2) {
			$qual_tally += $SUBSET_EDGES{$qual_scaffold}{$vertex1}{$vertex3};
		}
	}
	$qual_return .= "(".sprintf("%.2f", $qual_tally).":";
	
	if ($CONTIG_LENGTH[$qual_contig] <= $opt{'E'}) {
		$orient_quality = sprintf("%.2f", 0);
	} else {
		$orient_quality = sprintf("%.2f", (1-$PAIR_FRAC{$vertex1}));
	}
	$qual_return .= "$orient_quality:";
	
	$qual_tally = 0;
	foreach $vertex3 (keys %{$SUBSET_EDGES{$qual_scaffold}{$vertex2}}) {
		if ($vertex3 ne $vertex1) {
			$qual_tally += $SUBSET_EDGES{$qual_scaffold}{$vertex2}{$vertex3};
		}
	}
	$qual_return .= sprintf("%.2f", $qual_tally).")";
	
	return $qual_return;
}

foreach $contig (keys %ORPHAN) {
	if (defined $opt{'f'}) {
		print OUT ">original_scaffold_$finalScaffoldID\tINPUT_ID=$contig\t0(-1:-1:-1)\n";
	} else {
		print OUT ">original_scaffold_$finalScaffoldID\n";
		print QUAL "original_scaffold_$finalScaffoldID\tINPUT_ID=$contig\t0(-1:-1:-1)\n";
	}
	@S = split(//, $CONTIGS{$contig});
	$scaff_size = @S;
	$scaff_size_nonN = 0;
	$row_ct = 0;
	for ($i = 0; $i < @S; $i++) {
		print OUT "$S[$i]";
		if ($S[$i] !~ /N/i) {$scaff_size_nonN++}
		$row_ct++;
		if ($row_ct >= $opt{'L'}) {print OUT "\n"; $row_ct = 0}
	} print OUT "\n";
	$SCAFF_SIZES{$finalScaffoldID} = $scaff_size;
	$SCAFF_SIZES_nonN{$finalScaffoldID} = $scaff_size_nonN;
	$total_size += $scaff_size;
	$total_nonN += $scaff_size_nonN;
	$finalScaffoldID++;
}
close OUT;

if (!defined $opt{'f'}) {close QUAL}

fastaN50();

print N50OUT "
POST FRAGSCAFF ASSEMBLY:

	ASSEMBLY STATS (With N's):
	MAX SCAFFOLD = $max_size
	N10          = $N10
	N50          = $N50
	N90          = $N90

	ASSEMBLY STATS (no N's):
	MAX SCAFFOLD = $max_size_nonN
	N10          = $N10_nonN
	N50          = $N50_nonN
	N90          = $N90_nonN

";
close N50OUT;

## END MAKE NEW ASSEMBLY ##

$ts = localtime(time);
print LOG "$ts\tfragScaff Complete. Final Scaffold Count: $finalScaffoldID\n";
close LOG;

# BEGIN RUN QSUB #
sub run_qsub {
	my($qsubNode1,$qsubNode2);
	($dPFX,$nodeCT,$qsubNode1,$qsubNode2) = split(/,/, $opt{'Q'});
	pvalue_links($qsubNode1,$qsubNode2);
	system("date > $dPFX/node.$qsubNode1.$qsubNode2.complete");
}
## END RUN QSUB ##

