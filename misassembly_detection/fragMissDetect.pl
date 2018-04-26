#!/usr/bin/perl

use Getopt::Std; %opt = ();
getopts("d:B:s:W:O:G:f:F:", \%opt);

$die = "

======================== fragMissDetect =======================
                         Version 131202
            Andrew Adey, Shendure Lab (acadey\@uw.edu)

Calls fragments from bam file aligned to reference.

--- GENERAL OPTIONS -------------------------------------------
    -B   [FILE]   input bam file (must have RG tags) (REQ)
    -W   [FILE]   window bed file (REQ)
    -F   [FILE]   fasta window file was made from (REQ)
    -O   [STR]    output prefix (def = bam file name)
    -G   [R/N/H]  read group identifier:
                      R = RG:Z:group bam tag
                      N = read name (group:read_number) (def)
                      H = after hash (name#group)
    -q   [INT]    min mapping qual (def = 10)
    -f   [FLOAT]  fraction of windos to cut (def = 0.01)

--- DEPENDENCY OPTIONS ----------------------------------------
    -s   [STR]    samtools call (def = 'samtools')

";

# general defaults
if (!defined $opt{'B'}) {die $die}
if (!defined $opt{'W'}) {die $die}
if (!defined $opt{'F'}) {die $die}
if (!defined $opt{'O'}) {$opt{'O'} = $opt{'B'}};
if (!defined $opt{'s'}) {$opt{'s'} = "samtools"};
if (!defined $opt{'q'}) {$opt{'q'} = 10}
if (!defined $opt{'f'}) {$opt{'f'} = 0.01}
if (!defined $opt{'G'}) {$opt{'G'} = "N"} elsif ($opt{'G'} !~ /[RNH]/) {
	die "
    Can not determine group option (-G), please select from:
           R = RG:Z:group bam tag
           N = read name (group:read_number) (def)
           H = after hash (name#group)
";}

$PFX = $opt{'O'};

open LOG, ">$PFX.fragMissDetect.log";

select((select(LOG), $|=1)[0]); # make LOG hot

$ts = localtime(time);
print LOG "$ts\tProgram Called:\n";
foreach $option (keys %opt) {
	print LOG "\t-$option\t$opt{$option}\n";
}

$ts = localtime(time);
print LOG "$ts\tReading in bam header to identify read groups ... ";

open HEAD, "$opt{'s'} view -H $opt{'B'} |";
$groupID = 0;
while ($l = <HEAD>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] eq "\@RG") {
		$barc = $P[1]; $barc =~ s/^ID://;
		$lib = $P[2]; $lib =~ s/^LB://;
		$GROUP_NAME[$groupID] = $barc;
		$GROUP_ID{$barc} = $groupID;
		$LIB_ID{$lib} = $groupID;
		$prev_pos[$groupID] = -1;
		$prev_chr[$groupID] = "";;
		$groupID++;
	}
} close HEAD;
$groupCT = $groupID;

print LOG "done.\n";
$ts = localtime(time);
if ($groupCT < 1) {
	print LOG "\n!!!! ERROR !!!! NO GROUPS FOUND !!!! ERROR !!!!\n";
	die "\n!!!! ERROR !!!! NO GROUPS FOUND !!!! ERROR !!!!\n";
} else {
	print LOG "$ts\t$groupCT groups found.\n";
}

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
	}
}

$ts = localtime(time);
print LOG "$ts\tLoading in windows ... ";

$winID = 0;
$prevSeq = "";
$totalWin = 0;
open IN, "$opt{'W'}";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] ne $prevSeq) {$winID = 0}
	$WIN_START{$P[0]}[$winID] = $P[1];
	$WIN_END{$P[0]}[$winID] = $P[2];
	$WIN_COUNT{$P[0]}++;
	$winID++;
	$prevSeq = $P[0];
	$totalWin++;
} close IN;

$ts = localtime(time);
print LOG "done.\n$ts\tParsing reads ...\n";

$check = 0.05; $increment = 0.05; $parsed = 0; $winID = 0; $exclID = 0;
open BAM, "$opt{'s'} view -q $opt{'q'} $opt{'B'} |";
open OUT, ">$PFX.fragMissDetect.winHits.txt";
$chrom = "";
while ($l = <BAM>) {
	chomp $l;
	@P = split(/\t/, $l);
	get_rg();
	if ($group_found > 0.5) {
		if ($P[2] ne $chrom) {
			$winID = 0;
			$chrom = $P[2];
			$hitCount = 0;
			$hitList = "";
		}
		# CHECK THAT READ IS BEFORE THE END OF THE CURRENT WINDOW - OTHERWISE PARSE THROUGH
		while (defined $WIN_END{$chrom}[$winID] && $WIN_END{$chrom}[$winID] <= $P[3]) {
			# new read is after the end of the window
			# finish window
			$hitList =~ s/,$//;
			print OUT "$chrom\t$WIN_START{$chrom}[$winID]\t$WIN_END{$chrom}[$winID]\t$winID\t$hitCount\t$hitList\n";
			$HIT_COUNT{$chrom}[$winID] = $hitCount;
			%{$HIT_LIST{$chrom}[$winID]} = %HITS;
			$win_tally++;
			# progress check
			$frac_complete = sprintf("%.2f", $win_tally/$totalWin);
			if ($frac_complete >= $check) {
				$ts = localtime(time);
				print LOG "\t$ts\t$check frac complete.\n";
				$check += $increment;
			}
			# start new window
			$winID++;
			$hitList = "";
			$hitCount = 0;
			%HITS = ();
		}
		
		if ($prev_loc{$groupID} ne "$P[3]:$P[4]") {
			# read is before the end of the current window
			if ($P[3] > $WIN_START{$chrom}[$winID]) {
				# read is after the start of the current window & before the start of the current exclusion window & not a dup
				if (!defined $HITS{$groupID}) {
					$HITS{$groupID} = 1;
					$hitList .= "$groupID,";
					$hitCount++;
				}
			}
		}
		$prev_loc{$groupID} = "$P[3]:$P[4]";
	}
}
# print last window
$hitList =~ s/,$//;
print OUT "$chrom\t$WIN_START{$chrom}[$winID]\t$WIN_END{$chrom}[$winID]\t$winID\t$hitCount\t$hitList\n";
%{$HIT_LIST{$chrom}[$winID]} = %HITS;
$HIT_COUNT{$chrom}[$winID] = $hitCount;
close OUT;
close BAM;

$ts = localtime(time);
print LOG "$ts\tCalculating scores ... ";
open OUT, ">$PFX.fragMissDetect.shareScores.txt";

foreach $chrom (keys %HIT_COUNT) {
	if (@{$HIT_COUNT{$chrom}} >= 5) {
		for ($winID = 0; $winID < (@{$HIT_COUNT{$chrom}}-2); $winID++) {
			if ($HIT_COUNT{$chrom}[$winID]>10) {
				if ($HIT_COUNT{$chrom}[$winID+1]>10) {
					$shared1 = 0;
					foreach $groupID (keys %{$HIT_LIST{$chrom}[$winID]}) {
						if (defined $HIT_LIST{$chrom}[$winID+1]{$groupID}) {
							$shared1++;
						}
					}
					if ((($HIT_COUNT{$chrom}[$winID]+$HIT_COUNT{$chrom}[$winID+1])-$shared1)>0) {
						$frac1 = sprintf("%.2f", $shared1/(($HIT_COUNT{$chrom}[$winID]+$HIT_COUNT{$chrom}[$winID+1])-$shared1));
					} else {
						$frac1 = sprintf("%.2f", 1);
					}
					$FRAC_1_HIST{$frac1}++;
					$frac1_total++;
					$WIN_FRAC1{$chrom}[$winID] = $frac1;
				} else {
					$frac1 = -1;
				}
				if ($HIT_COUNT{$chrom}[$winID+2]>10) {
					$shared2 = 0;
					foreach $groupID (keys %{$HIT_LIST{$chrom}[$winID]}) {
						if (defined $HIT_LIST{$chrom}[$winID+2]{$groupID}) {
							$shared2++;
						}
					}
					if ((($HIT_COUNT{$chrom}[$winID]+$HIT_COUNT{$chrom}[$winID+2])-$shared2)>0) {
						$frac2 = sprintf("%.2f", $shared2/(($HIT_COUNT{$chrom}[$winID]+$HIT_COUNT{$chrom}[$winID+2])-$shared2));
					} else {
						$frac2 = sprintf("%.2f", 1);
					}
					$FRAC_2_HIST{$frac2}++;
					$frac2_total++;
					$WIN_FRAC2{$chrom}[$winID] = $frac2;
				} else {
					$frac2 = -1;
				}
				print OUT "$chrom\t$winID\t$HIT_COUNT{$chrom}[$winID]\t$frac1\t$frac2\n";
			} else {
				print OUT "$chrom\t$winID\t$HIT_COUNT{$chrom}[$winID]\t-1\t-1\n";
			}
		}
	} else {
		print OUT "$chrom\t-1\t-1\t-1\t-1\n";
	}
} close OUT;

%$HIT_LIST = (); %HIT_COUNT = ();

$ts = localtime(time);
print LOG "done.\n$ts\tPrinting histogram & determining fraction cuts ... ";

$frac1_count = 0; $frac2_count = 0; $frac1_cut = -1; $frac2_cut = -1;
open OUT, ">$PFX.fragMissDetect.sharedHist.txt";
for ($i = 0; $i <= 1; $i+=0.01) {
	$bin = sprintf("%.2f", $i);
	$outLine = "$bin\t";
	if (defined $FRAC_1_HIST{$bin}) {
		$outLine .= "$FRAC_1_HIST{$bin}\t";
		$frac1_count+=$FRAC_1_HIST{$bin};
		if (($frac1_count/$frac1_total)>=$opt{'f'}&&$frac1_cut<0) {
			$frac1_cut = $bin;
		}
	} else {
		$outLine .= "0\t";
	}
	if (defined $FRAC_2_HIST{$bin}) {
		$outLine .= "$FRAC_2_HIST{$bin}";
		$frac2_count+=$FRAC_2_HIST{$bin};
		if (($frac2_count/$frac2_total)>=$opt{'f'}&&$frac2_cut<0) {
			$frac2_cut = $bin;
		}
	} else {
		$outLine .= "0";
	}
	print OUT "$outLine\n";
} close OUT;

$ts = localtime(time);
print LOG "adjacent cut = $frac1_cut, skip 1 cut = $frac2_cut\n$ts\tDetermining windows to cut ... ";

open OUT, ">$PFX.fragMissDetect.windowCuts.txt";
foreach $chrom (keys %WIN_COUNT) {
	for ($winID = 0; $winID < $WIN_COUNT{$chrom}; $winID++) {
		if (defined $WIN_FRAC1{$chrom}[$winID] && defined $WIN_FRAC2{$chrom}[$winID]) {
			if ($WIN_FRAC1{$chrom}[$winID]<=$frac1_cut&&$WIN_FRAC2{$chrom}[$winID]<=$frac2_cut) {
				$fail_count++;
				print OUT "$chrom\t$WIN_START{$chrom}[$winID+1]\t$WIN_END{$chrom}[$winID+1]\t$WIN_FRAC1{$chrom}[$winID]\t$WIN_FRAC2{$chrom}[$winID]\n";
				$CUT_WIN{$chrom}{$WIN_START{$chrom}[$winID+1]} = $WIN_END{$chrom}[$winID+1];
			} else {
				$pass_count++;
			}
		}
	}
} close OUT;

$ts = localtime(time);
print LOG "pass count = $pass_count, fail count = $fail_count\n$ts\tCutting assembly ... ";











