#!/usr/bin/perl

$die = "
Usage:

fragMissDetect_window.pl <assemly_fasta.fa> <window_size>

STDOUT = window bed file

Windows will exclude all N's.

";

if (!defined $ARGV[1]) {die $die}

$winSize = $ARGV[1];

$currentSize = -1;
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^>/) {
		if ($currentSize > -0.5) {
			if (@WIN_SET>=4) {
				for ($winID = 0; $winID < (@WIN_SET-1); $winID++) {
					print "$seqName\t$WIN_SET[$winID]\t$WIN_SET[$winID+1]\n";
				}
			}
		}
		$seqName = $l; $seqName =~ s/\s.*$//; $seqName =~ s/^>//;
		@WIN_SET = ();
		$currentSize = 0;
		$position = 0;
		$winStarted = 0;
	} else {
		$posLength = length($l);
		$mapSeq = $l; $mapSeq =~ s/N//ig;
		$mapLength = length($mapSeq);
		if ($winStarted < 0.5) {
			if ($mapLength>0) {
				if (($mapLength == $posLength)||($posLength !~ /^N/i)) {
					push @WIN_SET, 0;
				} else {
					@S = split(//, $l);
					for ($basePos = 0; $basePos < @S; $basePos++) {
						if ($S[$basePos] !~ /N/i) {
							push @WIN_SET, $basePos;
							$basePos = @S+1;
						}
					}
				}
				$winStarted = 1;
			}
		}
		if (($currentSize + $mapLength)>=$winSize) {
			$basesToAdd = ($currentSize + $mapLength) - $winSize;
			@S = split(//, $l);
			$basesAdded = 0;
			for ($basePos = 0; $basePos < @S; $basePos++) {
				if ($S[$basePos] !~ /N/i) {
					$basesAdded++;
					if ($basesAdded >= $basesToAdd) {
						$winPos = $basePos+$position;
						push @WIN_SET, $winPos;
						$basePos = @S+1;
					}
				}
			}
			$currentSize = $mapLength-$basesToAdd;
		} else {
			$currentSize += $mapLength;
		}
		$position += $posLength;
	}
} close IN;
















