#!/usr/bin/perl

$die = "

ARGV0 = actual node alignments
ARGV1 = ordered contigs ouput

ARGV2 (opt) = Qualities

Auto generates output.

";

if (!defined $ARGV[1]) {die $die};

## query_ID      best_target     start_on_target stop_on_target  unique_alignability     target_specificity
#0       15      28229271        24732220        0.950712        0.999285
#1       13      97230998        93885731        0.955655        0.999969
#2       9       4904079 1811596 0.948586        0.999874


open IN, "$ARGV[0]";
while ($l = <IN>) {
	if ($l !~ /^#/) {
		chomp $l;
		@P = split(/\t/, $l);
		$P[0] =~ s/\D//g;
		if ($P[1] > -0.5) {
			if ($ARGV[0] =~ /NA12878_Assembly\.mapping\.txt/ || $ARGV[0] =~ /trueMapping/) {
				$chr = $P[1]+1;
			} else {
				$chr = $P[1];
			}
			if ($P[1] > -0.5) {
				if ($P[2]<$P[3]) {
					$CONTIG_CHR{$P[0]} = $chr;
					$CONTIG_START{$P[0]} = $P[2];
					$CONTIG_END{$P[0]} = $P[3];
					if (!defined $ORDERING{$chr}{$P[2]}) {
						$ORDERING{$chr}{$P[2]} = $P[0];
					} else {
						$pos = $P[2];
						while (defined $ORDERING{$chr}{$pos}) {
							$pos += 0.0001;
						}
						$ORDERING{$chr}{$pos} = $P[0];
					}
					$ORIENTATION{$P[0]} = "f";
				} else {
					$CONTIG_CHR{$P[0]} = $chr;
					$CONTIG_START{$P[0]} = $P[3];
					$CONTIG_END{$P[0]} = $P[2];
					if (!defined $ORDERING{$chr}{$P[3]}) {
						$ORDERING{$chr}{$P[3]} = $P[0];
					} else {
						$pos = $P[3];
						while (defined $ORDERING{$chr}{$pos}) {
							$pos += 0.0001;
						}
						$ORDERING{$chr}{$pos} = $P[0];
					}
					$ORIENTATION{$P[0]} = "r";
				}
			}
		}
	}
} close IN;

$order_position = 0;
foreach $chr (sort {$a<=>$b} keys %ORDERING) {
	foreach $pos (sort {$a<=>$b} keys %{$ORDERING{$chr}}) {
		$CONTIG_ORDER_TRUTH{$ORDERING{$chr}{$pos}} = $order_position;
		$order_position++;
	}
}

#2316    6091.f,8581.r   6091.L,6091.R,8581.R,8581.L
#2218    8329.f,9049.r   8329.L,8329.R,9049.R,9049.L
#1614    1444.f,11781.f  1444.L,1444.R,11781.L,11781.R


$global_orient_agree = 0;
$global_orient_total = 0;
$total_scaffolds = 0;
open IN, "$ARGV[1]";
$out = $ARGV[1]; $out =~ s/txt/check.txt/;
open OUT, ">$out";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	@CALL = split(/,/, $P[1]);
	$total_scaffolds++;
	$total_contigs_in_call = 0;
	%CHR_MAP = (); %CONTIG_ORDER = (); %CONTIG_ORIENT = ();
	$win_chr = ""; $win_chr_CT = 0; $total_chr_CT = 0;
	$same_orient_agree = 0; $switch_orient_agree = 0;
	$ordering_actual = ""; $orienting_actual = ""; $orienting_called = ""; $positions_actual = "";
	@ACTUAL_ORDERING = ();
	for ($i = 0; $i < @CALL; $i++) {
		($contig,$orientation) = split(/\./, $CALL[$i]);
		if (defined $CONTIG_ORDER_TRUTH{$contig}) {
			$ordering_actual .= "$CONTIG_ORDER_TRUTH{$contig},";
			$orienting_actual .= "$ORIENTATION{$contig},";
			$orienting_called .= "$orientation,";
			$positions_actual .= "$CONTIG_CHR{$contig}:$CONTIG_START{$contig}:$CONTIG_END{$contig},";
			if ($ORIENTATION{$contig} eq $orientation) {
				$same_orient_agree++;
			} else {
				$switch_orient_agree++;
			}
			push @ACTUAL_ORDERING, $CONTIG_ORDER_TRUTH{$contig};
			$total_contigs_in_call++;
		} else {
			$ordering_actual .= "X,";
			$orienting_actual .= "X,";
			$orienting_called .= "X,";
			$positions_actual .= "X:X:X,";
			$total_contigs_in_call++;
		}
	}
	if ($total_contigs_in_call>1) {
		$global_orient_total += $total_contigs_in_call;
		$ordering_actual =~ s/,$//; $orienting_actual =~ s/,$//; $orienting_called =~ s/,$//; $positions_actual =~ s/,$//;
		$properly_placed = 0;
		for ($i = 0; $i < @ACTUAL_ORDERING; $i++) {
			$properly_linked = 0;
			$improperly_linked = 0;
			for ($j = 0; $j < @ACTUAL_ORDERING; $j++) {
				if (abs($ACTUAL_ORDERING[$i]-$ACTUAL_ORDERING[$j])<100) {
					$properly_linked++;
				} else {
					$improperly_linked++;
				}
			}
			$properly_placed+=($properly_linked/($properly_linked+$improperly_linked));
		}
		$global_total_links += $total_contigs_in_call;
		$global_proper_links += $properly_placed;
		$frac_properly_linked = int(($properly_placed/$total_contigs_in_call)*100);
		
		if ($same_orient_agree>=$switch_orient_agree) {
			$orient_agree_frac = int(($same_orient_agree/$total_contigs_in_call)*100);
			$global_orient_agree += $same_orient_agree;
		} else {
			$orient_agree_frac = int(($switch_orient_agree/$total_contigs_in_call)*100);
			$global_orient_agree += $switch_orient_agree;
		}
		print OUT "$P[0]\t$ordering_actual\t$frac_properly_linked\t$orient_agree_frac\t$orienting_actual\t$orienting_called\t$P[1]\t$positions_actual\n";
	}

} close IN;
close OUT;

#fragScaff_scaffold_0    FSCF_ID=###	0(0.00:0.00:143.70)7826(143.70:0.00:0.00)
#fragScaff_scaffold_1    FSCF_ID=###    0(0.00:1.00:113.36)1285076(113.36:0.93:98.28)1825276(9...

if (defined $ARGV[2]) {
	$qual_all_ori = $out; $qual_all_ori =~ s/txt/qual_all_ori/;
	open OUTQO, ">$qual_all_ori";
	$qual_all = $out; $qual_all =~ s/txt/qual_all/;
	open OUTQ, ">$qual_all";
	open IN, $ARGV[2];
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$ID = $P[1]; $ID =~ s/^FSCF_ID=//;
		@QUALS = split(/[\(\)]/, $P[2]);
		for ($i = 0; $i < @QUALS; $i++) {
			$coord = $QUALS[$i];
			$i++;
			($leftQ,$oriQ,$rightQ) = split(/:/, $QUALS[$i]);
			#print STDERR "$i\t$leftQ,$oriQ,$rightQ\n";
			push @{$LEFT_QUALS{$ID}}, $leftQ;
			print OUTQ "$leftQ\n";
			push @{$ORIENT_QUALS{$ID}}, $oriQ;
			print OUTQO "$oriQ\n";
			push @{$RIGHT_QUALS{$ID}}, $rightQ;
			print OUTQ "$rightQ\n";
			
		}
	} close IN; close OUTQ, close OUTQO;
	$qual_long = $out; $qual_long =~ s/txt/qual_long/;
	open OUTQ, ">$qual_long";
	$qual_ori = $out; $qual_ori =~ s/txt/qual_ori/;
	open OUTQO, ">$qual_ori";
}
#exit;
#559     6811,6813       100     50      r,f     r,r     12577.r,11654.r,13332.r,13681.r,10602.r 7:56931906,7:56962507

$total = 0;
$correct = 0;
$correct_ori = 0;
$bases_total = 0;
$longest_perfect_size = 0;
$longest_perfect_stretch = 0;
$ori_total = 0;
$total_for_miss = 0;
open IN, "$out";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	@O = split(/,/, $P[1]);
	@O1 = split(/,/, $P[4]);
	@O2 = split(/,/, $P[5]);
	$pos = $O[0];
	$ori1 = $O1[0]; $ori2 = $O2[0];
	@TRUTH = split(/,/, $P[7]);
	%CHR_HITS = (); $misassembly_incl = 0;
	%BASE_HITS = (); $non_perfect_scaffold_check = 0;
	$total_scaff_bases = 0; $perfect_stretch = 0;
	$scaff_bases_total = 0;
	if (defined $ARGV[2]) {
		if (defined $ORIENT_QUALS{$P[0]}) {
			$length1 = @{$ORIENT_QUALS{$P[0]}}; $length2 = @O;
			if ($length1 != $length2) {
				print STDERR "LENGTH DISCREPANCY FOR $P[0]: $length1 ne $length2\n";
			}
		} else {
			print STDERR "NO QUALS FOR $P[0]!\n";
		}
	}
	for ($i = 1; $i < @O; $i++) {
		$contig_count_total++;
		if ($O1[$i] ne "X" && $O2[$i] ne "X") {
			$mapped_contig_count++;
			($chr,$aln1,$aln2) = split(/:/, $TRUTH[$i]);
			if (abs($O[$i]-$pos) <= 10) {
				$correct++;
				$perfect_stretch+=($aln2-$aln1);
			} else {
				if ($perfect_stretch>$longest_perfect_stretch) {
					$longest_perfect_stretch=$perfect_stretch;
				}
				$perfect_stretch = 0;
				$non_perfect_scaffold_check = 1;
			}
			if ($ori1 ne "X" && $ori2 ne "X") {
				$new_total++;
				if (abs($O[$i]-$pos) <= 10) {
					$new_correct++;
					$bases_prop_tally += ($aln2-$aln1);
					if (defined $ARGV[2]) {
						if (defined $LEFT_QUALS{$P[0]}) {
							$GOOD_LEFT{$LEFT_QUALS{$P[0]}[$i]}++;
							$GOOD_RIGHT{$RIGHT_QUALS{$P[0]}[$i-1]}++;
							print OUTQ "1\t$LEFT_QUALS{$P[0]}[$i]\t$RIGHT_QUALS{$P[0]}[$i-1]\n";
						}
					}
				} else {
					if (defined $ARGV[2]) {
						if (defined $LEFT_QUALS{$P[0]}) {
							$BAD_LEFT{$LEFT_QUALS{$P[0]}[$i]}++;
							$BAD_RIGHT{$RIGHT_QUALS{$P[0]}[$i-1]}++;
							print OUTQ "0\t$LEFT_QUALS{$P[0]}[$i]\t$RIGHT_QUALS{$P[0]}[$i-1]\n";
						}
					}
				}
				
				if (($ori1 eq $ori2 && $O1[$i] eq $O2[$i]) || ($ori1 ne $ori2 && $O1[$i] ne $O2[$i])) {
					$correct_ori++;
					if (defined $ARGV[2]) {
						if (defined $LEFT_QUALS{$P[0]}) {
							$GOOD_ORI{$ORIENT_QUALS{$P[0]}[$i]}++;
							print OUTQO "1\t$ORIENT_QUALS{$P[0]}[$i]\n";
						}
					}
				} else {
					if (defined $ARGV[2]) {
						if (defined $LEFT_QUALS{$P[0]}) {
							$BAD_ORI{$ORIENT_QUALS{$P[0]}[$i]}++;
							print OUTQO "0\t$ORIENT_QUALS{$P[0]}[$i]\n";
						}
					}
				}
				$ori_total++;
				$misassembly_incl++;
				$CHR_HITS{$chr}++;
				$BASE_HITS{$chr} += ($aln2-$aln1);
				$bases_total += ($aln2-$aln1);
				$scaff_bases_total += ($aln2-$aln1);
				$total_scaff_bases += ($aln2-$aln1);
			}
			$pos = $O[$i];
		}
		$ori1 = $O1[$i]; $ori2 = $O2[$i];
		$total++;
	}
	if ($scaff_bases_total>0) {
		$fraction_top = -1; $check_second = 0; $miss_check = 0;
		foreach $base_hit (sort {$BASE_HITS{$b}<=>$BASE_HITS{$a}} %BASE_HITS) {
			if ($fraction_top<-0.5) {
				$max_chr_count = $CHR_HITS{$base_hit};
				$fraction_top = $BASE_HITS{$base_hit}/$scaff_bases_total;
				if ($fraction_top<0.9) {$miss_check=1}
			} elsif ($check_second < 0.5) {
				$check_second = 1;
				if ($CHR_HITS{$base_hit} == 1 && $max_chr_count>=2) {
					$miss_check=0;
				}
			}
		}
		if ($miss_check > 0.5) {$total_miss++};
		$total_for_miss++;
	}
	if ($non_perfect_scaffold_check>0.5) {
		$non_perfect_scaffolds++;
	} elsif ($total_scaff_bases > $longest_perfect_size) {
		$longest_perfect_size = $total_scaff_bases
	}
	if ($perfect_stretch>$longest_perfect_stretch) {
		$longest_perfect_stretch=$perfect_stretch;
	}
} close IN; close OUTQ, close OUTQO;

$place_frac = sprintf("%.2f", ($new_correct/$new_total)*100);
$orient_frac = sprintf("%.2f", ($correct_ori/$ori_total)*100);
$miss_percent = sprintf("%.2f", ($total_miss/$total_for_miss)*100);
$proper_bases_percent = sprintf("%.2f", ($bases_prop_tally/$bases_total)*100);
$perfect_scaffold_percent = sprintf("%.2f", (($total_scaffolds-$non_perfect_scaffolds)/$total_scaffolds)*100);

$stats = $out; $stats =~ s/txt/stats/;
open OUT, ">$stats";

$ts = localtime(time);
print OUT "
RUN AT $ts

TOTAL CONTIGS SCAFFOLDED = $global_orient_total
TOTAL SCAFFOLDS          = $total_scaffolds
TOTAL JOINS              = $total
TOTAL TRUE MAPPED        = $mapped_contig_count

PLACEMENT ACCURACY       = $place_frac
ORIENTATION ACCURACY     = $orient_frac
PERFECT SCAFFOLDS        = $perfect_scaffold_percent
FUSED SCAFFOLDS          = $miss_percent
BASES PROP. PLACED       = $proper_bases_percent
LONGEST PERFECT SCAFF    = $longest_perfect_size
LONGEST PERFECT JOIN     = $longest_perfect_stretch

";

close OUT;


















