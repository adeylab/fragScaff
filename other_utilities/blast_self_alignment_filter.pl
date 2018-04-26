#!/usr/bin/perl

$die = "

ARGV0 = self-blast
ARGV1 = percent identity min (to be a repeat)

STDOUT = bed of repeats (need smerging)

After filtering run:

sortBed -i MY_OUT.bed \> MY_OUT.srt.bed
mergeBed -i MY_OUT.srt.bed -nms \> MY_REPEATS.bed

";

if (!defined $ARGV[1]) {die $die}

#super_0 super_0 91.13   1613    123     17      799569  801172  311520  313121  0.0     2169
#super_0 super_0 91.13   1613    123     17      311520  313121  799569  801172  0.0     2169
#super_0 super_0 87.64   1813    141     42      367052  368806  191053  192840  0.0     2030

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[2] >= $ARGV[1]) {
		if ($P[0] eq $P[1]) {
			if (abs($P[6]-$P[8])>100&&abs($P[7]-$P[9])>100) {
				print "$P[0]\t$P[6]\t$P[7]\t$P[1]:$P[8]-$P[9]\n";
#				print "$P[1]\t$P[8]\t$P[9]\t$P[0]:$P[6]-$P[7]\n";
			}
		} else {
			print "$P[0]\t$P[6]\t$P[7]\t$P[1]:$P[8]-$P[9]\n";
#			print "$P[1]\t$P[8]\t$P[9]\t$P[0]:$P[6]-$P[7]\n";
		}
	}
} close IN;