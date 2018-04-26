#!/usr/bin/perl

$die = "

ARGV0 = X.ordered.txt
ARGV1 = X.quals
ARGV2 = X.bam

STDOUT = scaffolds using original IDs

";

if (!defined $ARGV[2]) {die $die}

$id = 0;
open IN, "samtools view -H $ARGV[2] |";
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^\@SQ/) {
		@P = split(/\t/, $l);
		$name = $P[1]; $name =~ s/^SN://;
		$ID_2_NAME{$id} = $name;
		$id++;
	}
} close IN;

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$P[0] =~ s/^>//;
	($null,$ID) = split(/=/, $P[1]);
	if ($P[0] =~ /original/) {
		$ORIGINAL{$P[0]} = $ID;
	} else {
		$OUT_NAME{$ID} = $P[0];
	}
} close IN;

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	@S = split(/,/, $P[1]);
	($num,$ori) = split(/\./, $S[0]);
#	$out_scaff = "scaffold_$num";
	$out_scaff = $ID_2_NAME{$num};
	$out_ori = "$ori";
	for ($i = 1; $i < @S; $i++) {
		($num,$ori) = split(/\./, $S[$i]);
#		$out_scaff .= ",scaffold_$num";
		$out_scaff .= ",".$ID_2_NAME{$num};
		$out_ori .= ",$ori";
	}
	print "$OUT_NAME{$P[0]}\t$out_scaff\t$out_ori\n";
} close IN;

foreach $original (keys %ORIGINAL) {
	print "$original\t$ORIGINAL{$original}\tf\n";
}