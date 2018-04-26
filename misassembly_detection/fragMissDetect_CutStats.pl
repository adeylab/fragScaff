#!/usr/bin/perl

# INSTEAD FOR CHECKING JUST USE BEDTOOLS W/ FASTA

$die = "

ARGV0 = map.txt from bedtools check

STDOUT = stats

";

open IN, "$ARGV[0]";
while ($site1 = <IN>) {
	chomp $site1;
	@P1 = split(/\t/, $site1);
	$site2 = <IN>;
	chomp $site2;
	@P2 = split(/\t/, $site2);
	($name1,$pos1) = split(/:/, $P1[0]);
	($name2,$pos2) = split(/:/, $P2[0]);
	while ($name1 ne $name2) {
		($name1,$pos1) = ($name2,$pos2);
		@P1 = @P2;
		$site1 = $site2;
		$site2 = <IN>;
		chomp $site2;
		@P2 = split(/\t/, $site2);
		($name2,$pos2) = split(/:/, $P2[0])
	}
	$total_joins++;
	if ($P1[1] =~ /-1/ || $P2[1] =~ /-1/) {
		$input_misassemblies++;
	} elsif ($P1[1] eq $P2[1] && abs($P1[2]-$P2[2]) < 100000) {
		$proper_joins++;
	} else {
		$fragScaff_misassemblies++;
	}
} close IN;

$prop_pct = sprintf("%.2f", ($proper_joins/$total_joins)*100);
$input_pct = sprintf("%.2f", ($input_misassemblies/$total_joins)*100);
$fragScaff_pct = sprintf("%.2f", ($fragScaff_misassemblies/$total_joins)*100);

print "

TOTAL CUTS           =	$total_joins
CUTS AT:
   CORRECT REGIONS   =	$proper_joins	$prop_pct
   INPUT MISS.       =	$input_misassemblies	$input_pct
   FRAGSCAFF MISS.   =	$fragScaff_misassemblies	$fragScaff_pct

";
