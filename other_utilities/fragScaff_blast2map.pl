#!/usr/bin/perl

#example:
#~/tools/ncbi-blast-2.2.27+/bin/makeblastdb -in P_HALLI_MAIN_1M.ASSEMBLY.fa  -dbtype nucl
#~/tools/ncbi-blast-2.2.27+/bin/blastn -outfmt 6 -query P_HALLI_NO_FES_BES.ASSEMBLY.fa -db P_HALLI_MAIN_1M.ASSEMBLY.fa -word_size 80 -num_threads 20 -evalue 100 2> P_HALLI_NO_FES_BES.ASSEMBLY.blast.err > P_HALLI_NO_FES_BES.ASSEMBLY.blast.out &

$die = "

ARGV0 = dict file for query
ARGV1 = blast output

STDOUT = outmap

";

if (!defined $ARGV[1]) {die $die}

#@SQ     SN:contig_0     LN:7453 UR:file:/net/shendure/vol10/nobackup/acadey/FRAGSCAFF/HUMAN/NA12878_Assembly_Contigs.fa M5:0b20bd520d8968ced96672ce36a99eaf
#@SQ     SN:contig_1     LN:116352       UR:file:/net/shendure/vol10/nobackup/acadey/FRAGSCAFF/HUMAN/NA12878_Assembly_Contigs.fa M5:e56c7babb804ca1b815f0f9956adf266
#@SQ     SN:contig_2     LN:5263 UR:file:/net/shendure/vol10/nobackup/acadey/FRAGSCAFF/HUMAN/NA12878_Assembly_Contigs.fa M5:684c9db66978f882b926a47da1c7b702

open DICT, "$ARGV[0]";
while ($l = <DICT>) {
	chomp $l;
	if ($l =~ /^\@SQ/) {
		@P = split(/\t/, $l);
		$name = $P[1]; $name =~ s/^SN://;
		$length = $P[2]; $length =~ s/^LN://;
		$LENGTH{$name} = $length;
		$WINFRAC{$name} = 0;
	}
} close DICT;

#contig_403      13      89.93   1361    118     12      155998  157342  82436843        82435486        0.0     1748
#contig_403      13      90.84   1310    99      9       156053  157341  48781089        48779780        0.0     1744

open BLAST, "$ARGV[1]";
while ($l = <BLAST>) {
	chomp $l;
	@P = split(/\t/, $l);
	$check = $P[3]/$LENGTH{$P[0]};
	if (!defined $WINFRAC{$P[0]}) {
		$WINFRAC{$P[0]} = $check;
		$WINNER{$P[0]} = "$P[1]\t$P[8]\t$P[9]"
	} else {
		if ($check > $WINFRAC{$P[0]}) {
			$WINFRAC{$P[0]} = $check;
			$WINNER{$P[0]} = "$P[1]\t$P[8]\t$P[9]";
		}
	}
} close BLAST;

foreach $contig (keys %LENGTH) {
	if ($WINFRAC{$contig}>0.8) {
		print "$contig\t$WINNER{$contig}\t$WINFRAC{$contig}\n";
	} else {
		print "$contig\t-1\t-1\t-1\t$WINFRAC{$contig}\n";
	}
}