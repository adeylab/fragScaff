
============================= fragScaff ==============================

   Version 140324.1 (180112 10X support update - beta -G X option)

                   Andrew Adey (adey\@ohsu.edu)

UPDATE: Date: Jan. 12, 2018. Updated to allow 10X genomics bam files
        by detecting the CR:Z field and building the RG info in
        pre-processing. Note: This has not been fully tested.
        For 10X data, the @RG lines in the header are NOT required.

Attempts to scaffold contigs using fragment pool data. Performs best
using CPT-Seq (Contiguity Preserving Transposase) libraries, but will
work using fosmid (decent) or LFR (Long Fragment Read) (poor) fragment
data, though these methods are far less tested in fragScaff. It is
required that the input bam has the proper @SQ lines, as well as @RG
lines and are aligned to the reference assembly provided as the -F
option. The bam (-B) must be reads aligned to the input assembled fasta
(-F).

fragScaff.pl -B <BAM/bamParse> -F <Input fasta> <Other options>

NOTE: It is recommended to run as:
'fragScaff.pl -B <BAM> -b 1 -E <node size>'

to first generate the bamParse file. This file can then be used for all
subsequent stages to eliminate the need to re-parse the bam. Similarly
it is recommended to run the node-node calculations once to generate a
'.links.txt' file which can then be used as option '-K' for all
subsequent runs thus eliminating the need to rerun the time-consuming
node-node calculations.

RECOMMENDED RUNS:

(1) 'fragScaff.pl -B <myBAM> -b 1 -E <node_size>
     (optional: -J <repeats.bed> -N <Nbases.bed> -o <max_node_size>)'


(2) 'fragScaff.pl -B <myBAM.E#.bamParse> -A -O <myOUT>
     -t <myTHREAD#/Q>'

(3) 'fragScaff.pl -B <myBAM.E#.bamParse>
     -K <myOUT.fragScaff.links.txt> -F <myIN.fa> -O <myOUT2>
     <other options>'

This will allow you to run all pre-processing in steps (1) and (2) and
then run the last stage (3) multiple times to compare outputs and find
the optimum setting for your input. Though defaults should give an
optimum output. Also note that step (3) can be run multiple times
without the -F option and with -I so that it will only perform edge
filtering and graph manipulations which require the most optimization.
Once the final otimized parameters are decided based on the stats
printed in the log file, it can then be run with -F to specify the
input assembly fasta and without -I to produce the final output
assembly.

-H   Display the following detailed options description.

----------------- INPUT OPTIONS --------------------------------------

-B   [FILE]   bam file (must have SQ (w/LN) & RG tags) (REQ)
                (can be replaced with bamParse file)
-O   [STR]    output prefix (def = -B)
-F   [FILE]   fasta file from assembler output
-P   [T/F/L]  platform: T = CPT-Seq / 10X (def), F = Fosmid Pool,
                L = Long Fragment Read
-G   [R/N/H]  read group identifier: N = read name (group:read_number)
                (def), R = RG:Z:group bam tag, H = after hash
                (name#group), X = CR:Z group (10X bam)
-K   [FILE]   previously generated link file.
-N   [FILE]   Nbase bed file for input scaffolds
-J   [FILE]   Repeatmasker bed file to exclude reads in windows

----------------- PRUNING OPTIONS ------------------------------------

-q   [INT]    min mapping qual (def = 10)
-m   [INT]    min contig size to include (def = 1)
-d   [FLOAT]  min fraction group hit cutoff (def = 0.05)
-l   [INT]    max number of links to use (def = 5)
-D   [FLOAT]  max fraction group hit cutoff (def = 0.95)
-a   [INT]    max number of links to allow (def = 20)

-U   [INT]    min num group hits per node (def = 0)
                (will remove nodes with more)

-p   [INT/A]  -log10 score use minimum (def = A) (set to A for auto
                determine)
-j   [FLOAT]  if -p A: mean links per p-bin (def = 1.25)
-c   [INT/D]  size of brach to check for split (def = D)
 (set to D to
                disable, experimental)
-r   [INT]    -log10 score report minimum (def = 1)
-M   [INT]    max score (def = 200)
-u   [FLOAT]  combined reciprocation score multiplier (def = 2)
-g   [INT/X]  max size to make join (ie. one scaff must be less than
                this, def = X (null))

----------------- PLATFORM-DEPENDENT OPTIONS -------------------------

-E   [INT]    contig end node size (T/L=5000, F=1000) (Important for
                bamParse creation)
-C   [INT]    read count threshold (T=1, F/L=20)
-o   [INT]    max E (if dynamic) (T=10000, F/L=5000)

----------------- OUTPUT OPTIONS -------------------------------------

-b   [0/1]    generate bamParse file to remove need to parse bam each
                run (...bam.bamParse): 0 = make bamParse & continue
                1 = make bamParse then exit (can be run with just -b 1
                -B <bam>)
-n   [n/N]    case of N to add (def = N)
-L   [INT]    wrap length in output fasta (def = 100)
-v            print/keep intermediate files (def = no)
-f            print assembly qualities in fasta file (def = separate)
-R            print out link fraction overlaps (def = no)
-z            gzip bamParse (gzip must be command line callable)
-x   [INT]    min N spacer size (rec = max MP size, def = 3000)
-X   [INT]    max N spacer size (def = 8000)
-V   [INT/A]  print [INT] graph files (def = 0, A for all)
-w            exit after printing cytoscape graph (def = no, req -V)

----------------- PERFORMANCE OPTIONS --------------------------------

-t   [INT/Q]  threads (def = 1); runs with 1 during bam parsing then
                expands to multiple. The mutithread portion may
                require a large amount of memory. Prepare for 2G each.
                assign to Q to qsub jobs. (best)
-T   [INT]    if -t Q: number of qsub jobs at a time (def = 100)
                set to 'N' for no limit
-e   [STR]    if -t Q: memory per job (def = 2G)
-Q   [STR]    for qsub - do not manually toggle
-S   [INT]    if -t Q/>1: number of nodes per job (def = 100)
-A            exit after node all-by-all calculations
-I            exit after graph manipulation

----------------- DEPENDENCY OPTIONS ---------------------------------

-k   [STR]    fragScaff call (def = fragScaff.pl) (only necessary if
                -t >1 or Q)
-s   [STR]    samtools call (def = samtools)

