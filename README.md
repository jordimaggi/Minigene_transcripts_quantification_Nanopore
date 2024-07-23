# Splice junctions quantifications Nanopore

Repository for scripts used to clean NanoSplicer output data. The script uses a `.gff3` file containing coordinates of known exons included in the minigene constructs to identify and annotate transcripts based on the splice junctions detected by the `JWR_checker.py` script from [NanoSplicer](https://github.com/shimlab/NanoSplicer).

The script's input include a gff3 file containing coordinates of known exons included in the minigene constructs, and a .csv file created by the `JWR_checker.py` script from [NanoSplicer](https://github.com/shimlab/NanoSplicer) with the following command:



The script's outputs include tables and gff3 for reference and variant minigenes that lists the different transcripts identified and the number of reads representing them.
