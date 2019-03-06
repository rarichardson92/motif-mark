# motif-mark

Reads in an **RNA/DNA FASTA file or UCSC coordinates** along with desired RNA/DNA motifs and outputs one of the following: visually projects where motifs are in the sequence, provides a tabular list of indexes for motif and exon locations in the sequence on the command line, tsv with output similar to commandline printout. 

**First base index position = 1**, reverse positions based on reverse compliment with the same orientation as the initial sequence (From the left of the reverse compliment). 

Motif matches are insensitive to RNA/DNA file types (i.e., an RNA motif will match to a DNA FASTA and vice-versa). Defaults to FASTA input and an image output. 
**Exons are read as capitalized sequence from FASTA. Motifs should be seperated by newline characters.** 

Use of UCSC options requires download of UCSC 2bit file and twoBitToFa function (if not already present), thus requires internet access. Any pulled UCSC sequence can be found in UCSC.fa file after program completion. 

Imports argparse, re, and cairo, please **install packages appropriately before use**.



For **custom colors and image parameters**, please output as a tsv and access the jupyter notebook.

.svg and .tsv files are example outputs using motifs in Fig_1_motifs.txt and Figure_1.fasta.
