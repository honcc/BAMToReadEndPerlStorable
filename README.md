BAMToReadEndPerlStorable
========================

This is a perl script to read a BAM file (using samtools) and generate the pileup in format of perl storables. Users can choose to pileup only the read ends, mid-point or full read. It uses multi-thread to process the chromosome independently to increase the speed.
