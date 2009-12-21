package edu.mit.csail.cgs.datasets.alignments;

import edu.mit.csail.cgs.datasets.species.Genome;

/* From the MAQ Manpage:
For reads aligned before the Smith-Waterman alignment, each
line consists of 
1) read name
2) chromosome
3) position
4) strand
5) insert size from the outer coorniates of a pair
6) paired flag
7) mapping quality
8) single-end mapping quality
9) alternative mapping quality,
10) number of mismatches of the best hit
11) sum of qualities of mismatched bases of the best hit
12) number of 0-mismatch hits of the first 24bp
13) number of 1-mismatch hits of the first 24bp on the reference
14) length of the read
15) read sequence
16) its quality

Alternative mapping quality always equals to mapping
quality if the reads are not paired. If reads are paired, it
equals to the smaller mapping quality of the two ends. This
alternative mapping quality is actually the mapping quality of an
abnormal pair.

The fifth column, paired flag, is a bitwise flag. Its lower 4
bits give the orientation: 1 stands for FF, 2 for FR, 4 for RF,
and 8 for RR, where FR means that the read with smaller
coordinate is on the forward strand, and its mate is on the
reverse strand. Only FR is allowed for a correct pair. The higher
bits of this flag give further information. If the pair meets the
paired end requirement, 16 will be set. If the two reads are
mapped to different chromosomes, 32 will be set. If one of the
two reads cannot be mapped at all, 64 will be set. The flag for a
correct pair always equals to 18.

For reads aligned by the Smith-Waterman alignment afterwards, the
flag is always 130. A line consists of read name, chromosome,
position, strand, insert size, flag (always 130), position of the
indel on the read (0 if no indel), length of the indels (positive
for insertions and negative for deletions), mapping quality of
its mate, number of perfect hits, number of 1-mismatch hits on
the reference, length of the read, read sequence and its quality.

*/

/* Since MAQ is not a general purpose alignment program, some of the fields in 
   BLATAlignmentHitRegion are a little wierd.  For example, MAQ just gives a 
   single position so we fudge the actual position from the read length.
*/

public class MAQHitRegion extends BLATAlignmentHitRegion {

    private int mappingQ, singleEndedMappingQ, alternateMappingQ;

    public MAQHitRegion(Genome genome,
                        String readName,
                        String chrom,
                        int position,
                        char strand,
                        int insertSize,
                        boolean paired,
                        int mappingQuality,
                        int singleEndedMappingQuality,
                        int alternateMappingQuality,
                        int numMismatches,
                        int readLength) {
        super(genome,
              chrom,
              position,
              position+readLength,
              readName,
              strand,
              (readLength - numMismatches) / ((double)readLength),
              readLength,
              numMismatches,
              0,
              1,
              readLength);
        mappingQ = mappingQuality;
        singleEndedMappingQ = singleEndedMappingQuality;
        alternateMappingQ = alternateMappingQuality;
    }

    public int getMappingQuality() {return mappingQ;}
    public int getSingleEndedMappingQuality() {return singleEndedMappingQ;}
    public int getAlternateMappingQuality() {return alternateMappingQ;}

}