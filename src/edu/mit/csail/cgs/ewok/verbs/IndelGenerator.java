package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * IndelGenerator uses multizalign-style alignments (stored as blocks
 * from some alignment algorithm between a source genome and a target genome) and
 * dynamically computes the set of indels in a region.
 *
 * Alex should add a better description of how the alignment algorithm works
 */
public class IndelGenerator<X extends Region> implements Expander<X, Indel> {

    private Genome compareTo, compareFrom;
    private MultiZAlignGenerator<X> multiz;
    private String alignPrefix;
    /* indels bigger than this don't count.  This is to help
       filter out places where the a section of one genome aligns to
       several sections of the other.  By not allowing huge jumps,
       we try to retain more continuity in the alignement */
    private int maxPlausibleChange = 10000;

    private int splitPenalty = -1, chromPenalty = -1;
    
    public IndelGenerator(Genome compareTo) {
        this.compareTo = compareTo;
        multiz = null;
        alignPrefix = null;
    }
    public String getAlignPrefix() {return alignPrefix;}
    public void setAlignPrefix(String p) {
        alignPrefix = p;
        if (multiz != null) {
            multiz.setAlignPrefix(p);
        }
    }
    
    private void createMultiZ (Genome compareFrom) {
        multiz = new MultiZAlignGenerator<X>(compareFrom, compareTo);
        if (alignPrefix != null) {
            multiz.setAlignPrefix(alignPrefix);
        }
        this.compareFrom = compareFrom;
    }
    public void setSplitPenalty(int s) {splitPenalty = s;}
    public void setChromPenalty(int c) {chromPenalty = c;}

    public static List<Indel> getIndels(Genome genome, String chrom, int maxPlausibleChange, List<? extends MultiZAlignRegion> alignedto) {

        ArrayList<Indel> out = new ArrayList<Indel>();

        if (alignedto.size() == 1) {
            // can't tell anything, so return nothing.
            return out;
        }

        Collections.sort(alignedto);
        for (MultiZAlignRegion b : alignedto) {
            System.err.println(b);
        }

        /* You can't just run through alignedto in order; you'll mess
           up inversions.  So the BlockChunks are chunks of alignedto
           and the orientation in which to process them.  
        */
        ArrayList<BlockChunk> chunks = new ArrayList<BlockChunk>();
        int index = 0;
        int laststart = -1;
        char laststrand = ' ';
        while (index < alignedto.size()) {
            char strand = alignedto.get(index).getStrand();
            if (strand != laststrand) {
                if (laststart >= 0) {
                    chunks.add(new BlockChunk(laststart,
                                              index - 1,
                                              laststrand));
                }
                laststart = index;
            }
            laststrand = strand;
            index++;
        }
        chunks.add(new BlockChunk(laststart,
                                  index - 1,
                                  laststrand));

        MultiZAlignRegion prevregion = null;
        int prevend = 0, prevotherend = 0;

        for (BlockChunk chunk : chunks) {
            System.err.println(String.format("CHUNK %d %d %c\n", chunk.start, chunk.end, chunk.direction));
            if (chunk.direction == '-') {
                int s = Math.min(alignedto.get(chunk.start).getStart(),
                                 alignedto.get(chunk.end).getEnd());
                int e = Math.max(alignedto.get(chunk.start).getStart(),
                                 alignedto.get(chunk.end).getEnd());
                out.add(new Indel(genome,
                                  chrom,
                                  s,
                                  e,
                                  e - s,
                                  Indel.INVERSION));
            }

            for (int i = chunk.start; i != chunk.end; i += (chunk.direction == '+' ? 1 : -1)) {
                if (prevregion == null) {
                    prevregion = alignedto.get(i);
                    prevend = prevregion.getEnd();
                    prevotherend = prevregion.getOtherStart() < prevregion.getOtherEnd() ?
                        prevregion.getOtherEnd() : prevregion.getOtherStart();
                } else {
                    MultiZAlignRegion thisregion =alignedto.get(i);
                    int thisstart = thisregion.getStart();
                    int thisotherstart = thisregion.getOtherStart() < thisregion.getOtherEnd() ?
                        thisregion.getOtherStart() : thisregion.getOtherEnd();
                    int gapsize  = thisstart - prevend;
                    int othergapsize = thisotherstart - prevotherend;
                    System.err.println("region " + thisregion + "   " + gapsize +", " + othergapsize);
                    System.err.println("  prevend=" + prevend +" thisstart=" + thisstart + 
                                       " prevotherend="+ prevotherend + " thisotherstart=" + thisotherstart);
                    if ((Math.abs(gapsize - othergapsize) < maxPlausibleChange) &&
                        (thisstart > prevend) &&
                        (gapsize != othergapsize)) {
                        if (gapsize > othergapsize) {
                            out.add(new Indel(genome,
                                              chrom,
                                              prevend,
                                              thisstart,
                                              gapsize - othergapsize,
                                              Indel.DELETION));
                        } else {
                            out.add(new Indel(genome,
                                              chrom,
                                              prevend,
                                              thisstart,
                                              othergapsize - gapsize,
                                              Indel.INSERTION));
                        }                        
                    }
                    if ((gapsize < 0 || othergapsize < 0) &&
                        (Math.signum(gapsize) * Math.signum(othergapsize) * (chunk.direction == '+' ? 1 : -1) < 0)) {
                        out.add(new Indel(genome,
                                          chrom,
                                          Math.min(prevend,thisstart),
                                          Math.max(prevend,thisstart),
                                          Math.abs(Math.min(gapsize,othergapsize)),
                                          Indel.REARRANGEMENT));
                    }

                    
                    prevregion = thisregion;
                    prevend = prevregion.getEnd();
                    prevotherend = prevregion.getOtherStart() < prevregion.getOtherEnd() ?
                        prevregion.getOtherEnd() : prevregion.getOtherStart();
                    
                }                               
            }
        }
        return out;
    }

    /**
     * Aligns this region to the genome compareTo.  Finds the best
     * matching chromosome as that which has the highest sum of
     * scores.  Returns indels of region relative to the aligned
     * regions on that chromosome.
     */
    public Iterator<Indel> execute(X region) {
        if (!region.getGenome().equals(compareFrom)) {
            createMultiZ(region.getGenome());
        }
        Iterator<MultiZAlignRegion> iter = multiz.execute(region);
        ArrayList<MultiZAlignRegion> list = new ArrayList<MultiZAlignRegion>();
        while (iter.hasNext()) {
            list.add(iter.next());
        }
        AlignmentStitcher<MultiZAlignRegion> stitcher = new AlignmentStitcher<MultiZAlignRegion>(list);
        if (splitPenalty > 0) {
            stitcher.splitPenalty(splitPenalty);
        }
        if (chromPenalty > 0) {
            stitcher.chromSwitchPenalty(chromPenalty);
        }

        stitcher.stitch();
        List<MultiZAlignRegion> alignedto = stitcher.getBestAlignment(stitcher.getBestAlignedChrom());

        if (alignedto.size() == 0) {
            ArrayList<Indel> out = new ArrayList<Indel>();
            out.add(new Indel(region.getGenome(),
                              region.getChrom(),
                              region.getStart(),
                              region.getEnd(),
                              region.getWidth(),
                              Indel.DELETION));
            return out.iterator();
        }

        return getIndels(region.getGenome(), region.getChrom(), maxPlausibleChange, alignedto).iterator();
    }

    public static int center(MultiZAlignRegion r) {
        return (int)((r.getOtherEnd() + r.getOtherStart()) / 2);
    }
    
    public static void main(String args[]) throws NotFoundException {
        String prefix = "chain";
        ArrayList<Region> regions = new ArrayList<Region>();
        ArrayList<Genome> genomes = new ArrayList<Genome>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--prefix")) {
                prefix = args[++i];
            }
            if (args[i].equals("--species")) {
                String[] pieces = args[++i].split(";");
                if (pieces.length == 2) {
                    Organism org = new Organism(pieces[0]);
                    Genome genome = org.getGenome(pieces[1]);
                    genomes.add(genome);
                }
            }
        }
        if (genomes.size() < 2) {
            System.err.println("Must supply two --species arguments");
            System.exit(1);
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--region")) {
                regions.add(Region.fromString(genomes.get(0),args[++i]));
            }
        }
        IndelGenerator<Region> gen = new IndelGenerator<Region>(genomes.get(1));
        gen.setAlignPrefix(prefix);
        for (Region r : regions) {
            Iterator<Indel> iter = gen.execute(r);
            while (iter.hasNext()) {
                Indel indel = iter.next();
                System.err.println(indel.toString());
            }

        }

    }
}

class BlockChunk {
    public int start, end;
    public char direction;
    public BlockChunk( int s, int e, char d) {
        direction = d;
        if (direction == '+') {
            start = Math.min(s,e);
            end = Math.max(s,e);
        } else {
            start = Math.max(s,e);
            end = Math.min(s,e);
        }
    }
}
    