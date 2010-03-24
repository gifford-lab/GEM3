package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ChipSeqExpander implements Expander<Region, ChipSeqHit>, Closeable {

    private ChipSeqLoader loader;
    private Genome lastGenome;
    private LinkedList<ChipSeqAlignment> alignments;
    private ChipSeqLocator locator;
    private boolean closeLoader;


    public ChipSeqExpander(ChipSeqLocator loc) throws SQLException, IOException {
        loader = new ChipSeqLoader();
        closeLoader = true;
        locator = loc;
        alignments = null;
        lastGenome = null;
    }
    private void getAligns(Genome genome) throws SQLException {
        if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<ChipSeqAlignment>();
        lastGenome = genome;
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    public ChipSeqExpander(ChipSeqLoader l, ChipSeqAlignment a, boolean closeLoader) {
        loader = l;
        alignments = new LinkedList<ChipSeqAlignment>();
        alignments.add(a);
        this.closeLoader = closeLoader;
    }


    public Iterator<ChipSeqHit> execute(Region a) {
        try {
            getAligns(a.getGenome());
            Collection<ChipSeqHit> hits = loader.loadByRegion(alignments, a);
            return hits.iterator();
        }
        catch (Exception e) {
            e.printStackTrace();
            return new LinkedList<ChipSeqHit>().iterator();
        }
    }


    public int getHitCount(Region a) {
        int hits = 0;
        try {
            getAligns(a.getGenome());
            hits = loader.countByRegion(alignments, a);
            return hits;
        }
        catch (Exception e) {
            e.printStackTrace();
            return 0;
        }
    }


    public Collection<Genome> alignedGenomes() {
        LinkedList<Genome> genomes = new LinkedList<Genome>();
        if (alignments != null) {
            for (ChipSeqAlignment align : alignments) {
                genomes.add(align.getGenome());
            }
        }
        return genomes;
    }


    public void close() {
        if (closeLoader) {
            loader.close();
        }
        loader = null;
        if (alignments != null) {
            alignments.clear();
        }
    }


    public boolean isClosed() {
        return loader == null;
    }
}
