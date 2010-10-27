package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.projects.readdb.PairedHitLeftComparator;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class PairedEndModel extends WarpModel implements RegionModel, Runnable {

    private Client client;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    private Region region;
    private boolean newinput;
    private List<PairedHit> results;
    private Comparator<PairedHit> comparator;
    private PairedEndProperties props;

    public PairedEndModel (Collection<ChipSeqAlignment> alignments) throws IOException, ClientException{
        client = new Client();
        comparator = new PairedHitLeftComparator();
        this.alignments = new HashSet<ChipSeqAlignment>();
        this.alignments.addAll(alignments);
        ids = new HashSet<String>();
        for (ChipSeqAlignment a : alignments) {
            ids.add(Integer.toString(a.getDBID()));
        }
        results = null;
        props = new PairedEndProperties();
    }
    public PairedEndProperties getProperties() {return props;}

    public void clearValues() {
        results = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public List<PairedHit> getResults () {return results;}
    public boolean isReady() {return !newinput;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) { }
            if (newinput) {
                try {
                    results = new ArrayList<PairedHit>();
                    double mindist = getProperties().MinimumDistance;
                    boolean showself = getProperties().ShowSelfLigation;
                    if (mindist < 1) {
                        mindist = mindist * region.getWidth();
                    }
                    for (String alignid : ids) {
                        List<PairedHit> r = client.getPairedHits(alignid,
                                                                 region.getGenome().getChromID(region.getChrom()),
                                                                 true,
                                                                 region.getStart(),
                                                                 region.getEnd(),
                                                                 null,
                                                                 null);
                        for (PairedHit h : r) {
                            if (h.leftChrom == h.rightChrom && 
                                h.rightPos >= region.getStart() &&
                                h.rightPos <= region.getEnd() &&
                                Math.abs(h.leftPos - h.rightPos) > mindist &&
                                (showself || !isSelfLigation(h))) {
                                results.add(h);
                            }
                        }
                    }
                    Collections.sort(results, comparator);
                    if (getProperties().PrintData) {
                        for (PairedHit h : results) {
                            System.out.println(h.toString());
                        }

                    }

                    if (getProperties().DeDuplicateByPosition && results.size() > 0) {
                        ArrayList<PairedHit> deduped = new ArrayList<PairedHit>();
                        deduped.add(results.get(0));
                        for (int i = 1; i < results.size(); i++) {
                            PairedHit a = results.get(i);
                            PairedHit b = deduped.get(deduped.size() - 1);
                            if (a.leftPos != b.leftPos ||
                                a.rightPos != b.rightPos ||
                                a.leftStrand != b.leftStrand ||
                                a.rightStrand != b.rightStrand) {
                                deduped.add(a);
                            }
                        }
                        results = deduped;
                    }
                    if (getProperties().LeftAlwaysLesser) {
                        for (PairedHit h : results) {
                            if (h.leftPos > h.rightPos) {
                                h.flipSides();
                            }
                        }
                        Collections.sort(results, comparator);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    // assign empty output.  This is useful because Client
                    // throws an exception for non-existant chromosomes, such
                    // as those for which there were no alignment results
                    results = new ArrayList<PairedHit>();
                }
                newinput = false;
                notifyListeners();
            }
        }
        client.close();
    }                     

    public boolean isSelfLigation(PairedHit p) {
    	if (getProperties().RightFlipped) {
    		return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= getProperties().SelfLigationCutoff) && (p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
    		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    	} else {
    		return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= getProperties().SelfLigationCutoff) && !(p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
    		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    	}
    }

}
