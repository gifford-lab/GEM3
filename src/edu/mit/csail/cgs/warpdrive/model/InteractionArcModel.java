package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.projects.readdb.PairedHitLeftComparator;

public class InteractionArcModel extends WarpModel implements RegionModel,
Runnable {

	private Client client;
	private Set<ChipSeqAlignment> alignments;
	private Set<String> ids;
	private Region region;
	private boolean newinput;
	private List<PairedHit> results, otherchrom;
	private Comparator<PairedHit> comparator;
	private InteractionArcProperties props;

	public InteractionArcModel(Collection<ChipSeqAlignment> alignments) throws IOException, ClientException {
		client = new Client();
		comparator = new PairedHitLeftComparator();
		this.alignments = new HashSet<ChipSeqAlignment>();
		this.alignments.addAll(alignments);
		ids = new HashSet<String>();
		for (ChipSeqAlignment a : alignments) {
			ids.add(Integer.toString(a.getDBID()));
		}
		results = null;
		otherchrom = null;
		props = new InteractionArcProperties();
	}

	public InteractionArcProperties getProperties() {return props;}

	public void clearValues() {
		results = null;
		otherchrom = null;
	}

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
					otherchrom = new ArrayList<PairedHit>();
					for (String alignid : ids) {
						List<PairedHit> r = client.getPairedHits(alignid,
								region.getGenome().getChromID(region.getChrom()),
								true,
								region.getStart(),
								region.getEnd(),
								null,
								null);
						for (PairedHit h : r) {
							if (h.leftChrom == h.rightChrom) { 
								if (h.rightPos >= region.getStart() &&
										h.rightPos <= region.getEnd()) {
									results.add(h);
								}
							} else {
								otherchrom.add(h);
							}
						}
					}
					for (PairedHit h : results) {
						if (h.leftPos > h.rightPos) {
							h.flipSides();
						}
					}
					Collections.sort(results, comparator);
					Collections.sort(otherchrom, comparator);
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

	public Region getRegion() {return region;}

	public List<PairedHit> getResults () {return results;}
	public List<PairedHit> getOtherChromResults() {return otherchrom;}

}
