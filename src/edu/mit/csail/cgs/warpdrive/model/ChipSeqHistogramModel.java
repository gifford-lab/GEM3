package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * Data model for chipseq histogram.  Separate methods for retrieving
 * plus and minus strand results
 */
public class ChipSeqHistogramModel extends WarpModel implements RegionModel, Runnable {
    
    private Client client;
    private TreeMap<Integer,Float> resultsPlus, resultsMinus;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    private ChipSeqHistogramProperties props;

    private Region region;
    private boolean newinput;

    public ChipSeqHistogramModel (ChipSeqAlignment a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.add(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        ids = new HashSet<String>();
        ids.add(Integer.toString(a.getDBID()));
    }
    public ChipSeqHistogramModel (Collection<ChipSeqAlignment> a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.addAll(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        ids = new HashSet<String>();
        for (ChipSeqAlignment align : alignments) {
            ids.add(Integer.toString(align.getDBID()));
        }
    }    
    public ChipSeqHistogramProperties getProperties() {return props;}
    
    
    public void clearValues() {
        resultsPlus = null;
        resultsMinus = null;
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
    public boolean isReady() {return !newinput;}
    public Map<Integer,Float> getPlus() {return resultsPlus;}
    public Map<Integer,Float> getMinus() {return resultsMinus;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    int width = props.BinWidth;
                    boolean extension = props.ReadExtension;
                    // for GaussianKernel, get 1bp resolution data
                    if (props.GaussianKernelWidth!=0 && region.getWidth()<=1000){ 
                    	width = 1;
                    }
                    resultsPlus = null;
                    resultsMinus = null;
                    if (props.ShowSelfLigationOverlap) {
                    	System.err.println("Computing Self Overlap");
                    	resultsPlus = getSelfHistogram();
                    	resultsMinus = new TreeMap<Integer,Float>();
                    	System.err.println(resultsPlus.size()+" in model plus set");
                    } else if (props.UseWeights) {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                            resultsPlus = Aggregator.mergeHistogramsFF(resultsPlus,
                                                                       client.getWeightHistogram(ids,
                                                                                                 region.getGenome().getChromID(region.getChrom()),
                                                                                                 false,
                                                                                                 extension,
                                                                                                 width,
                                                                                                 (int)props.DeDuplicate,
                                                                                                 region.getStart(),
                                                                                                 region.getEnd(),
                                                                                                 null,
                                                                                                 true));
                            
                            resultsMinus = Aggregator.mergeHistogramsFF(resultsMinus,
                                                                        client.getWeightHistogram(ids,
                                                                                                  region.getGenome().getChromID(region.getChrom()),
                                                                                                  false,
                                                                                                  extension,
                                                                                                  width,
                                                                                                  (int)props.DeDuplicate,
                                                                                                  region.getStart(),
                                                                                                  region.getEnd(),
                                                                                                  null,
                                                                                                  false));
                        }
                        if (props.ShowPairedReads) {
                            resultsPlus = Aggregator.mergeHistogramsFF(resultsPlus,
                                                                    client.getWeightHistogram(ids,
                                                                                              region.getGenome().getChromID(region.getChrom()),
                                                                                              true,
                                                                                              extension,
                                                                                              width,
                                                                                              (int)props.DeDuplicate,
                                                                                              region.getStart(),
                                                                                              region.getEnd(),
                                                                                              null,
                                                                                              true));
                            
                            resultsMinus = Aggregator.mergeHistogramsFF(resultsMinus,
                                                                     client.getWeightHistogram(ids,
                                                                                               region.getGenome().getChromID(region.getChrom()),
                                                                                               true,
                                                                                               extension,
                                                                                               width,
                                                                                               (int)props.DeDuplicate,
                                                                                               region.getStart(),
                                                                                               region.getEnd(),
                                                                                               null,
                                                                                               false));

                        }

                    } else {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                            resultsPlus = Aggregator.mergeHistogramsIF(client.getHistogram(ids,
                                                                                           region.getGenome().getChromID(region.getChrom()),
                                                                                           false,
                                                                                           extension,
                                                                                           width,
                                                                                           (int)props.DeDuplicate,
                                                                                           region.getStart(),
                                                                                           region.getEnd(),
                                                                                           null,
                                                                                           true),
                                                                    resultsPlus);
                            
                            resultsMinus = Aggregator.mergeHistogramsIF(client.getHistogram(ids,
                                                                                            region.getGenome().getChromID(region.getChrom()),
                                                                                            false,
                                                                                            extension,
                                                                                            width,
                                                                                            (int)props.DeDuplicate,
                                                                                            region.getStart(),
                                                                                            region.getEnd(),
                                                                                            null,
                                                                                            false),
                                                                     resultsMinus);
                        }
                        if (props.ShowPairedReads) {
                            resultsPlus = Aggregator.mergeHistogramsIF(
                                                                    client.getHistogram(ids,
                                                                                        region.getGenome().getChromID(region.getChrom()),
                                                                                        true,
                                                                                        extension,
                                                                                        width,
                                                                                        (int)props.DeDuplicate,
                                                                                        region.getStart(),
                                                                                        region.getEnd(),
                                                                                        null,
                                                                                        true),
                                                                       resultsPlus);
                            
                            resultsMinus = Aggregator.mergeHistogramsIF(
                                                                     client.getHistogram(ids,
                                                                                         region.getGenome().getChromID(region.getChrom()),
                                                                                         true,
                                                                                         extension,
                                                                                         width,
                                                                                         (int)props.DeDuplicate,
                                                                                         region.getStart(),
                                                                                         region.getEnd(),
                                                                                         null,
                                                                                         false),
                                                                       resultsMinus);

                        }

                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    // assign empty output.  This is useful because Client
                    // throws an exception for non-existant chromosomes, such
                    // as those for which there were no alignment results
                    resultsPlus = new TreeMap<Integer,Float>();
                    resultsMinus = resultsPlus;
                }

                newinput = false;
                notifyListeners();
            }
        }
        client.close();
    }
    
    public TreeMap<Integer,Float> getSelfHistogram() throws IOException, ClientException {
    	int[] profile = new int[region.getWidth()];
		Set<PairedHit> hitset = getHitSet(region);
		for (PairedHit p : hitset) {
			if (p==null) System.err.println("null PairedHit");
			if (isSelfLigation(p)) {
				int startpos = p.leftPos < p.rightPos ? p.leftPos : p.rightPos;
				int endpos = p.leftPos < p.rightPos ? p.rightPos : p.leftPos;
				startpos = Math.max(startpos-region.getStart(), 0);
				endpos = Math.min(endpos-region.getStart(), region.getWidth());
				for (int i=startpos; i<endpos; i++) {
					profile[i]++;
				}
			}
		}
		TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
		for (int i=0; i<profile.length; i++) {
			output.put(i+region.getStart(), (float)profile[i]);
		}
		return output;
    }
    
    public boolean isSelfLigation(PairedHit p) {
    	if (props.RightFlipped) {
    		return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && (p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
    		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    	} else {
    		return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && !(p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
    		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    	}
    }
    
    public Set<PairedHit> getHitSet(Region r) throws IOException, ClientException {
		int chrom = r.getGenome().getChromID(r.getChrom());
		Set<PairedHit> hitset = new HashSet<PairedHit>();
		/*
		if (chrom==9848 || chrom==9853 || chrom==9845) {
			System.err.println(r.getChrom());
			return hitset;
		}
		*/
		for (String s : ids) {
			List<PairedHit> hits = client.getPairedHits(s,
					chrom,
					true,
					r.getStart()-props.SelfLigationCutoff,
					r.getEnd(),
					null,
					null);
			for (PairedHit p : hits) {
				if (!(hitset.contains(p))) {
					p.flipSides();
					if (!hitset.contains(p)) {
						p.flipSides();
						hitset.add(p);
						
					}
				}
			}
			
		}
		System.err.println(hitset.size()+" hits");
		return hitset;
	}
    
 }