package edu.mit.csail.cgs.warpdrive.model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.ProfileRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.chiapet.Neighborhood;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * Data model for chipseq histogram.  Separate methods for retrieving
 * plus and minus strand results
 */
public class ChipSeqHistogramModel extends WarpModel implements RegionModel, Runnable {
    
    private Client client;
    private TreeMap<Integer,Float> resultsPlus, resultsMinus, resultsPval;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    private ChipSeqHistogramProperties props;

    private Region region;
    private boolean newinput;
    
    private double[] kernel;
    private double[] readdist;
	private double[] eventdist;
	private double[] condist;
	private int cutoff;
	private Map<Integer,String> revChromMap;
	private Poisson poisson = new Poisson(1, new DRand());

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
        resultsPval = null;
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
    public Map<Integer,Float> getPval() {return resultsPval;}
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
                    resultsPval = null;
                    if (props.ShowSelfLigationOverlap) {
                    	resultsPlus = getSelfHistogram();
                    	resultsMinus = new TreeMap<Integer,Float>();
                    } 
                    else if (props.UseWeights) {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                        	try{
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
                        	}catch (Exception ex) {
                                //Fail silently if there are no single read alignments
                            }
                        }
                        if (props.ShowPairedReads) {
                        	try{
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
                        	}catch (Exception ex) {
                                //Fail silently if there are no paired read alignments
                            }
                        }

                    } else {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                        	try{
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
                        	}catch (Exception ex) {
                                //Fail silently if there are no single read alignments
                            }
                        }
                        if (props.ShowPairedReads) {
                        	try{
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
                        	}catch (Exception ex) {
                                //Fail silently if there are no paired read alignments
                            }
                        }
                    }
                } catch (Exception ex) {
                    //ex.printStackTrace();
                }
                if(resultsPlus==null || resultsMinus==null){
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
    
    public static double[] readDoubleList(String file) throws IOException {
		List<Double> list = new ArrayList<Double>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			list.add(Double.parseDouble(split[0]));
		}
		double[] tor = new double[list.size()];
		for (int i=0; i<tor.length; i++) {
			tor[i] = list.get(i);
		}
		return tor;
	}
    
    /*
    public TreeMap<Integer,Float> getInteractionProbability() throws IOException, ClientException {
    	resultsPval = new TreeMap<Integer,Float>();
    	Region anchor = Region.fromString(region.getGenome(), props.Anchor);
    	int anchorchrom = anchor.getGenome().getChromID(anchor.getChrom());
    	int regionchrom = region.getGenome().getChromID(region.getChrom());
    	List<PairedHit> hits = new ArrayList<PairedHit>();
    	for (String s : ids) {
			hits.addAll(client.getPairedHits(s,
					anchorchrom,
					true,
					anchor.getStart(),
					anchor.getEnd(),
					null,
					true));
			hits.addAll(client.getPairedHits(s,
					anchorchrom,
					false,
					anchor.getStart(),
					anchor.getEnd(),
					null,
					true));
			hits.addAll(client.getPairedHits(s,
					anchorchrom,
					true,
					anchor.getStart(),
					anchor.getEnd(),
					null,
					false));
			hits.addAll(client.getPairedHits(s,
					anchorchrom,
					false,
					anchor.getStart(),
					anchor.getEnd(),
					null,
					false));
    	}
    	int[] probs = new int[region.getWidth()/props.BinWidth+1];
    	for (PairedHit p : hits) {
    		if (p.leftChrom==anchorchrom && anchor.getStart()<=p.leftPos && anchor.getEnd()>=p.leftPos) {
    			if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
    				probs[(p.rightPos-region.getStart())/props.BinWidth]++;
    			}
    		} else if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
    			probs[(p.leftPos-region.getStart())/props.BinWidth]++;
    		}
    	}
    	TreeMap<Integer,Float> tor = new TreeMap<Integer,Float>();
    	System.err.println(probs.length+" bins");
    	for (int i=0; i<probs.length; i++) {
    		if (probs[i]>0) {
    			tor.put(region.getStart()+i*props.BinWidth, (float)probs[i]);
    			Region probe = new Region(region.getGenome(),region.getChrom(),(region.getStart()+i*props.BinWidth),(region.getStart()+(i+1)*props.BinWidth));
    			//System.err.println(probe+"\t"+interactionPValue(anchor, probe, probs[i]));
    			float pval = (float)interactionPValue(anchor, probe, probs[i]);
    			if (pval <= (props.PValueCutoff/((float)probs.length))) {
    				resultsPval.put(probe.getStart(), pval);
    			}
    		}
    	}
    	return tor;
    }
    */
    
    /*
    private double interactionPValue(Region r1, Region r2, int count) throws IOException, ClientException {
    	int count1 = 0;
		int chrom1 = r1.getGenome().getChromID(r1.getChrom());
		for (String a : ids) {
			count1 += client.getCount(a, chrom1, true, r1.getStart(), r1.getEnd(), null, true, true);
			count1 += client.getCount(a, chrom1, true, r1.getStart(), r1.getEnd(), null, true, false);
			count1 += client.getCount(a, chrom1, true, r1.getStart(), r1.getEnd(), null, false, true);
			count1 += client.getCount(a, chrom1, true, r1.getStart(), r1.getEnd(), null, false, false);
		}
		int count2 = 0;
		int chrom2 = r2.getGenome().getChromID(r2.getChrom());
		for (String a : ids) {
			count2 += client.getCount(a, chrom2, true, r2.getStart(), r2.getEnd(), null, true, true);
			count2 += client.getCount(a, chrom2, true, r2.getStart(), r2.getEnd(), null, true, false);
			count2 += client.getCount(a, chrom2, true, r2.getStart(), r2.getEnd(), null, false, true);
			count2 += client.getCount(a, chrom2, true, r2.getStart(), r2.getEnd(), null, false, false);
		}
		double prob1 = ((double)count1) / ((double)props.TotalReads);
		double prob2 = ((double)count2) / ((double)props.TotalReads);
		poisson.setMean(props.ChimericReads*prob1*prob2);
		return 1.0d - poisson.cdf(count) + poisson.pdf(count);
    }
    
    public TreeMap<Integer,Float> getForwardInteractionKernel() throws IOException, ClientException {
    	if (kernel==null) {
    		kernel = readDoubleList(props.Kernel);
    	}
    	TreeMap<Integer,Float> tor = new TreeMap<Integer,Float>();
    	double[] arr = getForwardProfile(region);
    	double max = 0;
    	for (int i=0; i<arr.length; i++) {
    		if (arr[i]>max) max = arr[i];
    	}
    	for (int i=0; i<arr.length; i++) {
    		arr[i] = arr[i]*100d/max;
    		tor.put(region.getStart()+i, (float)arr[i]);
    	}
    	return tor;
    }
    
    public TreeMap<Integer,Float> getReverseInteractionKernel() throws IOException, ClientException {
    	if (kernel==null) {
    		kernel = readDoubleList(props.Kernel);
    	}
    	TreeMap<Integer,Float> tor = new TreeMap<Integer,Float>();
    	double[] arr = getReverseProfile(region);
    	double max = 0;
    	for (int i=0; i<arr.length; i++) {
    		if (arr[i]>max) max = arr[i];
    	}
    	for (int i=0; i<arr.length; i++) {
    		arr[i] = arr[i]*100d/max;
    		tor.put(region.getStart()+i, (float)arr[i]);
    	}
    	return tor;
    }
    
    public double[] getForwardProfile(Region r) throws IOException, ClientException {
    	int halfkern = kernel.length/2;
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.length, kernel.length);
		Set<PairedHit> subset = getHitSet(e);
		for (PairedHit p : subset) {
			if (p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel[i-pmkern];
				}
			}
			if (p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel[i-pmkern];
				}
			}
		}
		return tor;
	}
	*/
    
    public double[] getReverseProfile(Region r) throws IOException, ClientException {
    	int halfkern = kernel.length/2;
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.length, kernel.length);
		Set<PairedHit> subset = getHitSet(e);
		for (PairedHit p : subset) {
			if (!p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel[i-pmkern];
				}
			}
			if (!p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel[i-pmkern];
				}
			}
		}
		return tor;
	}
    
    /*
    public TreeMap<Integer,Float> getInteractionProfile() throws IOException, ClientException {
    	if (readdist==null) {
    		double[] readdista = readDoubleList(props.ReadDistribution);
    		if (readdista.length % 2 == 0) {
    			double[] tmparr = readdista;
    			readdista = new double[tmparr.length-1];
    			for (int i=0; i<readdista.length; i++) {
    				readdista[i] = tmparr[i];
    			}
    		}
    		double[] eventdista = readDoubleList(props.EventDistribution);
    		if (eventdista.length % 2 == 0) {
    			double[] tmparr = eventdista;
    			eventdista = new double[tmparr.length-1];
    			for (int i=0; i<eventdista.length; i++) {
    				eventdista[i] = tmparr[i];
    			}
    		}
    		this.readdist = readdista;
    		this.eventdist = eventdista;
    		this.cutoff = eventdist.length/2;
    		revChromMap = region.getGenome().getRevChromIDMap();
    	}
    	Neighborhood n = new Neighborhood();
    	Point tss = Point.fromString(region.getGenome(), props.TSS);
		Region tssregion = new Region(region.getGenome(), tss.getChrom(), tss.getLocation()-cutoff, tss.getLocation()+cutoff);
		Set<PairedHit> hitset = getHitSet(tssregion);
		Region rextend = region.expand(cutoff, cutoff);
		Region rexex = rextend.expand(cutoff, cutoff);
		double[] exprof = new double[rexex.getEnd()-rexex.getStart()+1];
		n.addProfile(new ProfileRegion(region.getGenome(), rexex.getChrom(), rexex.getStart(), rexex.getEnd(), exprof));
		for (PairedHit p : hitset) {
			double leftscore = 0;
			double rightscore = 0;
			String leftChrom = revChromMap.get(p.leftChrom);
			String rightChrom = revChromMap.get(p.rightChrom);
			if ((leftChrom.equals(tss.getChrom())) && (Math.abs(p.leftPos-tss.getLocation())<cutoff)) {
				leftscore = eventdist[p.leftPos-tss.getLocation()+cutoff];
			}
			if ((rightChrom.equals(tss.getChrom())) && (Math.abs(p.rightPos-tss.getLocation())<cutoff)) {
				rightscore = eventdist[p.rightPos-tss.getLocation()+cutoff];
			}
			if (leftscore>rightscore) {
				if (rextend.getChrom().equals(rightChrom) && p.rightPos>=rextend.getStart() && p.rightPos<=rextend.getEnd()) {
					double[] profile = new double[readdist.length];
					if (p.rightStrand) {
						for (int i=0; i<profile.length; i++) {
							profile[i] = readdist[i]*leftscore;
						}
					} else {
						for (int i=0; i<profile.length; i++) {
							profile[i] = readdist[readdist.length-i-1]*leftscore;
						}
					}
					ProfileRegion pr = new ProfileRegion(region.getGenome(), rightChrom, p.rightPos-cutoff, p.rightPos+cutoff, profile);
					
					n.addProfile(pr);
				}
			} else {
				if (rextend.getChrom().equals(leftChrom) && p.leftPos>=rextend.getStart() && p.leftPos<=rextend.getEnd()) {
					double[] profile = new double[readdist.length];
					if (p.leftStrand) {
						for (int i=0; i<profile.length; i++) {
							profile[i] = readdist[i]*rightscore;
						}
					} else {
						for (int i=0; i<profile.length; i++) {
							profile[i] = readdist[readdist.length-i-1]*rightscore;
						}
					}
					ProfileRegion pr = new ProfileRegion(region.getGenome(), leftChrom, p.leftPos-cutoff, p.leftPos+cutoff, profile);
					//System.err.println("max query: "+(pr.getEnd()-pr.getStart()));
					//System.err.println("array size: "+profile.length);
					n.addProfile(pr);
				}
			}
		}
		TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
		float[] profile = n.getProfile();
		int dcutoff = cutoff*2;
		for (int i=dcutoff; i<dcutoff+region.getWidth(); i++) {
			output.put(i-dcutoff+region.getStart(), profile[i]*300000000);
		}
		return output;
    }
    */
    
    public TreeMap<Integer,Float> getSelfHistogram() throws IOException, ClientException {
    	int[] profile = new int[region.getWidth()];
		int shift = props.SelfLigationCutoff+props.SmoothingWindowWidth;
		int halfWindow = props.SmoothingWindowWidth/2;
		Region expanded = region.expand(shift, shift);
		int[] eprofile = new int[expanded.getWidth()];
		Set<PairedHit> hitset = getHitSet(expanded);
		for (PairedHit p : hitset) {
			if (isSelfLigation(p)) {
				int startpos = p.leftPos < p.rightPos ? p.leftPos : p.rightPos;
				int endpos = p.leftPos < p.rightPos ? p.rightPos : p.leftPos;
				startpos = Math.max(startpos-expanded.getStart(), 0);
				endpos = Math.min(endpos-expanded.getStart(), expanded.getWidth());
				for (int i=startpos; i<endpos; i++) {
					eprofile[i]++;
				}
			}
		}
		for (int i=0; i<profile.length; i++) {
			int sum = 0;
			for (int j=shift+i-halfWindow; j<shift+i+1+halfWindow; j++) {
				sum += eprofile[j];
			}
			profile[i] = sum / (1 + 2*halfWindow);
		}
		TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
		for (int i=0; i<profile.length; i++) {
			output.put(i+region.getStart(), (float)profile[i]);
		}
		return output;
    }
    
    public TreeMap<Integer,Float> getOrientedHistogram(boolean left, boolean right) throws IOException, ClientException {
    	int count = 0;
    	int[] profile = new int[region.getWidth()];
		Set<PairedHit> hitset = getHitSet(region);
		for (PairedHit p : hitset) {
			if (p==null) System.err.println("null PairedHit");
			if ((left && right && isPlusPlus(p)) ||
					(left && !right && isPlusMinus(p)) ||
					(!left && right && isMinusPlus(p)) ||
					(!left && !right && isMinusMinus(p))) {
				int startpos = p.leftPos < p.rightPos ? p.leftPos : p.rightPos;
				int endpos = p.leftPos < p.rightPos ? p.rightPos : p.leftPos;
				startpos = Math.max(startpos-region.getStart(), 0);
				endpos = Math.min(endpos-region.getStart(), region.getWidth());
				for (int i=startpos; i<endpos; i++) {
					profile[i]++;
				}
				count++;
			}
		}
		TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
		for (int i=0; i<profile.length; i++) {
			output.put(i+region.getStart(), (float)profile[i]);
		}
		System.out.println(left+"\t"+right+"\t"+count);
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
    
    public boolean isPlusMinus(PairedHit p) {
    	return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && (p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
		&& !(p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    }
    
    public boolean isMinusPlus(PairedHit p) {
    	return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && !(p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    }
    
    public boolean isPlusPlus(PairedHit p) {
    	return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && (p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
		&& (p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
    }
    
    public boolean isMinusMinus(PairedHit p) {
    	return (p.leftChrom == p.rightChrom) && (Math.abs(p.leftPos-p.rightPos) <= props.SelfLigationCutoff) && !(p.leftPos < p.rightPos ? p.leftStrand : p.rightStrand)
		&& !(p.leftPos < p.rightPos ? p.rightStrand : p.leftStrand);
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
			
			hits = client.getPairedHits(s,
					chrom,
					false,
					r.getStart(),
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