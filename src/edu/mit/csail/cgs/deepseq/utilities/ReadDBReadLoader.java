package edu.mit.csail.cgs.deepseq.utilities;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqExpt;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLoader;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.ExtReadHit;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.projects.readdb.*;

public class ReadDBReadLoader extends ReadLoader{
	
	private Client client=null;
	private List<String> exptNames =new ArrayList<String>();
	private List<ChipSeqAlignment> aligns = new ArrayList<ChipSeqAlignment>();
	protected int currID=0;
	
	public ReadDBReadLoader(Genome g, List<ChipSeqLocator> locs, int rLen){
		super(g, rLen);
		
		currID=new Random(System.currentTimeMillis()).nextInt();
		try {
			client = new Client();
		
			//Initialize
            ChipSeqLoader loader = new ChipSeqLoader(false); 
			for(ChipSeqLocator locator : locs){
				String exptName = locator.getExptName(); exptNames.add(exptName);
				if (locator.getAlignName() == null) {
		            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
		            	Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
		        		for(ChipSeqExpt expt : expts) { 
		                	Collection<ChipSeqAlignment> aligns;
							aligns = loader.loadAllAlignments(expt);
							for (ChipSeqAlignment currentAlign : aligns) {
		            			if (currentAlign.getGenome().equals(g)) { 
		            				aligns.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		            } else { //No alignment name, given replicate names
		                for(String repName : locator.getReplicates()) { 
		                    ChipSeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
		                    ChipSeqAlignment alignment = 
		                        loader.loadAlignment(expt, locator.getAlignName(), g);
		                    if(alignment != null) { 
		                        aligns.add(alignment);
		                        break;
		                    }
		                }
		            }
		        } else {
		        	if(locator.getReplicates().isEmpty()) {//Given alignment name, no replicate names
		        		Collection<ChipSeqExpt> expts = loader.loadExperiments(locator.getExptName());
		        		for(ChipSeqExpt expt : expts) { 
		                	Collection<ChipSeqAlignment> alignments;
							alignments = loader.loadAllAlignments(expt);
							for (ChipSeqAlignment currentAlign : alignments) {
		            			if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) { 
		            				aligns.add(currentAlign);
		    						break;
		    					}
		            		}
		    			}
		            }else{
		            	for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
		        			aligns.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
                                                                                  replicate), 
                                                            locator.getAlignName(),
                                                            g));
		        		}
		            }
		        }
			}countHits();
			
			//Error that doesn't seem to be caught by the exceptions
            //			if(totalHits==0){
            //				System.err.println("No reads found for these experiment names");
            //				System.exit(1);
            //			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	protected double countHits(){
		totalHits=0;
		totalWeight=0;
		try {
			for(ChipSeqAlignment alignment : aligns) { 
				double currHits = (double)client.getCount(Integer.toString(alignment.getDBID()), false, null, null);
				totalHits+=currHits;
				double currWeight =(double)client.getWeight(Integer.toString(alignment.getDBID()), false, null, null);
				totalWeight +=currWeight;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
			return(totalHits);
		}			
		return totalHits;
	}
	
	//Count total read hits on one strand 
	protected double countStrandedWeight(char strand) {
		double strandWeight=0;
		try {
			for(ChipSeqAlignment alignment : aligns) {
                strandWeight += (double)client.getWeight(Integer.toString(alignment.getDBID()), 
                                                         false, false, null);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
			return(strandWeight);
		}
        return strandWeight;
	}
	
	//Load reads in a region
	public List<ReadHit> loadHits(Region r) {
		ArrayList<ReadHit> total = new ArrayList<ReadHit>();
		try {
			for(ChipSeqAlignment alignment : aligns) {
                for (SingleHit h : client.getSingleHits(Integer.toString(alignment.getDBID()),
                                                        r.getGenome().getChromID(r.getChrom()),
                                                        r.getStart(),
                                                        r.getEnd(),
                                                        null,
                                                        null)) {
                    total.add(convertToReadHit(r.getGenome(),
                                               currID++,
                                               h));
                }
            }
			return total;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
			return(total);
		}
		return total;
	}
	
    // load paired read hit coordinates (sorted) and counts
    public Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedBaseCounts(Region r, char strand){
        Collection<String> alignids = new ArrayList<String>();
        for(ChipSeqAlignment alignment : aligns) {
            alignids.add(Integer.toString(alignment.getDBID()));
        }
        TreeMap<Integer,Float> allHits = null;
        ArrayList<Integer> coords = new ArrayList<Integer>();
        ArrayList<Float> counts = new ArrayList<Float>();
        try {
            allHits = client.getWeightHistogram(alignids,
                                                r.getGenome().getChromID(r.getChrom()),
                                                false,
                                                false,
                                                1,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                strand == '+');
            coords.addAll(allHits.keySet());
            counts.addAll(allHits.values());
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ClientException e) {
            //Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
        }
        return new Pair<ArrayList<Integer>,ArrayList<Float>>(coords, counts);
    }
	
	// get BED-format reads
	// each hit count as 1 read, ignore the weights
	public String getBED_StrandedReads(Region r, char strand, double probability){
		DRand random = new DRand();
		String head="chr"+r.getChrom()+"\t";
		String tail="\tU\t0\t"+strand;
		StringBuilder sb = new StringBuilder();
		int start = 0;
		int end = 0;
		try {
			for(ChipSeqAlignment alignment : aligns) {
                for (SingleHit h : client.getSingleHits(Integer.toString(alignment.getDBID()),
                                                        r.getGenome().getChromID(r.getChrom()),
                                                        r.getStart(),
                                                        r.getEnd(),
                                                        null,
                                                        strand == '+')) {
                    sb.append(head).append(h.pos).append("\t").append(strand == '+' ? (h.pos+h.length) : (h.pos - h.length)).append(tail).append("\n");                    
                }
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
		}
		return sb.toString();
	}
	//Load extended reads in a region
	public List<ExtReadHit> loadExtHits(Region r, int startShift, int fivePrimeExt, int threePrimeExt) {
		ArrayList<ExtReadHit> total = new ArrayList<ExtReadHit>();
		try {
			for(ChipSeqAlignment alignment : aligns) {
                for (SingleHit h : client.getSingleHits(Integer.toString(alignment.getDBID()),
                                                        r.getGenome().getChromID(r.getChrom()),
                                                        r.getStart(),
                                                        r.getEnd(),
                                                        null,
                                                        null)) {
                    total.add(convertToExtReadHit(r.getGenome(),
                                                  currID++,
                                                  h,
                                                  startShift, fivePrimeExt, threePrimeExt));
                }
			}
        } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
		}
		return total;	
	}

	// count number of reads in region
	public int countHits (Region r) {
        int count = 0;
		try {
			for(ChipSeqAlignment alignment : aligns) { 
                count += client.getCount(Integer.toString(alignment.getDBID()),
                                         r.getGenome().getChromID(r.getChrom()),
                                         false,
                                         r.getStart(),
                                         r.getEnd(),
                                         null,
                                         null,
                                         null);
                
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
		}
        return count;
	}// sum weight of reads in region
	public double sumWeights (Region r) {
        double sum = 0;
		try {
			for(ChipSeqAlignment alignment : aligns) { 
                sum += client.getWeight(Integer.toString(alignment.getDBID()),
                                        r.getGenome().getChromID(r.getChrom()),
                                        false,
                                        r.getStart(),
                                        r.getEnd(),
                                        null,
                                        null,
                                        null);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
		}
        return sum;
	}

    public ReadHit convertToReadHit(Genome g, int id, SingleHit h) {
        return new ReadHit(g,
                           id,
                           g.getChromName(h.chrom),
                           h.pos,
                           h.pos + h.length,
                           h.strand ? '+' : '-',
                           h.weight);
    }
    public ExtReadHit convertToExtReadHit(Genome g, int id, SingleHit h, int startshift, int fiveprime, int threeprime) {
        return new ExtReadHit(g,
                              id,
                              g.getChromName(h.chrom),
                              h.pos,
                              h.pos + h.length,
                              h.strand ? '+' : '-',
                              h.weight,
                              startshift, fiveprime, threeprime);
    }

	
		//Convert ChipSeqHits to ExtReads
	private Collection<ExtReadHit> weights2extreads(Region reg, int [] hits, float [] weights,char strand, int startShift, int fivePrimeExt, int threePrimeExt){
		ArrayList<ExtReadHit> r = new ArrayList<ExtReadHit>();
		for(int i = 0; i< hits.length; i++){
			currID++;
			int start = strand=='+' ? hits[i] : hits[i]-readLength+1;
			int end = strand=='+' ? hits[i]+readLength-1 : hits[i];
			r.add(new ExtReadHit(reg.getGenome(), currID, reg.getChrom(), start, end, strand, weights[i], startShift, fivePrimeExt, threePrimeExt));
		}
        return(r);
	}

	public void setGenome(Genome g) {
		gen=g;
	}
	
	//Close the loaders
	public void cleanup(){
		if(client!=null){
			client.close();
		}
	}
}