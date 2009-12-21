package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;
import java.util.LinkedList;

import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.SimpleDomain;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class SimpleDomainFinder 
	implements Expander<Region,SimpleDomain> {
	
	private ChipChipData data;
	private int minSize;
	private int maxJump;
	private double meanThreshold;
	private double probeThreshold;
	
	public SimpleDomainFinder(ChipChipData ccd) { 
		data = ccd;
		minSize = 1000;
		meanThreshold = 3.0;
		maxJump = 500;
		probeThreshold = 0.75;
	}

	public void setMeanThreshold(double mt) { meanThreshold = mt; }
	
	public double getMeanThreshold() { return meanThreshold; }
	public int getMaxJump() { return maxJump; }
	public int getMinSize() { return minSize; }

	public Iterator<SimpleDomain> execute(Region a) {
		try {
			LinkedList<SimpleDomain> doms = new LinkedList<SimpleDomain>();
			data.window(a.getChrom(), a.getStart(), a.getEnd());
			int np = data.getCount();
			
			for(int pi = 0; pi < np-1; pi++) { 
				int f = findNextDomainIndex(pi, a);
				if(f != -1) { 
					SimpleDomain dom = buildDomain(pi, f, a);
					doms.addLast(dom);
					pi = f;
				}
			}
			
			return doms.iterator();
		} catch (NotFoundException e) {
			e.printStackTrace();
			return new EmptyIterator<SimpleDomain>();
		}
	}
	
	private int findNextDomainIndex(int pi, Region w) { 
		double sum = getRatio(pi);
		int farthest = pi;
		int start = data.getPos(pi);
		int lastHighProbe = pi;
		int probeCount = 1;
		
		if(sum >= meanThreshold) { 
			boolean searching = true;
			int ppos = start;
			
			for(int j = pi + 1; searching && j < data.getCount(); j++) { 
				double r = getRatio(j);
				int jpos = data.getPos(j);
				
				if(jpos - ppos <= maxJump) { 
					sum += r;

					if(r >= meanThreshold) { 
						//lastHighProbe = j; 
						probeCount += 1;
					
						double mean = sum / (double)(j-pi+1);
						int numProbes = j-pi+1;
						double probeThreshFrac = (double)probeCount / (double)numProbes;

						if(mean >= meanThreshold && probeThreshFrac >= probeThreshold) { 
							farthest = j;
						}
					}

				} else { 
					searching = false;
				}
				
				ppos = jpos;
			}
			
			//farthest = lastHighProbe;
		} else { 
			return -1;
		}
		
		if(farthest > pi && 
			data.getPos(farthest) - data.getPos(pi)
				> minSize) { 
			return farthest;
		} else { 
			return -1;
		}
	}
	
	private int getLeftBound(int i, Region w) { 
		int left = i == 0 ? w.getStart() : 
			Math.max(Math.max(1, data.getPos(i) - maxJump), (data.getPos(i) + data.getPos(i-1)) / 2);
		return left;
	}
	
	private int getRightBound(int i, Region w) { 
		int right = i >= data.getCount()-1 ? w.getEnd() : 
			Math.min(data.getPos(i) + maxJump, (data.getPos(i) + data.getPos(i+1))/2);
		return right;
	}
	
	private SimpleDomain buildDomain(int pi, int farthest, Region w) { 
		int left = getLeftBound(pi, w);
		int right = getRightBound(farthest, w);
		return new SimpleDomain(w.getGenome(), w.getChrom(), left, right);
	}
	
	private double getRatio(int i) { 
		double sum = 0.0;
		int c = 0;
		for(int j = 0; j < data.getReplicates(i); j++) { 
			sum += data.getRatio(i, j);
			c += 1;
		}
		return c > 0 ? sum / (double)c : 0.0; 
	}
}
