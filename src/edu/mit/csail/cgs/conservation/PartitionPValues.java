/*
 * Created on Sep 18, 2006
 */
package edu.mit.csail.cgs.conservation;

import java.io.*;

/**
 * @author tdanford
 */
public class PartitionPValues {
	
	public static void main(String[] args) { 
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		String line = null;
		//ConservationPredicate pred = new ConservationPredicate(0.0);
		
		try {
			System.out.print(">"); System.out.flush();
			while((line = br.readLine()) != null) { 
				String[] array = line.split("\\s+");
				if(array.length >= 2) { 
					int[] sizes = new int[array.length-1];
					int numTrue = Integer.parseInt(array[0]);
					for(int i = 0; i < sizes.length; i++) { sizes[i] = Integer.parseInt(array[i+1]); }
					ThresholdPredicate pred = new ThresholdPredicate(sizes, 0.6, 0.4);
					Results r = PartitionPValues.countArrangements(pred, sizes, numTrue);
					double frac = r.getFraction();
					System.out.println(r.getAccepted() + "/" + r.getTotal() + ": " + frac);
				}
				System.out.print(">"); System.out.flush();
			}
		} catch(IOException ie ) { 
			ie.printStackTrace(System.err);
		}
	}

    public PartitionPValues() { 
        
    }
    
    public static interface Predicate { 
        public boolean accepts(int[] counts, int[] sizes, int numTrue);
    }
    
    public static class ThresholdPredicate implements Predicate { 
    	private int[] threshs;
    	
    	public ThresholdPredicate(int[] sizes, double uthresh, double dthresh) {
    		threshs = new int[sizes.length];
    		for(int i = 0; i < threshs.length; i++) { 
    			if(i <= 1 ) { 
    				threshs[i] = (int)Math.ceil((double)sizes[i] * uthresh);
    			} else { 
    				threshs[i] = (int)Math.floor((double)sizes[i] * dthresh);
    			}
    		}
    	}

		public boolean accepts(int[] counts, int[] sizes, int numTrue) {
			if(counts[0] < threshs[0] || counts[1] < threshs[1]) { return false; }
			if(counts[2] > threshs[2] || counts[3] > threshs[3]) { return false; }
			return true;
		}
    	
    	
    }
    
    public static class ConservationPredicate implements Predicate {
        
        private double fracSeparation; 
        
        public ConservationPredicate(double fs) {
            fracSeparation = fs;
        }

        public boolean accepts(int[] counts, int[] sizes, int numTrue) {
            double f0 = (double)counts[0] / (double)sizes[0];
            double f1 = (double)counts[1] / (double)sizes[1];
            double f2 = (double)counts[2] / (double)sizes[2];
            double f3 = (double)counts[3] / (double)sizes[3];
            double fmin = Math.min(f0, f1), fmax = Math.max(f2, f3);
            return fmin - fmax > fracSeparation;
        } 
    }
    
    public static class Results {
       
        private long accepted, total;
        
        public Results() { 
            accepted = total = 0;
        }
        
        public long getAccepted() { return accepted; }
        public long getTotal() { return total; }
        
        public double getFraction() { 
            return (double)accepted / (double)total;
        }
        
        public void addAccept() { 
            accepted += 1;
            total += 1;
            
            if(accepted == Long.MAX_VALUE || total == Long.MAX_VALUE) { 
                throw new RuntimeException("Overflow Condition Reached.");
            }
        }
        
        public void addReject() { 
            total += 1;
            
            if(total == Long.MAX_VALUE) { 
                throw new RuntimeException("Overflow Condition Reached.");                
            }
        }
    }
    
    public static Results countArrangements(Predicate p, int[] sizes, int numTrue) { 
        Results r = new Results();
        int[] counts = new int[sizes.length];
        for(int i = 0; i < counts.length; i++) { counts[i] = 0; }
        recCountArrangements(p, sizes, numTrue, numTrue, counts, 0, r);
        return r;
    }
    
    public static void printArray(int[] array, PrintStream ps) { 
    	for(int i = 0; i < array.length; i++) { 
    		ps.print(array[i]);
    		if(i < array.length-1) { ps.print(" "); }
    	}
    }
    
    private static void recCountArrangements(Predicate p, int[] sizes, int numTrue, int trueLeft, 
            int[] counts, int pos, Results r) {
        
        if(pos >= counts.length) {
            if(trueLeft == 0) {
                if(p.accepts(counts, sizes, numTrue)) { 
                	//printArray(counts, System.out); System.out.println();
                    r.addAccept();
                } else { 
                    r.addReject();
                    //System.out.print(" .");
                }
                
            } else { 
            	//System.out.print(" --> None");
            }
            
            //System.out.println();
        } else { 
            for(int i = 0; i <= Math.min(sizes[pos], trueLeft); i++) { 
                counts[pos] = i;
                recCountArrangements(p, sizes, numTrue, trueLeft-i, counts, pos+1, r);
            }
        }
    }
}
