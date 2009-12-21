/*
 * Created on Feb 11, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;
import java.io.*;
import java.text.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class BoundaryDataset implements Saveable {

	public static void main(String[] args) { 
		NumberFormat nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(5);
		try { 
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			System.out.print(">"); System.out.flush();
			BoundaryPValues bpv = new BoundaryPValues();
			double thresh = Math.log(Double.parseDouble(args[0]));
			//BoundaryTester tester = new TrialReplacementBoundaryTester(10, 10, 4, new PValueBoundaryTester(bpv, thresh));
			String line;
			while((line = br.readLine()) != null && line.length() > 0) { 
				line = line.trim();
				BoundaryDataset ds = new BoundaryDataset(line);
				ds.findBoundary();
				int errors = ds.getError();
				int N = ds.size();
				int pos = ds.getNumPositive();


				System.out.println(ds.toString());
				System.out.println("N: " + N + ", pos: " + pos + ", err: " + errors);
				if(errors > Math.min(pos, N - pos)) { continue; }
				double logPValue = bpv.getLogPValue(N, pos, errors);
				System.out.println("logPvalue: " + nf.format(logPValue) + " --> " + nf.format(Math.exp(logPValue)));
				System.out.println();

                /*
				if(tester.testDataset(ds)) { 
					System.out.println("==> Accepted"); 
				} else { 
					System.out.println("==> Rejected"); 
				}
                */

				System.out.print(">"); System.out.flush();
			}
			
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}
    
    private static Random rand;
    
    static { 
        rand = new Random();
    }
    
    private int numPos, numNeg;
	private String stringRep;
    private TreeMap<Double,Integer> positive, negative;
    private Double boundary;
    private int boundDir;
    private int error;

    public BoundaryDataset(DataInputStream dis) throws IOException { 
        positive = new TreeMap<Double,Integer>();
        negative = new TreeMap<Double,Integer>();

        numPos = dis.readInt();
        numNeg = dis.readInt();
        stringRep = dis.readUTF();
        
        int num = dis.readInt();
        for(int i = 0; i < num; i++) { 
            double key = dis.readDouble();
            int value = dis.readInt();
            positive.put(key, value);
        }

        num = dis.readInt();
        for(int i = 0; i < num; i++) { 
            double key = dis.readDouble();
            int value = dis.readInt();
            positive.put(key, value);
        }
        
        boundary = dis.readDouble();
        boundDir = dis.readInt();
        error = dis.readInt();
    }

	public BoundaryDataset(String bitString) { 
		numPos = numNeg = 0;
		positive = new TreeMap<Double,Integer>();
		negative = new TreeMap<Double,Integer>();
		boundary = null;
		boundDir = 0;
		error = -1;

		double value = 0.0;
		for(int i = 0; i < bitString.length(); i++) { 
			char c = bitString.charAt(i);
			switch(c) { 
				case '1':
					positive.put(value, 1);
					numPos += 1;
					break;
				case '0':
					negative.put(value, 1);
					numNeg += 1;
					break;
				default:
					throw new IllegalArgumentException(String.valueOf(c));
			}

			value += 1.0;
		}

		assembleStringRep();
	}

    public BoundaryDataset(Collection<Double> p, Collection<Double> n) {
        numPos = p.size();
        numNeg = n.size();
        positive = new TreeMap<Double,Integer>();
        negative = new TreeMap<Double,Integer>();
        
        for(double pv : p) { 
            if(positive.containsKey(pv)) { 
                positive.put(pv, positive.get(pv) + 1);
            } else { 
                positive.put(pv, 1);
            }
        }
        
        for(double nv : n) { 
            if(negative.containsKey(nv)) { 
                negative.put(nv, negative.get(nv) + 1);
            } else { 
                negative.put(nv, 1);
            }
        }
        
        boundary = null;
        error = -1;
        boundDir = 0;
		assembleStringRep();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.psrg.utils.Saveable#save(java.io.DataOutputStream)
     */
    public void save(DataOutputStream dos) throws IOException {
        dos.writeInt(numPos);
        dos.writeInt(numNeg);
        dos.writeUTF(stringRep);
        
        dos.writeInt(positive.size());
        for(double key : positive.keySet()) { 
            dos.writeDouble(key);
            dos.writeInt(positive.get(key));
        }

        dos.writeInt(negative.size());
        for(double key : negative.keySet()) { 
            dos.writeDouble(key);
            dos.writeInt(negative.get(key));
        }
        
        dos.writeDouble(boundary);
        dos.writeInt(boundDir);
        dos.writeInt(error);
    }
    
	public String toString() { return stringRep; }

    public int size() { return numPos + numNeg; }
    public int getNumPositive() { return numPos; }
    public int getNumNegative() { return numNeg; }
    public double getBoundary() { return boundary; }
    public int getBoundaryDirection() { return boundDir; }
    public boolean hasBoundary() { return boundary != null; }
    public int getError() { return error; }
    
    public double[] getPositiveArray() { 
        double[] array = new double[numPos];
        int i = 0;
        for(Double v : positive.keySet()) { 
            for(int j = 0; j < positive.get(v); j++) { 
                array[i++] = v;
            }
        }
        return array;
    }

    public double[] getNegativeArray() { 
        double[] array = new double[numNeg];
        int i = 0;
        for(Double v : negative.keySet()) { 
            for(int j = 0; j < negative.get(v); j++) { 
                array[i++] = v;
            }
        }
        return array;
    }
    
    public BoundaryDataset subsampleDatasetNegatives(int newNumNeg, boolean replacement) { 
        LinkedList<Double> posValues = new LinkedList<Double>(positive.keySet());
        LinkedList<Double> negValues = new LinkedList<Double>();
        
        if(replacement) { 
            for(int i = 0; i < newNumNeg; i++) { 
                negValues.addLast(chooseRandomNegativeValue());
            }
        } else { 
            double[] na = chooseRandomRemovalNegative(newNumNeg);
            for(int i = 0; i < na.length; i++) { negValues.addLast(na[i]); }
        }
        return new BoundaryDataset(posValues, negValues);        
    }
    
    public BoundaryDataset skipsampleDatasetNegatives(int offset, int skip) { 
        LinkedList<Double> posValues = new LinkedList<Double>(positive.keySet());
        LinkedList<Double> negValues = new LinkedList<Double>();

        int negOffset = 0;
        for(double val : negative.keySet()) { 
            for(int i = 0; i < negative.get(val); i++) { 
                if((negOffset-offset) % skip == 0) { negValues.addLast(val); }
                negOffset++;
            }
        }
        
        return new BoundaryDataset(posValues, negValues);                
    }
    
    public BoundaryDataset skipsampleDataset(int offset, int skip) { 
        LinkedList<Double> posValues = new LinkedList<Double>();
        LinkedList<Double> negValues = new LinkedList<Double>();

        int negOffset = 0;
        for(double val : negative.keySet()) { 
            for(int i = 0; i < negative.get(val); i++) { 
                if((negOffset-offset) % skip == 0) { negValues.addLast(val); }
                negOffset++;
            }
        }
        
        int posOffset = 0;
        for(double val : positive.keySet()) { 
            for(int i = 0; i < positive.get(val); i++) { 
                if((posOffset-offset) % skip == 0) { posValues.addLast(val); }
                posOffset++;
            }
        }
        

        return new BoundaryDataset(posValues, negValues);                
    }
    
    public BoundaryDataset subsampleDataset(int np, int nn, boolean replacement) { 
        LinkedList<Double> posValues = new LinkedList<Double>();
        LinkedList<Double> negValues = new LinkedList<Double>();
        
        if(replacement) { 
            for(int i = 0; i < np; i++) { 
                posValues.addLast(chooseRandomPositiveValue());
            }
            for(int i = 0; i < nn; i++) { 
                negValues.addLast(chooseRandomNegativeValue());
            }
        } else { 
            double[] pa = chooseRandomRemovalPositive(np);
            double[] na = chooseRandomRemovalNegative(nn);
            for(int i = 0; i < pa.length; i++) { posValues.addLast(pa[i]); }
            for(int i = 0; i < na.length; i++) { negValues.addLast(na[i]); }
        }
        
        return new BoundaryDataset(posValues, negValues);
    }
    
    public boolean chooseRandomLabel() { 
        double posProb = (double)numPos / (double)(numPos + numNeg);
        double val = rand.nextDouble();
        return val <= posProb;
    }
    
    public double[] chooseRandomRemovalPositive(int num) { 
        if(num <= 0 || num > numPos) { throw new IllegalArgumentException(); }
        TreeMap<Double,Integer> copy = new TreeMap<Double,Integer>(positive);
        int copySize = numPos;
        double[] array = new double[num];
        
        for(int i = 0; i < array.length; i++) { 
            array[i] = chooseRandomMapValue(copySize, copy);
            if(copy.get(array[i]) == 0) { 
                copy.remove(array[i]);
            } else { 
                copy.put(array[i], copy.get(array[i]) - 1);
            }
            copySize -= 1;
        }
        
        return array;
    }

    public double[] chooseRandomRemovalNegative(int num) { 
        if(num <= 0 || num > numPos) { throw new IllegalArgumentException(); }
        TreeMap<Double,Integer> copy = new TreeMap<Double,Integer>(negative);
        int copySize = numNeg;
        double[] array = new double[num];
        
        for(int i = 0; i < array.length; i++) { 
            array[i] = chooseRandomMapValue(copySize, copy);
            if(copy.get(array[i]) == 0) { 
                copy.remove(array[i]);
            } else { 
                copy.put(array[i], copy.get(array[i]) - 1);
            }
            copySize -= 1;
        }
        
        return array;
    }
    

    public double chooseRandomMapValue(int size, SortedMap<Double,Integer> map) { 
        int valueIndex = rand.nextInt(size);
        for(double val : map.keySet()) { 
            valueIndex -= map.get(val);
            if(valueIndex <= 0) { 
                return val; 
            }
        }
        return map.lastKey();        
    }
    
    public double chooseRandomPositiveValue() { return chooseRandomMapValue(numPos, positive); }
    public double chooseRandomNegativeValue() { return chooseRandomMapValue(numNeg, negative); }

    public String createStringRep(int polarity) { return createStringRep(polarity, false); }
    public String createStringRep(int polarity, boolean bound) {  
        TreeSet<LabeledPoint> lp = new TreeSet<LabeledPoint>();
        int plabel = polarity == 1 ? 1 : 0;
        int nlabel = 1 - plabel;
        
        for(double v : positive.keySet()) { lp.add(new LabeledPoint(v, plabel)); }
        for(double v : negative.keySet()) { lp.add(new LabeledPoint(v, nlabel)); }
        StringBuilder sb = new StringBuilder();
    
        for(LabeledPoint lpt : lp) { 
            sb.append(String.valueOf(lpt.getLabel()));
        }

        if(bound && boundary != null) { 
            LabeledPoint bp = new LabeledPoint(boundary, 0);
            SortedSet<LabeledPoint> lower = lp.headSet(bp);
            int lowerSize = lower.size();
            String str = " << ";
            if(boundDir == 1) { str = " >> "; }
            sb.insert(lowerSize, str);
        }

        return sb.toString();        
    }
    
	private void assembleStringRep() { 
		TreeSet<LabeledPoint> lp = new TreeSet<LabeledPoint>();
		for(double v : positive.keySet()) { lp.add(new LabeledPoint(v, 1)); }
		for(double v : negative.keySet()) { lp.add(new LabeledPoint(v, 0)); }
		StringBuilder sb = new StringBuilder();
	
		for(LabeledPoint lpt : lp) { 
			sb.append(String.valueOf(lpt.getLabel()));
		}

		if(boundary != null) { 
			LabeledPoint bp = new LabeledPoint(boundary, 0);
			SortedSet<LabeledPoint> lower = lp.headSet(bp);
			int lowerSize = lower.size();
			String str = " << ";
			if(boundDir == 1) { str = " >> "; }
			sb.insert(lowerSize, str);
		}

		stringRep = sb.toString();
	}

	private static class LabeledPoint implements Comparable<LabeledPoint> { 
		private double value;
		private int label;

		public LabeledPoint(double d, int l) { 
			value = d;
			label = l;
		}

		public double getValue() { return value; }
		public int getLabel() { return label; }

		public boolean equals(Object o) { 
			if(!(o instanceof LabeledPoint)) { return false; }
			LabeledPoint lp = (LabeledPoint)o;
			return lp.value == value && lp.label == label;
		}

		public int hashCode() { 
			int code = 17;
			long bits = Double.doubleToLongBits(value);
			code += (int)(bits >> 32); code *= 37;
			code += label; code *= 37;
			return code;
		}

		public int compareTo(LabeledPoint lp) { 
			if(value < lp.value) { return -1; }
			if(value > lp.value) { return 1; }
			if(label < lp.label) { return -1; }
			if(label > lp.label) { return 1; }
			return 0;
		}
	}
    
    public double[] createSortedPointArray() { 
        double[] array = new double[positive.size() + negative.size()];
        int i = 0;
        for(double pv : positive.keySet()) { array[i++] = pv; }
        for(double nv : negative.keySet()) { array[i++] = nv; }
        Arrays.sort(array);
        return array;
    }
    
    public void findBoundary() { 
        if(boundary == null) {
            double[] array = createSortedPointArray();
            double currentBound = array[0] - 1.0;
            int currentDir = -1;
            int currentError = numPos;

            int firstErrors = countErrors(currentBound, 1);
			if(firstErrors < currentError) { 
				currentDir = 1; 
				currentError = firstErrors;
			}
            
            for(int i = 0; i < array.length-1; i++) { 
                double nextBound = (array[i] + array[i+1]) / 2.0;
                int upErrors = countErrors(nextBound, 1); 
                int downErrors = countErrors(nextBound, -1);
                
                if(upErrors < currentError) { 
                    currentError = upErrors;
                    currentBound = nextBound;
                    currentDir = 1;
                }

                if(downErrors < currentError) { 
                    currentError = downErrors;
                    currentBound = nextBound;
                    currentDir = -1;
                }
            }
            
            double lastBound = array[array.length-1] + 1.0;
            int lastErrors = countErrors(lastBound, -1);
            if(lastErrors < currentError) { 
                currentBound = lastBound;
                currentDir = -1;
                currentError = lastErrors;
            }
            
            boundary = currentBound;
            boundDir = currentDir;
            error = currentError;
           
		   	assembleStringRep();
        }
    }
    
    public int countMap(Map<Double,Integer> m) {
        int count = 0;
        for(double v : m.keySet()) { 
            count += m.get(v);
        }
        return count;
    }
    
    public int countErrors(double bound, int dir) { 
        if(dir == -1) { 
            SortedMap<Double,Integer> lesserNegative = negative.headMap(bound);
            SortedMap<Double,Integer> greaterPositive = positive.tailMap(bound);
            int negErrors = countMap(lesserNegative);
            int posErrors = countMap(greaterPositive);
            return negErrors + posErrors;
        } else { 
            SortedMap<Double,Integer> greaterNegative = negative.tailMap(bound);
            SortedMap<Double,Integer> lesserPositive = positive.headMap(bound);
            int negErrors = countMap(greaterNegative);
            int posErrors = countMap(lesserPositive);
            return negErrors + posErrors;
        }
    }
}
