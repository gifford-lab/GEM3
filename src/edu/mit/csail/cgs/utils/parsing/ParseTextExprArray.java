/*
 * Created on Mar 9, 2005
 */
package edu.mit.csail.cgs.utils.parsing;

import java.util.*;
import java.io.*;
import java.lang.*;

/**
 * @author tdanford
 */
public class ParseTextExprArray implements ParseExpr {

	public static void main(String[] args) { 
		try { 
			ParseTextExprArray array = new ParseTextExprArray(new File(args[0]));
			System.out.print(">"); System.out.flush();
			String line;
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			while((line = br.readLine()) != null) { 
				StringTokenizer st = new StringTokenizer(line);
				String orf = st.nextToken();
				int expt = Integer.parseInt(st.nextToken());
				if(!array.exptHasORF(expt, orf)) { 
					System.out.println("ORF " + orf + " has no value in " + expt);
					System.out.println("\tTotal Values: " + array.countExpts(orf));
					System.out.println("\tExpt Values: " + array.values.get(expt).size());
					for(String korf : array.values.get(expt).keySet()) { 
						System.out.print(korf + " ");
					}
					System.out.print("\n");
				} else { 
					System.out.println(orf + " ==> " + array.getValue(expt, orf));
				}
				System.out.print(">"); System.out.flush();
			}
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}

    private Vector<String> exptNames;
    private Vector<String> orfNames;
    private Vector<Map<String,Double>> values;
    private Map<Integer,Double> maxMap, minMap;
    
    public ParseTextExprArray(File f) throws IOException {
        parse(f, null);
    }
    
    public ParseTextExprArray(File f, Set<Integer> exptInds) throws IOException { 
        parse(f, exptInds);
    }
    
    private void parse(File f, Set<Integer> exptInds) throws IOException { 
        int exptNum = 0;
        BufferedReader br = new BufferedReader(new FileReader(f));
        
        String line = br.readLine();
		StringTokenizer st = new StringTokenizer(line);
		int orfCount = Integer.parseInt(st.nextToken());
		int exptCount = Integer.parseInt(st.nextToken());

        orfNames = readORFLine(br);
        exptNames = new Vector<String>();
        values = new Vector<Map<String,Double>>();
        maxMap = new HashMap<Integer,Double>();
        minMap = new HashMap<Integer,Double>();

		for(int i = 0; i < exptCount; i++) { 
			exptNames.add(br.readLine().trim());
			values.add(new HashMap<String,Double>());
			minMap.put(i, Double.MAX_VALUE); maxMap.put(i, -Double.MAX_VALUE);
		}

		System.out.print("[" + orfCount + "] "); System.out.flush();
		for(int i = 0; i < orfCount; i++) { 
			parseValueLine(orfNames.get(i), br, exptInds);
			if(i % 100 == 0) { System.out.print("."); System.out.flush(); }
			if(i % 1000 == 0) { System.out.print("(" + i + ")"); System.out.flush(); }
		}
		System.out.print("\n");
        
        br.close();
    }
    
    private Vector<String> readORFLine(BufferedReader br) throws IOException { 
        String line = br.readLine();
        StringTokenizer st = new StringTokenizer(line);
        Vector<String> orfs = new Vector<String>();

		int count = Integer.parseInt(st.nextToken());
        
		int added = 0;
        while(st.hasMoreTokens()) { 
            String orf = st.nextToken();
            orfs.add(orf);
			added += 1;
        }
		System.out.println(added + "/" + count);
        
        return orfs;
    }
    
    private void parseValueLine(String orf, BufferedReader br, Set<Integer> exptInds) throws IOException {
        Map<String,Double> valueMap = new TreeMap<String,Double>();
        
        String line = br.readLine();
		String[] vArray = line.split(" ");
		for(int i = 0; i < vArray.length; i++) { 
			if(exptInds == null || exptInds.contains(i)) { 
				try {
					String v = vArray[i];
					double d = Double.parseDouble(v);
					values.get(i).put(orf, d);
					if(d < minMap.get(i)) { minMap.put(i, d); }
					if(d > maxMap.get(i)) { maxMap.put(i, d); }
				} catch(NumberFormatException nfe) { 
					// do nothing.
				}
			}
        }
    }
    
    public String getExptName(int i) { return exptNames.get(i); }
    public String getORFName(int i) { return orfNames.get(i); }

	public int countExpts(String orf) { 
		int count = 0;
		for(int i = 0; i < values.size(); i++) { 
			if(values.get(i).containsKey(orf)) { count += 1; }
		}
		return count;
	}

	public Double[] getValues(int start, int end, String orf) { 
		Double[] array = new Double[end-start+1];
		for(int i = 0; i < array.length; i++) { 
			if(exptHasORF(start+i, orf)) { 
				array[i] = getValue(start+i, orf);
			} else { 
				array[i] = null;
			}
		}
		return array;
	}
    
    public Map<String,Double> getWholeExpt(int i) { return values.get(i); }
    public boolean exptHasORF(int i, String orf) { return values.get(i).containsKey(orf); }
    public double getValue(int i, String orf) { return values.get(i).get(orf); }
    public int numExpts() { return exptNames.size(); }
    public int numORFs(int i) { return values.get(i).size(); }
    public double getMax(int i) { 
		return maxMap.get(i); 
	}
    public double getMin(int i) { return minMap.get(i); }
}
