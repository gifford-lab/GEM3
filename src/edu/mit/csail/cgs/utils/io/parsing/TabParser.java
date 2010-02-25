/*
 * Created on May 25, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.ArgParser;

/**
 * @author tdanford
 * 
 * This is a parser for "expression" data that comes as a tab-separated,
 * tabular-formatted file (rows are probes, columns are experiments).  In
 * particular, it was created to handle the parsing of the Mnaimneh data for
 * input into the database. 
 */
public class TabParser implements ParseExpr {

	public static void main(String[] args) { 
		ArgParser ap = new ArgParser(args);
		File f = new File(ap.getKeyValue("expr"));
		TabParser tp = new TabParser(f);
		interactive(tp);
	}

	public static void interactive(ParseExpr pe) { 
		try { 
			System.out.print(">"); System.out.flush();
			String line;
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			while((line = br.readLine()) != null) { 
				line = line.trim();
				if(line.equals("expts")) { 
					for(int i = 0; i < pe.numExpts(); i++) { 
						System.out.println(i + ": " + pe.getExptName(i));
					}
				} else { 
					StringTokenizer st = new StringTokenizer(line);
					String orf = st.nextToken();
					int expt = Integer.parseInt(st.nextToken());
					if(pe.exptHasORF(expt, orf)) { 
						double value = pe.getValue(expt, orf);
						System.out.println(orf + "," + expt + " --> " + value);
					} else { 
						System.out.println(orf + " --> no value in " + expt);
					}
				}
				System.out.print(">"); System.out.flush();
			}
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
			throw new RuntimeException(ie);
		}
	}
    
    private Map<String,Integer> probeIndices, exptIndices;
    private String[] probeNames, exptNames;
    private Double[][] values;
    private int[] probeCounts;
    
    public TabParser(File f) { 
        parse(f);
    }
    
    private void parse(File f) {
        LinkedList<String> probeList = new LinkedList<String>();
        LinkedList<String> exptList = new LinkedList<String>();
        LinkedList<Double[]> valueList = new LinkedList<Double[]>();
        String line;
        
        try { 
            BufferedReader br = new BufferedReader(new FileReader(f));
            line = br.readLine();
            StringTokenizer st = new StringTokenizer(line, "\t");
            while(st.hasMoreTokens()) { exptList.addLast(st.nextToken()); }
            System.out.println("# Expts: " + exptList.size());
            
			int errorCount = 0, successCount = 0;
            while((line = br.readLine()) != null) { 
                st = new StringTokenizer(line, "\t");
				int tokCount = st.countTokens();
				String probe = st.nextToken();
                if(tokCount != (1 + exptList.size())) { 
                    System.out.println("ERROR: Probe Line [" + probe + "] has too few (" + tokCount + " != " + (1+probeList.size()) + ") tokens.");
					errorCount += 1;
                } else { 
                    Double[] array = new Double[exptList.size()];
                    for(int i = 0; i < exptList.size(); i++) { 
                        try { 
                            array[i] = new Double(st.nextToken());
                        } catch(NumberFormatException nfe) { 
                            array[i] = null;
                        }
                    }
                    valueList.addLast(array);
                    probeList.addLast(probe);
					successCount += 1;
                }
            }
            
            System.err.println("# Probes: " + probeList.size());
			if(errorCount > 0) { 
				System.err.println("\tFailure Count: " + errorCount);
			}
            br.close();
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
        probeNames = probeList.toArray(new String[probeList.size()]);
		//System.out.println("# Probe Names: " + probeNames.length);

        exptNames = exptList.toArray(new String[exptList.size()]);
		//System.out.println("# Expt Names: " + exptNames.length);

        values = valueList.toArray(new Double[valueList.size()][]);
		//System.out.println("# Value Rows: " + values.length);

        probeCounts = new int[exptNames.length];
        fillProbeCounts();

        probeIndices = new HashMap<String,Integer>();
        exptIndices = new HashMap<String,Integer>();
        for(int i = 0; i < probeNames.length; i++) { probeIndices.put(probeNames[i], i); }
        for(int i = 0; i < exptNames.length; i++) { exptIndices.put(exptNames[i], i); }
    }
    
    private void fillProbeCounts() { 
        for(int e = 0; e < exptNames.length; e++) {
            int count = 0;
            for(int p = 0; p < probeNames.length; p++) { 
                if(values[p][e] != null) { 
                    count += 1;
                }
            }
            probeCounts[e] = count;
        }
    }
    
    public String getExptName(int i) { 
        return exptNames[i];
    }
    
    public Map<String,Double> getWholeExpt(int i) { 
        Map<String,Double> vals = new HashMap<String,Double>();
        for(int j = 0; j < probeNames.length; j++) { 
            if(values[j][i] != null) { 
                vals.put(probeNames[j], values[j][i]);
            }
        }
        
        return vals;
    }
    
    public boolean exptHasORF(int i, String orf) { 
        return probeIndices.get(orf) != null && values[probeIndices.get(orf)][i] != null;
    }
    
    public double getValue(int i, String orf) { 
        return values[probeIndices.get(orf)][i].doubleValue();
    }

    public Double[] getValues(int start, int stop, String orf) { 
        Double[] array = new Double[stop-start+1];
        for(int i = 0; i < start-stop+1; i++) { 
            if(exptHasORF(start+i, orf)) { 
                array[i] = new Double(getValue(start+i, orf));
            } else { 
                array[i] = null;
            }
        }
        return array;
    }
    
    public int numExpts() { 
        return exptNames.length;
    }
    
    public int numORFs(int i) {
        return probeCounts[i];
    }
    
    public double getMax(int i) { 
        double extremum = -Double.MAX_VALUE;
        for(int p = 0; p < probeNames.length; p++) { 
            if(values[p][i] != null && values[p][i].doubleValue() > extremum) { 
                extremum = values[p][i].doubleValue();
            }
        }
        return extremum;
    }
    
    public double getMin(int i) {     
        double extremum = Double.MAX_VALUE;
        for(int p = 0; p < probeNames.length; p++) { 
            if(values[p][i] != null && values[p][i].doubleValue() < extremum) { 
                extremum = values[p][i].doubleValue();
            }
        }
        return extremum;        
    }
}
