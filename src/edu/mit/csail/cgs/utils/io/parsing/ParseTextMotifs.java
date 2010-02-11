package edu.mit.csail.cgs.utils.io.parsing;

import java.lang.*;
import java.io.*;
import java.util.*;

public class ParseTextMotifs { 

	public static void main(String[] args) { 
		try { 
			ParseTextMotifs ptm = new ParseTextMotifs(new File(args[0]));
            for(ParsedMotif m : ptm.motifList) { 
                System.out.println(m);
            }
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}
	
	private LinkedList<ParsedMotif> motifList;

	public ParseTextMotifs(File f) throws IOException { 
	    motifList = new LinkedList<ParsedMotif>();
	    BufferedReader br = new BufferedReader(new FileReader(f));
	    try { 
	        while(true) { 
	            ParsedMotif m = new ParsedMotif(br);
	            motifList.addLast(m);
	        }
	    } catch(EOFException e) { 
	        for(ParsedMotif m : motifList) { System.out.println(m.toString()); }
	    }
	    br.close();
	}
    
    public List<ParsedMotif> findMotifs(String name) { 
        LinkedList<ParsedMotif> motifs = new LinkedList<ParsedMotif>();
        for(ParsedMotif m : motifList) { 
            if(m.getSource().indexOf(name) != -1) { 
                motifs.addLast(m);
            }
        }
        return motifs;
    }
	
	public static class ParsedMotif { 
	    private String source;
	    private double max_score;
	    private Map<String,double[]> probs;
	    private int cols;
	    
	    public String toString() { 
	        return "[" + source + "] : score (" + max_score + "/" + calcCutoff(0.6, true) + ") --> " + getMaxMotifString();
	    }
	    
	    public String getSource() { return source; }
        public double getMaxScore() { return max_score; }
        public int numCols() { return cols; }
        
        public double[] scoreDNASequence(String s) { 
            double[] array = new double[s.length() - cols + 1];
            for(int start = 0; start < s.length()-cols; start++) { 
                array[start] = scoreDNAString(s, start);
            }
            return array;
        }
       
        public boolean containsMatch(String s) { 
            double cutoff = calcCutoff(0.6, true);
            return containsMatch(s, cutoff);
        }
        
        public boolean containsMatch(String s, double score) { 
            for(int start = 0; start < s.length()-cols; start++) { 
                double v = scoreDNAString(s, start);
                if(v >= score) { return true; }
            }
            return false;
        }
        
        public double scoreDNAString(String s, int offset) {
            double sum = 0.0;
            for(int i = 0; i < cols; i++) { 
				try { 
					sum += probs.get(String.valueOf(Character.toUpperCase(s.charAt(offset + i))))[i];
				} catch(NullPointerException npe) { 
					System.err.println("NULL Value: " + s.charAt(offset + i));
					System.err.print("Values (");
					for(String k : probs.keySet()) { System.err.print(k + " "); }
					System.err.print(")\n");
				}
            }
            return sum;
        }
        
        public double calcCutoff(double factor, boolean sum) { 
			return factor * calcMaxScore(sum);
        }
        
        public double calcMaxScore(boolean sum) { 
            double base;
            if(sum) { base = 0.0; } else { base = 1.0; }
            for(int i = 0; i < cols; i++) {
                double max = 0.0;
                String max_key = null;
                for(String key : probs.keySet()) {
                    if(max_key == null || probs.get(key)[i] > max) { 
                        max_key = key; max = probs.get(key)[i];
                    }
                }
                
                if(sum) { base += max; } else { base *= max; }
            }
            return base;
        }

        public double calcMinScore(boolean sum) { 
            double base;
            if(sum) { base = 0.0; } else { base = 1.0; }
            for(int i = 0; i < cols; i++) {
                double max = 0.0;
                String max_key = null;
                for(String key : probs.keySet()) {
                    if(max_key == null || probs.get(key)[i] < max) { 
                        max_key = key; max = probs.get(key)[i];
                    }
                }
                
                if(sum) { base += max; } else { base *= max; }
            }
            return base;
        }

	    public String getMaxMotifString() { 
	        StringBuilder sb = new StringBuilder();
	        for(int i = 0; i < cols; i++) { 
	            String maxString = null;
	            double maxScore = 0.0;
	            for(String key : probs.keySet()) { 
	                double s = probs.get(key)[i];
	                if(maxString == null || s > maxScore) { 
	                    maxString = key;
	                    maxScore = s;
	                }
	            }
	            sb.append(maxString);
	        }
	        return sb.toString();
	    }
        
	    public ParsedMotif(BufferedReader br) throws IOException { 
	        String line = br.readLine();
	        if(line==null) { throw new EOFException(); }
	        int firstColon = line.indexOf(":");
	        source = line.substring(firstColon+1, line.length());
	        line = br.readLine();
	        firstColon = line.indexOf(":");
	        max_score = Double.parseDouble(line.substring(firstColon+1, line.length()).trim());
	        line = br.readLine();
	        if(line.startsWith("#")) { line = line.substring(1, line.length()); }
	        StringTokenizer st = new StringTokenizer(line);
	        cols = st.countTokens();
	        probs = new HashMap<String,double[]>();
	        
	        line = br.readLine();
	        if(line.startsWith("#")) { line = line.substring(1, line.length()); }
	        st = new StringTokenizer(line);
	        double[] array = new double[st.countTokens()-1];
	        String name = st.nextToken();
	        for(int i = 1; i < cols; i++) { array[i-1] = Double.parseDouble(st.nextToken()); }
	        probs.put(name, array);
	        
	        line = br.readLine();
	        if(line.startsWith("#")) { line = line.substring(1, line.length()); }
	        st = new StringTokenizer(line);
	        array = new double[st.countTokens()-1];
	        name = st.nextToken();
	        for(int i = 1; i < cols; i++) { array[i-1] = Double.parseDouble(st.nextToken()); }
	        probs.put(name, array);
	        
	        line = br.readLine();
	        if(line.startsWith("#")) { line = line.substring(1, line.length()); }
	        st = new StringTokenizer(line);
	        array = new double[st.countTokens()-1];
	        name = st.nextToken();
	        for(int i = 1; i < cols; i++) { array[i-1] = Double.parseDouble(st.nextToken()); }
	        probs.put(name, array);
	        
	        line = br.readLine();
	        if(line.startsWith("#")) { line = line.substring(1, line.length()); }
	        st = new StringTokenizer(line);
	        array = new double[st.countTokens()-1];
	        name = st.nextToken();
	        for(int i = 1; i < cols; i++) { array[i-1] = Double.parseDouble(st.nextToken()); }
	        probs.put(name, array);
	        
	    }
	}

}
