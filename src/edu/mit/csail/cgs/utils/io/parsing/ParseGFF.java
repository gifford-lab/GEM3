package edu.mit.csail.cgs.utils.io.parsing;

import java.io.BufferedReader;
import java.io.*;
import java.util.*;
import java.net.URLEncoder;

public class ParseGFF implements Iterator { 
    
    public static void main(String[] args) { 
        try { 
            ParseGFF parser = new ParseGFF(new File(args[0]));
            parser.printLines(System.out);
        } catch(IOException ie) { 
            ie.printStackTrace();
        }
    }
    
    private BufferedReader br;    
    private int lineNum;
    private String line, filename;
    private boolean dirty;

    public ParseGFF(File f) throws IOException { 
        br = new BufferedReader(new FileReader(f));
        lineNum = 0;
        dirty = true;
        filename = f.getName();
    }

    public boolean hasNext() {
        if (dirty) {
            try {
                lineNum++;
                line = br.readLine();
            } catch (IOException ex) {
                throw new RuntimeException("Parsing Error, File \"" + filename + "\", line " + lineNum);
            }
            dirty = false;
        }
        if (line == null) {
            try {
                br.close();
            } catch (IOException ex) {
                throw new RuntimeException("Can't close " + filename);
            }
            return false;
        } else {
            if (line.startsWith("#")) {
                dirty = true;
                return hasNext();
            } else {
                return true;
            }
        }
    }
    public Line next() throws NoSuchElementException {
        if (line == null) {
            throw new NoSuchElementException("No more lines to parse");
        }
        dirty = true;
        line = line.trim();
        if(!line.startsWith("#")) { 
            try { 
                Line gffLine = new Line(line);
                return gffLine;
            } catch(NoSuchElementException e) { 
                throw new RuntimeException("Parsing Error, File \"" + filename + "\", line " + lineNum);
            }
        } else {
            return next();
        }
    }
    public void remove() throws UnsupportedOperationException {
        throw new UnsupportedOperationException("Can't remove lines from GFF file");
    }
    
    public void printLines() { printLines(System.out); }
    public void printLines(PrintStream ps) { 
        while(hasNext()) { 
            next().printLine(ps);
        }
    }
    
    public static class Line { 
        private String fSeq, fSource, fFeature;
        private int fStart, fEnd;
        private double fScore;
        private String fStrand, fFrame;
        private String fAttribString;
        private Set<String> fSimpleAttribs;
        private Map<String,String> fComplexAttribs;
        
        public Line(String line) { 
            StringTokenizer st = new StringTokenizer(line, "\t");
            fSeq = st.nextToken();
            fSource = st.nextToken();
            fFeature = st.nextToken();
            fStart = Integer.parseInt(st.nextToken());
            fEnd = Integer.parseInt(st.nextToken());
            try { 
                fScore = Double.parseDouble(st.nextToken());
            } catch(NumberFormatException nfe) { 
                fScore = 1.0;
            }
            
            fStrand = st.nextToken();
            fFrame = st.nextToken();
            String attribs = st.nextToken();
            fAttribString = attribs;
            st = new StringTokenizer(attribs, ";");
            fSimpleAttribs = new HashSet<String>();
            fComplexAttribs = new HashMap<String,String>();
            
            while(st.hasMoreTokens()) { 
                String aTok = st.nextToken();
                if(aTok.indexOf("=") != -1) { 
                    String key = aTok.substring(0, aTok.indexOf("="));
                    String value = aTok.substring(aTok.indexOf("=") + 1, aTok.length());
                    fComplexAttribs.put(key, value);
				} else if(aTok.startsWith("Site ")) { 
					String value = aTok.substring(5, aTok.length());
					fComplexAttribs.put("Site", value);
                } else {
                    fSimpleAttribs.add(aTok);
                }
            }
        }
        
        public void printLine() { printLine(System.out); }
        public void printLine(PrintStream ps) { 
            ps.println("[" + fFeature + "] (" + fSeq + "," + fSource + ")");
            ps.println(String.valueOf(fScore) + ": " + fStart + " - " + fEnd);
            Iterator itr = fSimpleAttribs.iterator();
            while(itr.hasNext()) { 
                ps.println("\t" + itr.next());
            }
            
            itr = fComplexAttribs.entrySet().iterator();
            while(itr.hasNext()) { 
                Map.Entry e = (Map.Entry)itr.next();
                ps.println("\t" + e.getKey() + " ==> " + e.getValue());
            }
        }
        
        public boolean hasComplexAttrib(String key) { return fComplexAttribs.containsKey(key); }
        public String getComplexAttrib(String key) { 
            return (String)fComplexAttribs.get(key);
        }
        
        public Set<String> getPrefixedSimpleAttribs(String pref) { 
            Set<String> hashSet = new HashSet<String>();
            Iterator itr = fSimpleAttribs.iterator();
            while(itr.hasNext()) { 
                String attrib = (String)itr.next();
                if(attrib.startsWith(pref)) { hashSet.add(attrib); }
            }
            
            return hashSet;
        }
        
        public int hashCode() { 
            int code = 17;
            code += fSeq.hashCode(); code *= 37;
            code += fSource.hashCode(); code *= 37;
            code += fFeature.hashCode(); code *= 37;
            code += fStart; code *= 37;
            code += fEnd; code *= 37;
            code += fFrame.hashCode(); code *= 37;
            code += fStrand.hashCode(); code *= 37;
            long bits = Double.doubleToLongBits(fScore);
            code += (int)(bits >> 32); code *= 37;
            Iterator itr = fSimpleAttribs.iterator();
            while(itr.hasNext()) { 
                code += itr.next().hashCode(); code *= 37;
            }
            
            itr = fComplexAttribs.keySet().iterator();
            while(itr.hasNext()) { 
                String key = (String)itr.next();
                code += key.hashCode();
                code += fComplexAttribs.get(key).hashCode();
                code *= 37;
            }
            
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Line)) { return false; }
            Line l = (Line)o;
            if(!fSeq.equals(l.fSeq)) { return false; }
            if(!fSource.equals(l.fSource)) { return false; }
            if(!fFeature.equals(l.fFeature)) { return false; }
            if(fStart != l.fStart || fEnd != l.fEnd) { return false; }
            if(!fStrand.equals(l.fStrand)) { return false; }
            if(!fFrame.equals(l.fFrame)) { return false; }
            if(fScore != l.fScore) { return false; }
            if(fSimpleAttribs.size() != l.fSimpleAttribs.size()) { return false; }
            if(fComplexAttribs.keySet().size() != l.fComplexAttribs.keySet().size()) { return false; }
            Iterator itr = fSimpleAttribs.iterator();
            while(itr.hasNext()) { 
                String key = (String)itr.next();
                if(!l.fSimpleAttribs.contains(key)) { return false; }
            }
            
            itr = fComplexAttribs.keySet().iterator();
            while(itr.hasNext()) { 
                String key = (String)itr.next();
                if(!l.fComplexAttribs.containsKey(key)) { return false; }
                String v1 = (String)fComplexAttribs.get(key);
                String v2 = (String)l.fComplexAttribs.get(key);
                if(!v1.equals(v2)) { return false; }
            }
            return true;
        }
        
        public String getSeq() { return fSeq; }
        public String getSource() { return fSource; }
        public String getFeature() { return fFeature; }
        public int getStart() { return fStart; }
        public int getEnd() { return fEnd; }
        public double getScore() { return fScore; }
        public String getStrand() { return fStrand; }
        public String getFrame() { return fFrame; }

        public String getAttribString() { 
			StringBuilder sb = new StringBuilder();
			for(String key : fComplexAttribs.keySet()) { 
				String value = fComplexAttribs.get(key);
				if(sb.length() > 0) { sb.append(";"); }
				try { 
					sb.append(URLEncoder.encode(key, "UTF-8"));
				} catch(UnsupportedEncodingException e) { 
					sb.append(key);
				}
				sb.append("=");
				try { 
					sb.append(URLEncoder.encode(value, "UTF-8"));
				} catch(UnsupportedEncodingException e) { 
					sb.append(value);
				}
			}
			return sb.toString();
		}
        
        public boolean hasAttrib(String v) { 
            return fSimpleAttribs.contains(v) || fComplexAttribs.containsKey(v); 
        }
        
        public String getAttrib(String v) { 	
            if(!fComplexAttribs.containsKey(v)) { return null; }
            return (String)fComplexAttribs.get(v);
        }
    }
}
