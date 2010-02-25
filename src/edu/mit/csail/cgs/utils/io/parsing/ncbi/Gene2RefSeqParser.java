/*
 * Created on Aug 15, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.ncbi;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 * 
 * Built to parse the "gene2refseq" file that comes from NCBI.
 * 
 * For instance, see this URL: "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"
 */
public class Gene2RefSeqParser {
    
    public static void main(String[] args) { 
        try {
            //Gene2RefSeqParser parser = args.length > 0 ? new Gene2RefSeqParser(new File(args[0])) : new Gene2RefSeqParser();
            DynamicLoader loader = new DynamicLoader(new File(args[0]));
            int count = 0;
            while(loader.hasNext()) { 
                loader.next();
                count += 1;
            }
            
            System.out.println("# Entries: " + count);
            
                
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public static Map<String,Set<String>> buildEntrez2RefSeqMap(File f) throws IOException { 
    	Map<String,Set<String>> map = new HashMap<String,Set<String>>();
    	DynamicLoader loader = new DynamicLoader(f);
    	
    	while(loader.hasNext()) { 
    		Entry e = loader.next();
    		String entrez = e.geneID;
    		String refseq = e.rnaAcc;
    		if(!refseq.equals("-")) { 
    			if(!map.containsKey(entrez)) { 
    				map.put(entrez, new HashSet<String>());
    			}
    			map.get(entrez).add(refseq);
    		}
    	}
    	
    	if(!loader.isClosed()) { loader.close(); } 
    	return map;
    }
    
    public static Map<String,Set<String>> buildEntrez2RefSeqMap(File f, Set<String> species) throws IOException { 
    	Map<String,Set<String>> map = new HashMap<String,Set<String>>();
    	DynamicLoader loader = new DynamicLoader(f);
    	
    	while(loader.hasNext()) { 
    		Entry e = loader.next();
    		String entrez = e.geneID;
    		String refseq = e.rnaAcc;
    		if(species.contains(e.taxID) && !refseq.equals("-")) { 
    			if(!map.containsKey(entrez)) { 
    				map.put(entrez, new HashSet<String>());
    			}
    			map.get(entrez).add(refseq);
    		}
    	}
    	
    	if(!loader.isClosed()) { loader.close(); } 
    	return map;
    }
    
    public static class DynamicLoader implements Iterator<Entry> {
        
        private Entry nextEntry;
        private BufferedReader br;
        
        public DynamicLoader(File f) throws IOException { 
            br = new BufferedReader(new FileReader(f));
            loadNext();
        }

        public boolean hasNext() {
            return nextEntry != null;
        }

        public Entry next() {
            Entry e = nextEntry;
            loadNext();
            return e;
        }
        
        private void loadNext() { 
            if(br != null) { 
                nextEntry = null;
                try {
                    String line = br.readLine();
                    if(line == null) { 
						//System.out.println("Null line.");
                        close(); 
                    } else { 
						//System.out.println("Entry parsed.");
                        nextEntry = new Entry(line);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                    close();
                }
            } else { 
				//System.out.println("Null BR.");
			}
        }
        
        public boolean isClosed() { return br == null; }
        
        public void close() { 
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            br = null;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }         
    }
    
    private Vector<Entry> entries;
    
    public Gene2RefSeqParser(File f) throws IOException { parse(f); }
    public Gene2RefSeqParser() throws IOException { parse(new File("gene2refseq")); }
    
    private void parse(File f) throws IOException { 
        entries = new Vector<Entry>();
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                Entry e = new Entry(line);
                entries.add(e);
            }
        }
        br.close();
        System.out.println("Loaded " + entries.size() + " entries from \"" + f.getName() + "\"");
    }
    
    public Set<Entry> selectEntries(String taxID, String geneID, String rnaAcc) { 
        HashSet<Entry> selected = new HashSet<Entry>();
        
        for(Entry e : entries) { 
            boolean isSel = true;
            if(isSel && taxID != null && !e.taxID.equals(taxID)) { isSel = false; }
            if(isSel && geneID != null && !e.geneID.equals(geneID)) { isSel = false; }
            if(isSel && rnaAcc != null && !e.rnaAcc.equals(rnaAcc)) { isSel = false; }
            
            if(isSel) { selected.add(e); }
        }
        
        return selected;
    }
    
    public static class Entry { 
        public String taxID, geneID;
        public String status;
        public String rnaAcc, rnaVersion, rnaGI;
        public String protAcc, protVersion, protGI;
        public String genAcc, genVersion, genGI;
        public int start, end, orientation;
        public String assembly;
        
        public Entry(String line) { 
            String[] array = line.split("\t");
            if(array.length != 13) { 
                throw new IllegalArgumentException("Need 13 tokens in line \"" + line + "\""); 
            }
            
            taxID = array[0]; 
            geneID = array[1];
            status = array[2];
            
			if(!array[3].equals("-")) { 
				String[] ra = array[3].split("\\.");
				rnaAcc = ra[0]; rnaVersion = ra.length > 1 ? ra[1] : "-";
			} else { 
				rnaAcc = rnaVersion = "-";
			}
            rnaGI = array[4];

			if(!array[5].equals("-")) { 
				String[] pa = array[5].split("\\.");
				protAcc = pa[0]; protVersion = pa.length > 1 ? pa[1] : "-";
			} else { 
				protAcc = protVersion = "-";
			}
            protGI = array[6];

			if(!array[7].equals("-")) { 
				String[] ca = array[7].split("\\.");
				genAcc = ca[0]; genVersion = ca.length > 1 ? ca[1] : "-";
			} else { 
				genAcc = genVersion = "-";
			}
            genGI = array[8];
            
			if(!array[9].equals("-")) { 
				start = Integer.parseInt(array[9]); 
			} else { 
				start = -1; 
			}

			if(!array[10].equals("-")) { 
				end = Integer.parseInt(array[10]);
			} else { 
				end = -1;
			}
            orientation = array[11].equals("+") ? 1 : -1;
            
            assembly = array[12];
        }
        
        public int hashCode() { 
            int code = 17;
            
            code += taxID.hashCode(); code *= 37;
            code += geneID.hashCode(); code *= 37;
            code += status.hashCode(); code *= 37;
            
            code += rnaAcc.hashCode(); code *= 37;
            code += rnaVersion.hashCode(); code *= 37;
            code += rnaGI.hashCode(); code *= 37;
            
            code += protAcc.hashCode(); code *= 37;
            code += protVersion.hashCode(); code *= 37;
            code += protGI.hashCode(); code *= 37;

            code += genAcc.hashCode(); code *= 37;
            code += genVersion.hashCode(); code *= 37;
            code += genGI.hashCode(); code *= 37;
            
            code += start; code *= 37;
            code += end; code *= 37;
            code += orientation; code *= 37;
            code += assembly.hashCode(); code *= 37;
            
            return code; 
        }
        
        public String toString() { 
            return "#" + geneID + " " + rnaAcc + "," + protAcc + "," + genAcc + " [" + taxID + "]";
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Entry)) { return false; }
            Entry e = (Entry)o;

            if(!taxID.equals(e.taxID)) { return false; }
            if(!geneID.equals(e.geneID)) { return false; }
            if(!status.equals(e.status)) { return false; }
            
            if(!rnaAcc.equals(e.rnaAcc)) { return false; }
            if(!rnaVersion.equals(e.rnaVersion)) { return false; }
            if(!rnaGI.equals(e.rnaGI)) { return false; }
            
            if(!genAcc.equals(e.genAcc)) { return false; }
            if(!genVersion.equals(e.genVersion)) { return false; }
            if(!genGI.equals(e.genGI)) { return false; }
            
            if(!protAcc.equals(e.protAcc)) { return false; }
            if(!protVersion.equals(e.protVersion)) { return false; }
            if(!protGI.equals(e.protGI)) { return false; }
            
            if(start != e.start) { return false; }
            if(end != e.end) { return false; }
            if(orientation != e.orientation) { return false; }
            if(!assembly.equals(e.assembly)) { return false; }
            
            return true;
        }
    }
}
