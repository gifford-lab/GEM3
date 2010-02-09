package edu.mit.csail.cgs.utils.io.parsing.alignment;

import java.io.*;
import java.util.*;

/**
 * This is a class that parses the BLAT commands PSL output.
 */
public class BlatPSLParser {
    
    public static void main(String[] args) { 
        //simple_parsing(args);
        iterative_parsing(args);
    }
    
    public static void iterative_parsing(String[] args) {
        File inputFile = new File(args[0]);
        BlatPSLParser parser = new BlatPSLParser();
        try {
            Iterator<BlatPSLEntry> entries = parser.parse(inputFile);
            
            Map<String,Set<BlatPSLEntry>> matchMap = new HashMap<String,Set<BlatPSLEntry>>();
            Comparator<BlatPSLEntry> comparator = new BlatPSLEntryMatchComparator();

			int threshMatch = 30;
            
            while(entries.hasNext()) { 
                BlatPSLEntry e = entries.next();
                if(e.getMatch() >= threshMatch) { 
                    if(!matchMap.containsKey(e.getQname())) { 
                        matchMap.put(e.getQname(), new TreeSet<BlatPSLEntry>(comparator));
                    }
                    matchMap.get(e.getQname()).add(e);
                }
            }
            
            for(String query : matchMap.keySet()) { 
                System.out.println("\n" + query);
                for(BlatPSLEntry e : matchMap.get(query)) { 
                    System.out.println("\t" + summarizeBlatEntry(e));
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

	public static void simple_parsing(String[] args) { 
		BlatPSLParser parser = new BlatPSLParser();
		try { 
			Collection<BlatPSLEntry> entries = parser.parseAll(new File(args[0]));
			
			Map<String,Set<BlatPSLEntry>> summary = 
				new HashMap<String,Set<BlatPSLEntry>>();
			for(BlatPSLEntry e : entries) { 
				if(!summary.containsKey(e.getQname())) { 
					summary.put(e.getQname(), new HashSet<BlatPSLEntry>());
				}
				summary.get(e.getQname()).add(e);
			}
			
			Comparator<BlatPSLEntry> comparator = new BlatPSLEntryMatchComparator();
			for(String query : summary.keySet()) { 
				TreeSet<BlatPSLEntry> qentries = new TreeSet<BlatPSLEntry>(comparator);
				qentries.addAll(summary.get(query));

				System.out.println("\n" + query + " ");
				for(BlatPSLEntry e : qentries) {
                    System.out.println("\t" + summarizeBlatEntry(e));
				}
			}
			
		} catch(IOException e) { 
			e.printStackTrace(System.err);
		}
	}
    
    private static String summarizeBlatEntry(BlatPSLEntry e) {  
        String str =  e.getMatch() + "/" + e.getMismatch() + " (" + e.getQsize() + 
                ")  bc:" + e.getBlockCount() + 
				" qgap#:" + e.getQgapCount() + 
                " tgap#:" + e.getTgapCount();
		if(e.getTgapCount() > 0) { 
			str += " (";
			for(int i = 0; i < e.getTgapCount(); i++) { 
				if(i > 0) { str += ","; }	
				str += String.valueOf(e.getTGapSize(i));
			}
			str += ")";
		}
		str += " " + e.getTname() + ":" + e.getTstart() + "-" + e.getTend();
		return str;
    }

	private static class BlatPSLEntryMatchComparator implements Comparator<BlatPSLEntry> { 
		public int compare(BlatPSLEntry e1, BlatPSLEntry e2) { 
			int m1 = e1.getMatch(), m2 = e2.getMatch();
			if(m1 > m2) { return -1; }
			if(m1 < m2) { return 1; }
            int c = e1.getTname().compareTo(e2.getTname());
			if(c != 0) { return c; }
			if(e1.getTstart() < e2.getTstart()) { return -1; }
			if(e1.getTstart() > e2.getTstart()) { return 1; }
			if(e1.getTend() < e2.getTend()) { return -1; }
			if(e1.getTend() > e2.getTend()) { return 1; }
			return e1.getQname().compareTo(e2.getQname());
		}
	}
    
    private boolean skipHeader;
	
	public BlatPSLParser() { 
	    skipHeader = true;
	}

	public BlatPSLParser(boolean skiph) { 
		skipHeader = skiph;
	}
    
    public Collection<BlatPSLEntry> parseAll(File f) throws IOException { 
        return parseAll(f, null);
    }
	
	public Collection<BlatPSLEntry> parseAll(File f, BlatPSLEntryPredicate p) 
		throws IOException { 

		LinkedList<BlatPSLEntry> entries = 
			new LinkedList<BlatPSLEntry>();
		
		String line;
		BufferedReader br = new BufferedReader(new FileReader(f));
		
		if(skipHeader) { 
		    for(int i = 0; i < 5; i++) { 
		        br.readLine();
		    }
		}
		
		while((line = br.readLine()) != null) { 
			BlatPSLEntry entry = new BlatPSLEntry(line);
            if(p == null || p.acceptsEntry(entry)) { 
                entries.addLast(entry);
            }
		}

		br.close();
		return entries;
	}
	
    public Iterator<BlatPSLEntry> parse(File f) throws IOException { 
        return new BlatParsingIterator(f, null);
    }
    
	public Iterator<BlatPSLEntry> parse(File f, BlatPSLEntryPredicate p) throws IOException { 
		return new BlatParsingIterator(f, p);
	}
	
	private class BlatParsingIterator 
		implements Iterator<BlatPSLEntry> {
        
        private BufferedReader br;
        private BlatPSLEntry nextEntry;
        private BlatPSLEntryPredicate pred;
		
		public BlatParsingIterator(File f, BlatPSLEntryPredicate p) throws IOException {
            pred = p;
            br = new BufferedReader(new FileReader(f));
            if(skipHeader) { 
                for(int i = 0; i < 5; i++) { 
                    br.readLine();
                }
            }
            findNextEntry();
		}
        
        private void findNextEntry() throws IOException { 
            String line;
            nextEntry = null;
            
            while((line = br.readLine()) != null) { 
                BlatPSLEntry e = new BlatPSLEntry(line);
                if(pred == null || pred.acceptsEntry(e)) {
                    nextEntry = e;
                    break;
                }
            }
        }

		public boolean hasNext() {
            return nextEntry != null;
		}

		public BlatPSLEntry next() {
            BlatPSLEntry e = nextEntry;
            try {
                findNextEntry();
            } catch (IOException e1) {
                e1.printStackTrace();
                nextEntry = null;
            }
            return e;
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
