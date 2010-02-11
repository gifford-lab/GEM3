package edu.mit.csail.cgs.tools.chippet;

import java.util.*;
import java.sql.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.alignments.parsing.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLEntry;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLEntryPredicate;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLParser;

/**
 * Reads all the PSL files in a single directory, and analyzes them as a total set of 
 * hits.
 * 
 * @author tdanford
 */
public class PSLDirectoryAnalysis {
    
    public static void main(String[] args) { 
		if(args.length == 0) { 
			System.err.println("USAGE: PSLDirectoryAnalysis <psl-dir> [db_version]");
			System.exit(1);
		}

        File dir = new File(args[0]);
		int matches = 34;
		int gap = 10;
		int max_matches = 1;

        try {
            Genome g = Organism.findGenome("mm6");
            ChipPetBlatPredicate pred = new ChipPetBlatPredicate(matches, gap);
            PSLDirectoryAnalysis pda = new PSLDirectoryAnalysis(dir, g, pred, max_matches);
            
            File[] pslfiles = null;
			pslfiles = pda.getPSLFiles();
            
            for(int i = 0; i < pslfiles.length; i++) { 
                System.out.println("Adding " + pslfiles[i].getName());
                pda.addFile(pslfiles[i]);
            }
            
            System.out.println("Finished directory creation.");
            
            System.out.println("Optimizing...");
            pda.optimize();
            System.out.println("Finished optimization.");

            if(args.length > 1) { 
                String vname = args[1];
                System.out.println("Inserting into db \"" + vname + "\" (" + g.getVersion() + ")");
                java.sql.Connection cxn = DatabaseFactory.getConnection("chippet");
                pda.insertIntoDB(cxn, vname);
                DatabaseFactory.freeConnection(cxn);
                System.out.println("Finished.");
            }

        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

	private File directory;
	private Genome genome;
	private Map<String,ChromosomalBlatSummary> posChromSummaries, negChromSummaries;
    private BlatPSLEntryPredicate pred;
    
    private Map<String,Integer> counts;
    private int thresholdCount;
	
	public PSLDirectoryAnalysis(File d, Genome g, BlatPSLEntryPredicate p, int tc) {
        thresholdCount = tc;
        counts = new HashMap<String,Integer>();
        pred = p;
		directory = d;
		genome = g;
        posChromSummaries = new HashMap<String,ChromosomalBlatSummary>();
        for(String chrom : genome.getChromList()) { 
            posChromSummaries.put(chrom, new ChromosomalBlatSummary(genome, chrom, "+"));
        }
        negChromSummaries = new HashMap<String,ChromosomalBlatSummary>();
        for(String chrom : genome.getChromList()) { 
            negChromSummaries.put(chrom, new ChromosomalBlatSummary(genome, chrom, "-"));
        }
	}

	public void setThresholdCount(int tc) { 
		thresholdCount = tc;
	}
    
    public void insertIntoDB(java.sql.Connection cxn, String exptName) throws SQLException {
        cxn.setAutoCommit(false);
        PreparedStatement ps = 
            cxn.prepareStatement("insert into chippetexpt (id, name) values (chippetexpt_id.nextval, ?)");
        ps.setString(1, exptName);
        ps.executeUpdate();
        ps.close();
        
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "chippetexpt_id"));
        if(!rs.next()) { throw new IllegalStateException(); }
        int exptID = rs.getInt(1);
        rs.close();
        s.close();

        ps = cxn.prepareStatement("insert into chippetToGenome (expt, genome) values (?, ?)");
        ps.setInt(1, exptID);
        ps.setInt(2, genome.getDBID());
        ps.executeUpdate();
        ps.close();
        
        for(String c : posChromSummaries.keySet()) { 
            posChromSummaries.get(c).insertIntoDB(cxn, exptID);
        }
        
        for(String c : negChromSummaries.keySet()) { 
            negChromSummaries.get(c).insertIntoDB(cxn, exptID);
        }
        
        cxn.commit();
        cxn.setAutoCommit(true);
    }

	public long getSize() { 
		long s = 0;
        for(String c : posChromSummaries.keySet()) { 
            s += posChromSummaries.get(c).getSize();
        }
        for(String c : negChromSummaries.keySet()) { 
            s += negChromSummaries.get(c).getSize();
        }
		return s;
	}
    
    public void optimize() { 
        for(String c : posChromSummaries.keySet()) { 
            posChromSummaries.get(c).optimize();
        }
        for(String c : negChromSummaries.keySet()) { 
            negChromSummaries.get(c).optimize();
        }
    }
    
    public void printSummary(PrintStream ps) { 
        TreeSet<String> chroms = new TreeSet<String>(posChromSummaries.keySet());
        chroms.addAll(negChromSummaries.keySet());
        
        for(String c : chroms) { 
            System.out.println("Chrom: \"" + c + "\"");
            posChromSummaries.get(c).printSummary(ps);
            negChromSummaries.get(c).printSummary(ps);
        }
    }
    
    public Set<String> getKeys() { 
        Set<String> keys = new HashSet<String>();
        for(String chrom : posChromSummaries.keySet()) { 
            keys.addAll(posChromSummaries.get(chrom).getKeys());
        }
        for(String chrom : negChromSummaries.keySet()) { 
            keys.addAll(negChromSummaries.get(chrom).getKeys());
        }
        return keys;
    }
    
    public int getKeyCount(String k) { 
        int c = 0;
        for(String chrom : posChromSummaries.keySet()) { 
            c += posChromSummaries.get(chrom).getCount(k);
        }
        for(String chrom : negChromSummaries.keySet()) { 
            c += negChromSummaries.get(chrom).getCount(k);
        }
        return c;
    }
	
	public void addFile(File pslfile) throws IOException { 
		BlatPSLParser parser = new BlatPSLParser();
		Iterator<BlatPSLEntry> entries = parser.parse(pslfile, pred);
        
		Set<String> unknownChroms = new TreeSet<String>();
		long c = 0;
        long r = 0;
        
		while(entries.hasNext()) { 
			BlatPSLEntry entry = entries.next();
            String key = entry.getQname();
            
            if(!counts.containsKey(key)) { counts.put(key, 0); }
            counts.put(key, counts.get(key) + 1);

            if(counts.get(key) == thresholdCount + 1) { 
                for(String chrom : posChromSummaries.keySet()) { 
                    r += posChromSummaries.get(chrom).removeKeyedValue(key);
                }
                for(String chrom : negChromSummaries.keySet()) { 
                    r += negChromSummaries.get(chrom).removeKeyedValue(key);
                }
                
            } else if(counts.get(key) <= thresholdCount) { 
                String chrom = entry.getTname();
                if(chrom.startsWith("chr")) { 
                    chrom = chrom.substring(3, chrom.length());
                }
                char strand = entry.getStrand();
                
                if(strand == '-') { 
                    if(negChromSummaries.containsKey(chrom)) {
                        negChromSummaries.get(chrom).addEntry(entry);
                        c += 1;
                    } else { 
                        unknownChroms.add(chrom);
                    }                                        
                } else { 
                    if(posChromSummaries.containsKey(chrom)) {
                        posChromSummaries.get(chrom).addEntry(entry);
                        c += 1;
                    } else { 
                        unknownChroms.add(chrom);
                    }                    
                }
            }
		}

		if(unknownChroms.size() > 0) { 
		    System.out.println("File " + pslfile.getName() + " bad chroms:");
		    for(String ch : unknownChroms) { 
		        System.out.println("\t" + ch);
		    }
		}

		System.out.println("\tAdded: " + c + " entries.");
        System.out.println("\tRemoved: " + r + " entries.");
		System.out.println("\tTotal: " + getSize());
	}
	
	public File[] getPSLFiles() throws IOException {
		FilenameFilter filter = new FilenameFilter() {
			public boolean accept(File f, String n) {
				return n.toUpperCase().endsWith(".PSL");
			} 
			
		};
		String[] names = directory.list(filter);
		File[] files = new File[names.length];
		for(int i = 0; i < names.length; i++) { files[i] = new File(directory, names[i]); }
		return files;
	}
}
