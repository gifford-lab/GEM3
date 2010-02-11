package edu.mit.csail.cgs.datasets.chipseq;

import java.util.*;
import java.util.regex.*;
import java.util.logging.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.io.parsing.*;
import edu.mit.csail.cgs.utils.io.parsing.alignment.PSL;
import edu.mit.csail.cgs.utils.io.parsing.alignment.PSLHit;

/**
 * Usage:
 * java edu.mit.csail.cgs.datasets.chipseq.PSLImport reads.fasta reads.against_SGDv1.psl "Sc Gcn4, sigma, YPD;7/2/08" "PSL > 32 matches" "SGDv1" \
 * "Gcn4 (myc)" "YPD" "sigma 1278b"
 *
 */
public class PSLImporter {
	
	public static void main(String[] args) { 
		System.out.println("Importing...");
		File fasta = new File(args[0]);
		File psl = new File(args[1]);
		String exptName = args[2];
        String replicateName;
        String pieces[] = exptName.split(";");
        exptName = pieces[0];
        replicateName = pieces[1];
		String alignmentName = args[3];
		String genome = args[4];
		Genome g = null;
		try {
			g = Organism.findGenome(genome);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		
		if(g==null) { System.exit(1); }
		
		String fname = args[5];
		String condname = args[6];
		String cellsname = args[7];
		
		PSLImporter importer = new PSLImporter();
		try {
			importer.loadFastaFile(fasta);
			importer.loadPSLFile(psl);
			importer.insertIntoDatabase(exptName, replicateName,                                        
                                        alignmentName, g, 
                                        fname, condname, cellsname);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		System.out.println("Done: " + importer.seqMap.size() + " reads.");
	}

	private static Pattern chromPattern = Pattern.compile("chr(.*)");
	
	public static String cleanChromName(String c) { 
		Matcher m = chromPattern.matcher(c);
		if(m.matches()) { 
			return m.group(1);
		} else { 
			return c;
		}
	}
	
	private Logger logger;
	private Map<String,Set<PSLHit>> hitMap;
	private Map<String,String> seqMap;

	public PSLImporter() {
		logger = Logger.getLogger("edu.mit.csail.cgs.datasets.chipseq.PSLImporter");
	}
	
	public void insertIntoDatabase(String exptName, String replicateName, String alignmentName, 
                                   Genome g, String factorName, String condName, String cellsName) 
		throws SQLException { 
		
		java.sql.Connection cxn = null;
		Statement s = null;
		PreparedStatement ps = null;
		ResultSet rs = null;
		MetadataLoader ml = new MetadataLoader();

		try { 
			Organism species = null;
			try {
				species = Organism.getOrganism(g.getName());
			} catch (NotFoundException e) {
				e.printStackTrace();
				throw new IllegalStateException("Couldn't find organism: " + g.getName());
			}
			Factor factor = ml.getFactor(factorName);
			Condition condition = ml.getCondition(condName);
			Cells cells = ml.getCells(cellsName);

			cxn = DatabaseFactory.getConnection("chipseq");		
			cxn.setAutoCommit(false);
			s = cxn.createStatement();

			logger.log(Level.INFO, "Creating Experiment...");
			ps = cxn.prepareStatement(String.format(
					"insert into chipseqexpts " +
					"(id, name, replicate, species, cells, condition, factor) values " +
					"(%s, ?, ?, ?, ?, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqexpt_id")));
			ps.setString(1, exptName);
            ps.setString(2,replicateName);
			ps.setInt(3, species.getDBID());
			ps.setInt(4, cells.getDBID());
			ps.setInt(5, condition.getDBID());
			ps.setInt(6, factor.getDBID());
			ps.executeUpdate();
			ps.close(); ps = null;

			int exptID = -1;
			rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "chipseqexpt_id"));
			rs.next();
			exptID = rs.getInt(1);
			rs.close(); rs = null;

			logger.log(Level.INFO, "Creating Alignment...");
			ps = cxn.prepareStatement(String.format(
					"insert into chipseqalignments (id, expt, name, genome) " +
					"values (%s, ?, ?, ?)", 
					Sequence.getInsertSQL(cxn, "chipseqalignment_id")));
			ps.setInt(1, exptID);
			ps.setString(2, alignmentName);
			ps.setInt(3, g.getDBID());
			ps.executeUpdate();
			ps.close(); ps = null;

			int alignID = -1;
			rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "chipseqalignment_id"));
			rs.next();
			alignID = rs.getInt(1);
			rs.close(); rs = null;

			logger.log(Level.INFO, "Inserting Reads...");

			ps = cxn.prepareStatement(String.format(
					"insert into chipseqreads (id, expt, name, sequence) values " +
					"(%s, ?, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqread_id")));
			int c = 0;
			for(String qname : seqMap.keySet()) { 
				String seq = seqMap.get(qname);
				ps.setInt(1, exptID);
				ps.setString(2, qname);
				ps.setString(3, seq);
				ps.executeUpdate();
				c += 1;
			}
			ps.close(); ps = null;
			logger.log(Level.INFO, String.format("Loaded %d reads.", c));

			logger.log(Level.INFO, "Building READ ID map...");

			Map<String,Integer> readIDs = new HashMap<String,Integer>();
			rs = s.executeQuery("select id, name from chipseqreads where expt=" + exptID);
			while(rs.next()) { 
				int readID = rs.getInt(1);
				String readName = rs.getString(2);
				readIDs.put(readName, readID);
			}
			rs.close(); rs = null;
			logger.log(Level.INFO, String.format("%d Reads have assigned IDs.", readIDs.size()));
			
			logger.log(Level.INFO, "Inserting Hits...");
			
			ps = cxn.prepareStatement(String.format(
					"insert into chipseqhits " +
					"(read, expt, alignment, chromosome, startpos, stoppos, strand) values" +
                    " (?, ?, ?, ?, ?, ?, ?)"));
			for(String qname : hitMap.keySet()) { 
				if(readIDs.containsKey(qname)) { 
					for(PSLHit hit : hitMap.get(qname)) { 
						ps.setInt(1, readIDs.get(qname));
						ps.setInt(2, exptID);
						ps.setInt(3, alignID);

						String chrom = cleanChromName(hit.tname);
						int chromID = g.getChromID(chrom);

						if(chromID != -1) { 
							ps.setInt(4, chromID);

							ps.setInt(5, hit.tstart-1);
							ps.setInt(6, hit.tend-2);

							ps.setString(7, String.valueOf(hit.strand));
                            try {
                                ps.executeUpdate();
                            } catch (SQLException e) {
                                if (e.getErrorCode() == 1) {
                                    System.err.println("Ignoring duplicate primary key " + e.toString());
                                } else {
                                    System.err.println(e.toString());
                                    e.printStackTrace();
                                }
                            }
						} else { 
							logger.log(Level.SEVERE, String.format(
									"Couldn't find chromosome %s for hit %s", 
									chrom, qname));
						}
					}
				} else { 
					logger.log(Level.SEVERE, String.format(
							"Couldn't find hits for query %s", qname));
				}
			}
			ps.close(); ps = null;

			s.close(); s = null;
			cxn.commit();
		} catch(SQLException excep) { 
			cxn.rollback();
			throw excep;
		} finally { 
			if(s != null) { s.close(); }
			if(ps != null) { ps.close(); } 
			if(rs != null) { rs.close(); }
			ml.close();
			DatabaseFactory.freeConnection(cxn);
		}
	}
	
	public void loadFastaFile(File f) throws IOException {
		seqMap = new HashMap<String,String>();
		FASTAStream fasta = new FASTAStream(f);
		while(fasta.hasNext()) { 
			Pair<String,String> p = fasta.next();
			String header = p.getFirst();
			String[] array = header.split("\\s+");
			String qname = array[0];
			if(qname.startsWith(">")) { qname = qname.substring(1, qname.length()); }
			String qseq = p.getLast();
			seqMap.put(qname, qseq);
		}
		fasta.close();
	}
	
	public void loadPSLFile(File f) throws IOException { 
		hitMap = new HashMap<String,Set<PSLHit>>();
		PSL psl = new PSL(f);
		while(psl.hasNext()) { 
			PSLHit hit = psl.next();
			if(hit.blockSizes.length==1 && 
					(double)hit.matchSize() / (double)hit.qsize >= 0.9) {
				
				if(!hitMap.containsKey(hit.qname)) { 
					hitMap.put(hit.qname, new HashSet<PSLHit>());
				}
				hitMap.get(hit.qname).add(hit);
			}
		}
		
		for(String qname : hitMap.keySet()) { 
			PSLHit maxHit = null;
			double maxIdent = 0.0;
			int i = 0;
			for(PSLHit hit : hitMap.get(qname)) { 
				double ident = hit.getQueryIdentity();
				if(ident > maxIdent) { 
					maxIdent = ident;
					maxHit = hit;
				}
			}
			if(maxIdent >= 0.95) { 
				Iterator<PSLHit> hits = hitMap.get(qname).iterator();
				while(hits.hasNext()) { 
					PSLHit hit = hits.next();
					if(hit.getQueryIdentity() < 0.95) { 
						hits.remove();
					}
				}
			} else { 
				Set<PSLHit> newhits = new HashSet<PSLHit>();
				newhits.add(maxHit);
				hitMap.put(qname, newhits);
			}
		}
	}

	
}
