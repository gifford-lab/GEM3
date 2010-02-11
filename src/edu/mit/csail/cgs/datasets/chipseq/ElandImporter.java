/*
 * Created on Jan 7, 2008
 */
package edu.mit.csail.cgs.datasets.chipseq;

import java.util.*;
import java.util.logging.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.io.parsing.alignment.*;

public class ElandImporter {

    public static void main(String[] args) { 
        ElandImporter importer = new ElandImporter();
        
        try {
            Genome genome = Organism.findGenome(args[0]);
            //importer.loadElandFile(new File(args[1]));
            
            String exptName = args[2];
            String repName = args[3];
            String alignmentName = args[4];
            String factorName = args[5];
            String cellsName = args[6];
            String condName = args[7];
            boolean append = (args.length>=9 && args[8].equals("append")) ? true : false;
            int startRead = args.length>=10 ? new Integer(args[9]).intValue() : 0;
            
            //importer.checkSequences(genome);
            
            if(append)
            	importer.appendToDatabase(new File(args[1]), exptName, repName, alignmentName, genome, factorName, condName, cellsName, startRead);
            else
            	importer.insertIntoDatabase(new File(args[1]), exptName, repName, alignmentName, genome, factorName, condName, cellsName);
         
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

    private Logger logger;
    private Set<ElandHit.Code> admissibleCodes;
    //private LinkedList<ElandHit> hits;

    public ElandImporter() {
        logger = Logger.getLogger("edu.mit.csail.cgs.datasets.chipseq.ElandImporter");
        
        admissibleCodes = new HashSet<ElandHit.Code>();
        admissibleCodes.add(ElandHit.Code.U0);
        admissibleCodes.add(ElandHit.Code.U1);
    //    hits = new LinkedList<ElandHit>();
    }

    /*public void loadElandFile(File elf) throws IOException { 
        ElandFile eland = new ElandFile(elf, admissibleCodes);
        while(eland.hasNext()) { 
            ElandHit hit = eland.next();
            hits.addLast(hit);
        }
        
        logger.log(Level.INFO, String.format("Loaded %d Eland Hits.", hits.size()));
    }
    
    public void checkSequences(Genome g) throws SQLException {
        for(ElandHit hit : hits) { 
            int start = hit.getCoordinate();
            int end = start + hit.getQuerySequence().length() - 1;
            boolean strand = hit.getStrand();
            String chrom = hit.getTargetSequence();
            Genome.ChromosomeInfo info = g.getChrom(chrom);
            String seq = g.getChromosomeSequence(info, start, end);
            if(!strand) { seq = revcomp(seq); }
            
            System.out.println(String.format("ELAND: %s (%s,%s)", 
                    hit.getQuerySequence(), hit.getCode().toString(), String.valueOf(strand)));
            System.out.println(String.format("DB:    %s\n", seq));
        }
    }*/
    
    private String revcomp(String str) { 
        StringBuilder sb = new StringBuilder();
        for(int i = str.length()-1; i >= 0; i--) { 
            sb.append(complement(str.charAt(i)));
        }
        return sb.toString();
    }
    
    private char complement(char c) { 
        switch(c) { 
        case 'a':
            return 't';
        case 'A':
            return 'T';
        case 't':
            return 'a';
        case 'T':
            return 'A';
        case 'g':
            return 'c';
        case 'G':
            return 'C';
        case 'c':
            return 'g';
        case 'C':
            return 'G';
        default:
            return c;
        }
    }

    public void insertIntoDatabase(File elf, String exptName, String repName, 
    		String alignmentName, 
            Genome g, String factorName, String condName, String cellsName) 
    throws SQLException, IOException { 

    	ElandFile eland = new ElandFile(elf, admissibleCodes);
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
            ps.setString(2, repName);
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
            
            PreparedStatement ps2 = cxn.prepareStatement(String.format(
                    "insert into chipseqhits " +
                    "(read, expt, alignment, chromosome, startpos, stoppos, strand) values" +
            " (?, ?, ?, ?, ?, ?, ?)"));
            
            int c = 0;
            
            int hi = 0;
            while(eland.hasNext()) { 
                ElandHit hit = eland.next();
                ps.setInt(1, exptID);
                ps.setString(2, hit.getName());
                ps.setString(3, hit.getQuerySequence());
                ps.executeUpdate();
                                                
                rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "chipseqread_id"));
                rs.next();
                int readID = rs.getInt(1);
                rs.close(); rs = null;
                
    
                
                ps2.setInt(1, readID);
                ps2.setInt(2, exptID);
                ps2.setInt(3, alignID);

                String chrom = hit.getTargetSequence();
                int chromID = g.getChromID(chrom);
                
                String strand = hit.getStrand() ? "+" : "-";
                int start = hit.getCoordinate();
                int end = start + hit.getQuerySequence().length() - 1;

                if(chromID != -1) { 
                    ps2.setInt(4, chromID);

                    ps2.setInt(5, start);
                    ps2.setInt(6, end);

                    ps2.setString(7, strand);
                    ps2.executeUpdate();
                } else { 
                    logger.log(Level.SEVERE, String.format(
                            "Couldn't find chromosome %s for hit %s", 
                            chrom, hit.getName()));
                }
                hi += 1;
                if(hi % 1000000 == 0) { //Commits every 1M reads
                	cxn.commit();
                	System.out.print("(Committed: " + hi + ")"); System.out.flush();
                }else if(hi % 100000 == 0) { System.out.print("(" + hi + ")"); System.out.flush();}
                else if(hi % 10000 == 0) { System.out.print("."); System.out.flush(); }
                
                c += 1;
            }               
            System.out.println();
                       
            ps.close(); ps = null;
            ps2.close(); ps2 = null;

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
    
    public void appendToDatabase(File elf, String exptName, String repName, 
    		String alignmentName, 
            Genome g, String factorName, String condName, String cellsName, int startRead) 
    throws SQLException, IOException { 

    	ElandFile eland = new ElandFile(elf, admissibleCodes);
        java.sql.Connection cxn = null;
        Statement s = null;
        PreparedStatement ps = null;
        ResultSet rs = null;
        MetadataLoader ml = new MetadataLoader();
        logger.log(Level.INFO, String.format("Append starting from valid read number %d", startRead));
        
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

            /*************************************************
            logger.log(Level.INFO, "Creating Experiment...");
            ps = cxn.prepareStatement(String.format(
                    "insert into chipseqexpts " +
                    "(id, name, replicate, species, cells, condition, factor) values " +
                    "(%s, ?, ?, ?, ?, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqexpt_id")));
            ps.setString(1, exptName);
            ps.setString(2, repName);
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
			**************************/
            //Lookup experiment
            int exptID = -1; 
            ps = cxn.prepareStatement("select id from chipseqexpts " +
            		"where name=? and replicate=?");
            ps.setString(1, exptName);
            ps.setString(2, repName);
            rs = ps.executeQuery();
            if(rs.next()) { 
            	exptID = rs.getInt(1);
            } 
            rs.close(); rs = null;
            ps.close(); ps = null;

            int alignID = -1;
            ps = cxn.prepareStatement("select id from chipseqalignments " +
            		"where name=? and expt=? and genome=?");
            ps.setString(1, alignmentName);
            ps.setInt(2, exptID);
            ps.setInt(3, g.getDBID());
            rs = ps.executeQuery();
            if(rs.next()) { 
            	alignID = rs.getInt(1);
            }
            rs.close(); rs = null;
            ps.close(); ps = null;
            
            logger.log(Level.INFO, "Inserting Reads...");

            ps = cxn.prepareStatement(String.format(
                    "insert into chipseqreads (id, expt, name, sequence) values " +
                    "(%s, ?, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqread_id")));
            
            PreparedStatement ps2 = cxn.prepareStatement(String.format(
                    "insert into chipseqhits " +
                    "(read, expt, alignment, chromosome, startpos, stoppos, strand) values" +
            " (?, ?, ?, ?, ?, ?, ?)"));
            
            int c = 0;
            
            int hi = 0;
            while(eland.hasNext()) { 
                ElandHit hit = eland.next();
                if(hi>startRead){
                    ps.setInt(1, exptID);
	                ps.setString(2, hit.getName());
	                ps.setString(3, hit.getQuerySequence());
	                ps.executeUpdate();
	                                                
	                rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "chipseqread_id"));
	                rs.next();
	                int readID = rs.getInt(1);
	                rs.close(); rs = null;
	                
	    
	                
	                ps2.setInt(1, readID);
	                ps2.setInt(2, exptID);
	                ps2.setInt(3, alignID);
	
	                String chrom = hit.getTargetSequence();
	                int chromID = g.getChromID(chrom);
	                
	                String strand = hit.getStrand() ? "+" : "-";
	                int start = hit.getCoordinate();
	                int end = start + hit.getQuerySequence().length() - 1;
	
	                if(chromID != -1) { 
	                    ps2.setInt(4, chromID);
	
	                    ps2.setInt(5, start);
	                    ps2.setInt(6, end);
	
	                    ps2.setString(7, strand);
	                    ps2.executeUpdate();
	                } else { 
	                    logger.log(Level.SEVERE, String.format(
	                            "Couldn't find chromosome %s for hit %s", 
	                            chrom, hit.getName()));
	                }
                }
                hi += 1;
                if(hi % 1000000 == 0 && hi>startRead) { //Commits every 1M reads
                	cxn.commit();
                	System.out.print("(Committed: " + hi + ")"); System.out.flush();
                }else if(hi % 100000 == 0) { System.out.print("(" + hi + ")"); System.out.flush();}
                else if(hi % 10000 == 0) { System.out.print("."); System.out.flush(); }
                
                c += 1;
            }               
            System.out.println();
                       
            ps.close(); ps = null;
            ps2.close(); ps2 = null;

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
/*
    public void appendToDatabase(String exptName, String repName, String alignmentName, 
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

           

            int exptID = -1;
            
           
            
            ps = cxn.prepareStatement("select id from chipseqexpts " +
            		"where name=? and replicate=?");
            ps.setString(1, exptName);
            ps.setString(2, repName);
            rs = ps.executeQuery();
            if(rs.next()) { 
            	exptID = rs.getInt(1);
            } 
            rs.close(); rs = null;
            ps.close(); ps = null;

           

            int alignID = -1;
          
            ps = cxn.prepareStatement("select id from chipseqalignments " +
            		"where name=? and expt=? and genome=?");
            ps.setString(1, alignmentName);
            ps.setInt(2, exptID);
            ps.setInt(3, g.getDBID());
            rs = ps.executeQuery();
            if(rs.next()) { 
            	alignID = rs.getInt(1);
            }
            rs.close(); rs = null;
            ps.close(); ps = null;

            logger.log(Level.INFO, "Inserting Reads...");

            ps = cxn.prepareStatement(String.format(
                    "insert into chipseqreads (id, expt, name, sequence) values " +
                    "(%s, ?, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqread_id")));
            int c = 0;
            
            for(ElandHit hit : hits) { 
                ps.setInt(1, exptID);
                ps.setString(2, hit.getName());
                ps.setString(3, hit.getQuerySequence());
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
            
            int hi = 0;
            System.out.println(String.format("[%d]", hits.size()));
            for(ElandHit hit : hits) { 
                if(readIDs.containsKey(hit.getName())) { 
                    ps.setInt(1, readIDs.get(hit.getName()));
                    ps.setInt(2, exptID);
                    ps.setInt(3, alignID);

                    String chrom = hit.getTargetSequence();
                    int chromID = g.getChromID(chrom);
                    
                    String strand = hit.getStrand() ? "+" : "-";
                    int start = hit.getCoordinate();
                    int end = start + hit.getQuerySequence().length() - 1;

                    if(chromID != -1) { 
                        ps.setInt(4, chromID);

                        ps.setInt(5, start);
                        ps.setInt(6, end);

                        ps.setString(7, strand);
                        ps.executeUpdate();
                    } else { 
                        logger.log(Level.SEVERE, String.format(
                                "Couldn't find chromosome %s for hit %s", 
                                chrom, hit.getName()));
                    }

                } else { 
                    logger.log(Level.SEVERE, String.format(
                            "Couldn't find hits for query %s", hit.getName()));
                }

                hi += 1;
                if(hi % 10000 == 0) { System.out.print("."); System.out.flush(); }
                if(hi % 100000 == 0) { System.out.print("(" + hi + ")"); System.out.flush(); }
            }
            System.out.println();
                       
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
    }*/
}
