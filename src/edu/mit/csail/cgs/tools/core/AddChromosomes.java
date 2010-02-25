package edu.mit.csail.cgs.tools.core;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.Pair;
import java.sql.*;
import java.util.Properties;
import java.io.*;

/**
 * Reads a FASTA file and creates the chromosomes and their sequence in the database.
 * Leading "chr" is stripped off names in the FASTA file:
 *
 * AddChromosomes --species "$MM;mm8" --fasta mm8.fa
 * AddChromosomes --species "$MM;mm8" < mm8.fa
 *
 *  --incremental does writes in little chunks if your PreparedStatement.setString()
 * can't handle a full chromosome sequence at once.
 */

public class AddChromosomes {
	
	public static void main(String[] args) { 
		try { 
			boolean incrementalLoad = Args.parseFlags(args).contains("incremental");
			
			if(!incrementalLoad) {
				System.out.println("-> Bulk loading");
				bulk_load(args);
			} else { 
				System.out.println("-> Incremental Loading");
				incremental_load(args);
			}
			
		} catch(IOException e) { 
			e.printStackTrace(System.err);
		} catch(SQLException e) { 
			e.printStackTrace(System.err);
		} catch(NotFoundException e) { 
			e.printStackTrace(System.err);
		}
	}

    public static void incremental_load(String[] args) throws IOException, NotFoundException, SQLException {
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        Genome genome = pair.getLast();
        String fastaname = Args.parseString(args,"fasta",null);
        
        BufferedReader reader = null;
        FASTAStream fasta;
        if (fastaname == null) {
            reader = new BufferedReader(new InputStreamReader(System.in));
        } else {
            reader = new BufferedReader(new FileReader(new File(fastaname)));
        }
        java.sql.Connection cxn = DatabaseFactory.getConnection("core");        

        boolean ac = cxn.getAutoCommit();
        cxn.setAutoCommit(false);
        
        PreparedStatement insertChrom = 
        	cxn.prepareStatement("insert into chromosome(id, genome, name) values (" +
        			Sequence.getInsertSQL(cxn, "chromosome_id") + "," +
        			genome.getDBID() + ",?)");
        
        PreparedStatement getChromID = 
            cxn.prepareStatement("select id from chromosome where genome = ? and name = ?");
        
        PreparedStatement insertSequence = 
        	cxn.prepareStatement("insert into chromsequence(id,sequence) values (?,'')");
        
        PreparedStatement appendSequence = 
        	cxn.prepareStatement("update chromsequence " +
        			"set sequence=CONCAT(sequence,?) where id=?");

        String line = null;
        int lastChromID = -1;
        int length = 0;
        
        int hashTick = 100000;
        int nextHash = 0;
        
        while((line = reader.readLine()) != null) { 
        	line = line.trim();
        	if(line.length() > 0) {
        		if(line.startsWith(">")) { 
        			String name = line.substring(1, line.length()).trim();
        			name = name.replaceFirst("^chr","");        			
        			
                    try {
                        insertChrom.setString(1,name);
                        insertChrom.execute();
                    } catch (SQLException e) {
                        System.err.println("Didn't add chromosome " + name + ".  Going to try sequence");
                    }
                    getChromID.setInt(1,genome.getDBID());
                    getChromID.setString(2,name);
            		ResultSet rs = getChromID.executeQuery();
            		rs.next();
            		lastChromID = rs.getInt(1);
            		rs.close();

            		insertSequence.setInt(1,lastChromID);
            		insertSequence.execute();
            		
            		length = 0;
            		nextHash = hashTick;
            		System.out.println();
            		System.out.print(String.format("\"%s\"", name));
            		System.out.flush();
            		
        		} else { 
        			String seq = line;

        			appendSequence.setString(1, seq);
        			appendSequence.setInt(2, lastChromID);
        			appendSequence.executeUpdate();
        			
        			length += seq.length();
        			
        			if(length-seq.length() < nextHash && length >= nextHash) { 
        				int lower = length-(length%hashTick);
        				System.out.print(String.format(" %dk", lower/1000));
        				System.out.flush();
        				nextHash = lower + hashTick;
        			}
        		}
        	}
        }
        
        System.out.println("\n-> Finished.");
        System.out.println("-> Committing to DB.");

        cxn.commit();
        cxn.setAutoCommit(ac);
    }

    public static void bulk_load(String args[]) throws IOException, NotFoundException, SQLException {
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        Genome genome = pair.getLast();
        String fastaname = Args.parseString(args,"fasta",null);
        FASTAStream fasta;
        if (fastaname == null) {
            fasta = new FASTAStream(new BufferedReader(new InputStreamReader(System.in)));
        } else {
            fasta = new FASTAStream(new File(fastaname));
        }
        java.sql.Connection cxn = DatabaseFactory.getConnection("core");        

        boolean ac = cxn.getAutoCommit();
        cxn.setAutoCommit(false);
        
        PreparedStatement insertChrom = cxn.prepareStatement("insert into chromosome(id, genome, name) values (" +
                                                    Sequence.getInsertSQL(cxn, "chromosome_id") + "," +
                                                    genome.getDBID() + ",?)");
        PreparedStatement getChromID = cxn.prepareStatement("select id from chromosome where genome = ? and name = ?");
        PreparedStatement insertSequence = cxn.prepareStatement("insert into chromsequence(id,sequence) values (?,?)");
        PreparedStatement getSequence = cxn.prepareStatement("select sequence from chromsequence where id = ?");
        while (fasta.hasNext()) {
            Pair<String,String> next = fasta.next();
            String name = next.car();
            name = name.replaceFirst("^chr","");
            String seq = next.cdr();
            try {
                insertChrom.setString(1,name);
                insertChrom.execute();
            } catch (SQLException e) {
                System.err.println("Didn't add chromosome " + name + ".  Going to try sequence");
            }
            getChromID.setInt(1,genome.getDBID());
            getChromID.setString(2,name);
            ResultSet rs = getChromID.executeQuery();
            rs.next();
            int id = rs.getInt(1);
            rs.close();
            insertSequence.setInt(1,id);
            edu.mit.csail.cgs.utils.database.ClobHandler.setClob(cxn, insertSequence, 2, seq);
            insertSequence.execute();
        }                
        cxn.commit();
        cxn.setAutoCommit(ac);
    }

}