/*
 * Created on Aug 27, 2007
 */
package edu.mit.csail.cgs.tools.chippet;

import java.io.*;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.util.regex.*;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import java.sql.Connection;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class ElandParser implements Iterator<ElandRecord> {
    
    public static void main(String[] args) { 
        File f = new File(args[0]);
        Filter<ElandRecord,ElandRecord> filter = new Filter<ElandRecord,ElandRecord>() {
            public ElandRecord execute(ElandRecord a) {
                return a.matchCode.equals("U0") ? a : null;
            }
        };
        
        try {
            ElandParser parser = new ElandParser(f);
            Iterator<ElandRecord> itr = new FilterIterator<ElandRecord,ElandRecord>(filter, parser);
            
            if(args.length == 1) { 
                int c = 0;
                while(itr.hasNext()) { 
                    ElandRecord rec = itr.next();
                    c += 1;
                }
                System.out.println("# Unique Matches: " + c);
            } else { 
                String genomeName = args[1];
                String exptName = args[2];
                Genome genome = Organism.findGenome(genomeName);
                Connection cxn = DatabaseFactory.getConnection("chippet");
                
                insertIntoDB(itr, genome, cxn, exptName);
                
                DatabaseFactory.freeConnection(cxn);
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public static void insertIntoDB(Iterator<ElandRecord> parser,
            Genome genome, 
            java.sql.Connection cxn, 
            String exptName) throws SQLException {
        
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

        PreparedStatement ips = cxn.prepareStatement("insert into chippetdata " +
                "(expt, chromosome, startpos, stoppos, strand, peakOverlap) values (?, ?, ?, ?, ?, ?)");
        
        int c = 0;
        while(parser.hasNext()) { 
            ElandRecord rec = parser.next();
            rec.insertIntoDB(genome, exptID, ips);
            c += 1;
            if(c % 1000 == 0) { System.out.print("."); System.out.flush(); }
            if(c % 100000 == 0) { System.out.print("[" + (c / 1000) + "k]"); System.out.flush(); }
        }
        System.out.println();
        
        cxn.commit();
        cxn.setAutoCommit(true);
    }
    
    private File file;
    private BufferedReader br;
    private String nextLine;
    
    public ElandParser(File f) throws IOException { 
        file = f;
        br = new BufferedReader(new FileReader(f));
        retrieveNextLine();
    }
    
    private void retrieveNextLine() { 
        try {
            nextLine = br.readLine();
        } catch (IOException e1) {
            e1.printStackTrace();
            nextLine = null;
        }
        
        if(nextLine == null) { 
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public boolean hasNext() {
        return nextLine != null;
    }

    public ElandRecord next() {
        ElandRecord rec = new ElandRecord(nextLine);
        retrieveNextLine();
        return rec;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}

class ElandRecord {
    
    public static Pattern chromPatt;
    
    static { 
        chromPatt = Pattern.compile("chr([\\w\\d]+)\\..*");
    }
    
    public String name;
    public String seq;
    public String matchCode;

    public String targetName;
    public int start;
    public boolean strand;
    public Vector<String> modifiers;
    
    public ElandRecord(String line) { 
        String[] array = line.split("\\s+");
        name = array[0];
        if(name.charAt(0) == '>') { name = name.substring(1, name.length()); }
        seq = array[1];
        matchCode = array[2];

        modifiers = new Vector<String>();
		strand = true;
		start = -1;
		targetName = null;

		if(matchCode.startsWith("U")) { 
			Matcher m = chromPatt.matcher(array[6]);
			if(m.matches()) { targetName = m.group(1); } else { 
				throw new IllegalArgumentException("Illegal target: " + array[6]);
			}
			strand = array[7].equals("F");
			for(int i = 8; i < array.length; i++) { 
				if(!array[i].equals("..")) { 
					modifiers.add(array[i]);
				}
			}
		}
    }
    
    public void insertIntoDB(Genome g, int exptID, PreparedStatement ps) throws SQLException { 
        //PreparedStatement ips = cxn.prepareStatement("insert into chippetdata " +
        //"(expt, chromosome, startpos, stoppos, strand, peakOverlap) values (?, ?, ?, ?, ?, ?)");
        
        int chromid = g.getChromID(targetName);
        ps.setInt(1, exptID);
        ps.setInt(2, chromid);
        if(strand) { 
            ps.setInt(3, start);
            ps.setInt(4, start+seq.length()-1);
        } else { 
            ps.setInt(3, start-seq.length()+1);
            ps.setInt(4, start);
        }
        ps.setString(5, strand ? "+" : "-");
        ps.setInt(6, 1);
        
        ps.executeUpdate();
    }
}
