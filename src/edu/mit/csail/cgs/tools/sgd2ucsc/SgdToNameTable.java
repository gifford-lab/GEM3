package edu.mit.csail.cgs.tools.sgd2ucsc;

import java.util.*;
import java.io.*;
import java.sql.*;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.*;
import org.biojava.utils.ParserException;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;

public class SgdToNameTable {
    
    public static void main(String[] args) { 
        File inputFile = new File(args[0]);
        SgdToNameTable tableFiller = new SgdToNameTable(inputFile);

        try {
            Genome g = Organism.findGenome(args[1]);
            tableFiller.populateTable(g);
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private File inputFile;
    
    public SgdToNameTable(File ifile) { 
        inputFile = ifile;
    }
    
    public void populateTable(Genome g) throws SQLException, IOException { 
        java.sql.Connection cxn = g.getUcscConnection();
        try { 
            populateTable(cxn);
        } catch(IOException ie) { 
            throw ie;
        } catch(SQLException ie) { 
            throw ie;
        } finally { 
            DatabaseFactory.freeConnection(cxn);
        }
    }
    
    public void populateTable(Connection cxn) throws SQLException, IOException { 
        try {
            GFFEntrySet gffEntries = GFFTools.readGFF(inputFile);

            Statement s = cxn.createStatement();
            s.executeUpdate("delete from sgdToName");
            s.close();

            cxn.setAutoCommit(false);
            PreparedStatement ps = cxn.prepareStatement("insert into sgdToName (name, value) values (?, ?)");
            
            Iterator itr = gffEntries.lineIterator();
            int count = 0;
            while(itr.hasNext()) { 
                Object val = itr.next();
                if(val instanceof GFFRecord) { 
                    GFFRecord rec = (GFFRecord)val;
                    if(rec.getFeature().endsWith("gene")) { 
                        Map<String,List<String>> attrs = SGDGFFParser.decodeAttrMap(rec);
                        
                        if(attrs.containsKey("ID") && attrs.containsKey("Name")) { 
                            String id = attrs.get("ID").get(0);
                            String name = attrs.get("Name").get(0);

                            ps.setString(1, id);
                            ps.setString(2, name);
                            ps.executeUpdate();
                            count += 1;
                        }
                    }
                }
            }

            cxn.commit();
            ps.close();
            cxn.setAutoCommit(true);
            
            System.out.println("Entered " + count + " name/value pairs into sgdToName.");
            
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ParserException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }
        
    }
}


