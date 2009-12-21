/*
 * Created on Jan 26, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.sgd2ucsc;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class ExamineUCSC {
    
    public static BinCalculator bincalc = new BinCalculator();
    
    public static void main(String[] args) {
        Genome g = null;
        try {
            g = Organism.findGenome("SGDv1");
            examineSGDTables(g);
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
        //printExonCounts(g, "refGene");
    }
    
    public static void printExonCounts(Genome g, String tblName) {
        Map<Integer,Integer> counts = new HashMap<Integer,Integer>();
        int maxCount = 0;
        RefGeneGenerator<Region> gen = new RefGeneGenerator<Region>(g, tblName);
        Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
        while(chroms.hasNext()) { 
            Region chrom = chroms.next();
            Iterator<Gene> genes = gen.execute(chrom);
            while(genes.hasNext()) {
                Gene gene = genes.next();
                if(gene instanceof ExonicGene) { 
                    ExonicGene exg = (ExonicGene)gene;
                    int count = exg.getNumExons();
                    if(!counts.containsKey(count)) { counts.put(count, 0); }
                    counts.put(count, counts.get(count) + 1);
                    maxCount = Math.max(count, maxCount);
                }
            }
        }
        
        for(int i = 0; i <= maxCount; i++) { 
            if(!counts.containsKey(i)) { 
                System.out.println(i + ": 0");
            } else { 
                System.out.println(i + ": " + counts.get(i));
            }
        }
    }
    
    public static void examineSGDTables(Genome g) throws SQLException { 
        java.sql.Connection cxn = g.getUcscConnection();
        
        int sgdGeneCount = getCount(cxn, "select count(*) from sgdGene");
        int sgdGenePlus = getCount(cxn, "select count(*) from sgdGene where strand='+'");
        int sgdOtherCount = getCount(cxn, "select count(*) from sgdOther");
        int sgdOtherPlus = getCount(cxn, "select count(*) from sgdOther where strand='+'");
        
        System.out.println("# sgdGene entries: " + sgdGeneCount);
        System.out.println("\t# Positive Strand: " + sgdGenePlus);
        System.out.println("# sgdOther entries: " + sgdOtherCount);
        System.out.println("\t# Positive Strand: " + sgdOtherPlus);
        
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select distinct(type) from sgdOther");
        System.out.println("SGDOther Types:");
        while(rs.next()) { 
            String t = rs.getString(1);
            System.out.println("\t" + t);
        }
        rs.close();
        s.close();
        
        DatabaseFactory.freeConnection(cxn);
    }
    
    public static int getCount(java.sql.Connection cxn, String query) throws SQLException {
        int count = -1;
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery(query);
        if(rs.next()) { count = rs.getInt(1); }
        rs.close();
        s.close();
        return count;
    }
}
