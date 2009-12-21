/**
 * 
 */
package edu.mit.csail.cgs.tools.chipseq;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;

/**
 * @author rca
 * This class was written to copy the chipseqhits table to a new table which uses
 * double-precision values for the weight, in order to fix floating point 
 * arithmetic problems that occur when a sql query sums the weights of a large
 * number of hits
 *
 */
public class ChipSeqHitCopier {
  
  public static PreparedStatement prepareCopyHitsByAlign(java.sql.Connection cxn) throws SQLException {
    String insert = "insert into chipseqhits_new (read, expt, alignment, chromosome, startpos, stoppos, strand, weight) "
      + "select read, expt, alignment, chromosome, startpos, stoppos, strand, TO_BINARY_DOUBLE(weight) from chipseqhits where alignment = ?";
    return cxn.prepareStatement(insert);    
  }

  public static PreparedStatement prepareGetAlignments(java.sql.Connection cxn) throws SQLException {
    String query ="select id from chipseqalignments order by id";
    return cxn.prepareStatement(query);
  }

  
  public static void main(String args[]) {
    java.sql.Connection cxn = null;
    PreparedStatement getAlignments = null;
    PreparedStatement copyHits = null;
    
    try {      
      cxn = DatabaseFactory.getConnection("chipseq", "schema goes here", "password goes here");

      //prepare statements
      getAlignments = ChipSeqHitCopier.prepareGetAlignments(cxn);
      getAlignments.setFetchSize(1000);
    
      copyHits = ChipSeqHitCopier.prepareCopyHitsByAlign(cxn);
      
      //Get the IDS of all the alignments
      ArrayList<Integer> alignIDs = new ArrayList<Integer>();
      System.out.println("Getting all Alignment IDs");
      ResultSet alignRS = getAlignments.executeQuery();
      while(alignRS.next()) {
        alignIDs.add(alignRS.getInt(1));
      }
      alignRS.close();
      getAlignments.close(); getAlignments = null;
      System.out.println("Done getting " + alignIDs.size() + " Alignment IDs. Min: " + alignIDs.get(0) + " Max: " + alignIDs.get(alignIDs.size() - 1));
      
      //Turn off auto-commit so this can be done with transactions
      cxn.setAutoCommit(false);
      
      
      //Copy over the hits
      System.out.println("Starting to copy hits");
      for (int i = 0; i < alignIDs.size(); i++) {
        int currentID = alignIDs.get(i);
        copyHits.setInt(1, currentID);
        copyHits.execute();
        cxn.commit();
        System.out.println(i + ": " + currentID);
      }    
      System.out.println("\n\nDone");
    }
    catch (SQLException ex) {
      // TODO Auto-generated catch block
      ex.printStackTrace();
    }
    finally {
      try {
        if (getAlignments != null) {
          getAlignments.close();
        }
        if (copyHits != null) {
          copyHits.close();
        }
        if (cxn != null) {
          cxn.close();
        }
      }
      catch (SQLException ex) {
        // TODO Auto-generated catch block
        ex.printStackTrace();
      }
    }
  }
}
