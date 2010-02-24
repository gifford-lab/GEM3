/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;
import java.text.ParseException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.PWMParser;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author rca
 *
 */
public class BackgroundModelImport {

    
  /** Imports a weight matrix from a TAMO formatted file.
   *  Only reads the first WM in the file. 
   *
   * Usage:
   * java edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport --species "Saccharomyces cerevisiae" --wmname HSF1 --wmversion MacIsaac06 --wmtype TAMO --wmfile v1.HSF1.wm
   */
//  public static void main(String args[]) {
//    String species = null;
//    String wmname = null, wmversion = null, wmtype = null;
//    String wmfile = null;
//    for (int i = 0; i < args.length; i++) {
//      if (args[i].equals("--species")) { 
//        species = args[++i];
//        if (species.indexOf(';') != -1) {
//          String[] pieces = species.split(";");
//          species = pieces[0];
//        }
//      }
//      if (args[i].equals("--wmname")) {
//        wmname = args[++i];
//        if (wmname.indexOf(';') != -1) {
//          String[] pieces = wmname.split(";");
//          wmname = pieces[0];
//          wmversion = pieces[1];
//          if (pieces.length >= 3) {
//            wmtype = pieces[2];
//          }
//        }
//      }
//      if (args[i].equals("--wmversion")) {
//        wmversion = args[++i];
//      }
//      if (args[i].equals("--wmtype")) {
//        wmtype = args[++i];
//      }
//      if (args[i].equals("--") ||
//          args[i].equals("--wmfile")) {
//        wmfile = args[++i];
//      }
//    }
//
//    //        if (species == null) {
//    //            System.err.println("Must supply a --species"); System.exit(1);
//    //        }
//    if (wmfile == null) {
//      System.err.println("Must supply a --wmfile"); System.exit(1);
//    } 
//    try {
//      if(wmname==null) { 
//        insertMultiWMFromFile(species,wmtype,wmfile, wmversion);
//      } 
//      else { 
//        if (wmversion == null) {
//          System.err.println("Must supply a --wmversion"); System.exit(1);
//        }
//        insertWMFromFile(species,wmname,wmversion,wmtype,wmfile);
//      }
//    } 
//    catch (SQLException ex) {
//      ex.printStackTrace();
//    } 
//    catch (NotFoundException ex) {
//      ex.printStackTrace();
//      System.err.println("Must supply a valid species and genome");
//    } 
//    catch (UnknownRoleException ex ){
//      ex.printStackTrace();
//      System.err.println("Couldn't connect to role annotations");
//    } 
//    catch (FileNotFoundException ex) {
//      ex.printStackTrace();
//      System.err.println("Couldn't find the input file");
//    } 
//    catch (ParseException ex) {
//      ex.printStackTrace();
//    } 
//    catch (IOException ex) {
//      ex.printStackTrace();
//    }
//  }
      

  public static int insertMarkovModelIntoDB(MarkovBackgroundModel model) throws SQLException, NotFoundException {
//FIXME
//    java.sql.Connection cxn = DatabaseFactory.getConnection("annotations");
//    int modelID = -1;
//    PreparedStatement exists = cxn.prepareStatement("select id from background_model where name = ? and model_order = ?");
//    exists.setString(1, model.name);
//    exists.setInt(2, model.getMaxKmerLen());
//    ResultSet rs = exists.executeQuery();
//    if (!rs.next()) {
//      rs.close();
//      PreparedStatement insertwm = cxn.prepareStatement("insert into background_model(id,name,model_order) values (weightmatrix_id.nextval,?,?)");
//      insertwm.setString(1, model.name);
//      insertwm.setInt(2, model.getMaxKmerLen());
//      insertwm.execute();
//      insertwm.close();
//      PreparedStatement getModelId = cxn.prepareStatement("select weightmatrix_id.currval from dual");
//      rs = getModelId.executeQuery();
//      if (rs.next()) {
//        modelID = rs.getInt(1);
//      }
//      else {
//        System.err.println("No Background Model ID");
//        System.exit(1);
//      }
//      rs.close();
//      getModelId.close();
//    }
//    else {
//      modelID = rs.getInt(1);
//      PreparedStatement deleteold = cxn.prepareStatement("delete from weightmatrixcols where weightmatrix = ?");
//      deleteold.setInt(1, modelID);
//      deleteold.execute();
//      deleteold.close();
//    }
//    rs.close();
//    exists.close();
//
//
//    PreparedStatement insertcol = cxn
//        .prepareStatement("insert into weightmatrixcols(weightmatrix,position,letter,weight) values(?,?,?,?)");
//    insertcol.setInt(1, wmid);
//    for (int col = 0; col < matrix.length(); col++) {
//      // System.err.println(String.format("  Column %d %f %f %f %f", col,
//      // matrix.matrix[col]['A'],
//      // matrix.matrix[col]['C'],
//      // matrix.matrix[col]['T'],
//      // matrix.matrix[col]['G']));
//      insertcol.setInt(2, col);
//
//      insertcol.setString(3, "A");
//      insertcol.setFloat(4, matrix.matrix[col]['A']);
//      insertcol.execute();
//
//      insertcol.setString(3, "C");
//      insertcol.setFloat(4, matrix.matrix[col]['C']);
//      insertcol.execute();
//
//      insertcol.setString(3, "T");
//      insertcol.setFloat(4, matrix.matrix[col]['T']);
//      insertcol.execute();
//
//      insertcol.setString(3, "G");
//      insertcol.setFloat(4, matrix.matrix[col]['G']);
//      insertcol.execute();
//    }
//    insertcol.close();
//    DatabaseFactory.freeConnection(cxn);
//    return wmid;
    return 0;
  }
}
