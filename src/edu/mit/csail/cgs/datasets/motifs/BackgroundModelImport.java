/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;
import java.text.ParseException;

import org.apache.log4j.Logger;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.Points2RegionsConverter;
import edu.mit.csail.cgs.utils.io.parsing.PWMParser;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author rca
 *
 */
public class BackgroundModelImport {

  private static Logger logger = Logger.getLogger(BackgroundModelImport.class);
    
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
      

  private static void insertModelColumns(BackgroundModel model, int bggmID, Connection cxn) throws SQLException {
//    if ((model instanceof MarkovBackgroundModel) || (model instanceof FrequencyBackgroundModel)) {
//      PreparedStatement insertCol = null;
//      for (int i = 1; i <= model.getMaxKmerLen(); i++) {
//        for (String kmer : model.modelProbs[i].keySet()) {
//          String prevBases = kmer.substring(0, kmer.length() - 1);
//          String currBase = kmer.substring(kmer.length() - 1, kmer.length());
//          double prob = model.getModelProb(kmer);
//          int count = model.getModelCount(kmer);
//
//          if (count >= 0) {
//            insertCol = cxn.prepareStatement("insert into background_model_cols(bggm_id,prev_bases,curr_base,probability,count) values(?,?,?,?,?)");
//            insertCol.setInt(5, count);
//          }
//          else {
//            insertCol = cxn.prepareStatement("insert into background_model_cols(bggm_id,prev_bases,curr_base,probability) values(?,?,?,?)");
//          }
//          insertCol.setInt(1, bggmID);
//          insertCol.setString(2, prevBases);
//          insertCol.setString(3, currBase);
//          insertCol.setDouble(4, prob);
//          insertCol.execute();
//        }
//      }
//      if (insertCol != null) {
//        insertCol.close();
//      }
//    }
//    else {
//      throw new IllegalArgumentException("Background Model must either be a MarkovBackgroundModel or a FrequencyBackgroundModel.");
//    }
  }
  
  
  public static int insertMarkovModel(MarkovBackgroundModel model) throws SQLException, NotFoundException {
//FIXME
    java.sql.Connection cxn = DatabaseFactory.getConnection("annotations");
    int modelID = -1;
    int bggmID = -1;
    PreparedStatement modelExists = cxn.prepareStatement("select id from background_model where name = ? and max_kmer_len = ? and model_type = 'MARKOV'");
    modelExists.setString(1, model.name);
    modelExists.setInt(2, model.getMaxKmerLen());
    ResultSet rs = modelExists.executeQuery();
    
    /**
     * Check whether there is already an entry for a model with this name, kmerlen, and type. If so, resuse the model ID, 
     * otherwise create one.
     */
    if (rs.next()) {
      modelID = rs.getInt(1);
      rs.close();
    }
    else {
      rs.close();
      PreparedStatement insertBG = cxn.prepareStatement("insert into background_model(id,name,max_kmer_len,model_type) values (weightmatrix_id.nextval,?,?,'MARKOV')");
      insertBG.setString(1, model.name);
      insertBG.setInt(2, model.getMaxKmerLen());
      insertBG.execute();
      insertBG.close();
      PreparedStatement getModelID = cxn.prepareStatement("select weightmatrix_id.currval from dual");
      rs = getModelID.executeQuery();
      if (rs.next()) {
        modelID = rs.getInt(1);
      }
      else {
        logger.fatal("No Background Model ID");
        System.exit(1);
      }
      rs.close();
      getModelID.close();
    }
    
    /**
     * Check whether there is already an entry for this model's genome
     */
    PreparedStatement bgGenomeMapExists = cxn.prepareStatement("select id from background_genome_map where bg_model_id = " + modelID + " and genome_id = " + model.gen.getDBID());
    rs = bgGenomeMapExists.executeQuery();
    if (rs.next()) {
      throw new IllegalArgumentException("Model already exists. Select a different name or use updateMarkovModel() to update the existing model.");
    }
    else {
      rs.close();
      PreparedStatement insertBGGM = 
        cxn.prepareStatement("insert into background_genome_map(id, genome_id, bg_model_id) values (bggm_id.nextval," + model.gen.getDBID() + "," + modelID + ")");
      insertBGGM.execute();
      insertBGGM.close();
      PreparedStatement getBGGMID = cxn.prepareStatement("select bggm_id.currval from dual");
      rs = getBGGMID.executeQuery();
      if (rs.next()) {
        bggmID = rs.getInt(1);
      }
      else {
        logger.fatal("No Background Model Genome Map ID");
        System.exit(1);
      }
      rs.close();
      getBGGMID.close();
    }

    BackgroundModelImport.insertModelColumns(model, bggmID, cxn);
    return bggmID;
  }
  
  
  
}
