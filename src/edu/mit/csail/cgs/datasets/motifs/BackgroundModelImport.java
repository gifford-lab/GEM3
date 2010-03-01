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

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.Points2RegionsConverter;
import edu.mit.csail.cgs.utils.io.motifs.BackgroundModelIO;
import edu.mit.csail.cgs.utils.io.parsing.PWMParser;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.species.Genome;
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
   * java edu.mit.csail.cgs.datasets.motifs.BackgroundModelImport --genome "Mus musculus;mm8" --bgname "whole genome" --bgtype MARKOV --bgfile foo.back
   */
  public static void main(String args[]) {
  	try {
  		Genome gen = null;
  		String bgName = null;
  		String bgType = null;
  		String bgFilename = null;
  		
  		args = new String[] {"--species", "Mus musculus;mm8", "--bgname", "whole genome", "--bgtype", "MARKOV", "--bgfile", "mm8.back"};

  		gen = Args.parseGenome(args).cdr();
  		bgName = Args.parseString(args, "bgname", null);
  		bgType = Args.parseString(args, "bgtype", null);
  		bgFilename = Args.parseString(args, "bgfile", null);

  		if (gen == null) {
  			logger.fatal("Must specify a genome in --species"); System.exit(1);
  		}
  		if (bgFilename == null) {
  			logger.fatal("Must supply a --bgfile"); System.exit(1);
  		} 
  		if (bgName == null) {
  			logger.fatal("Must supply a --bgname"); System.exit(1);
  		}
  		if (bgType == null) {
  			logger.fatal("Must supply a --bgtype"); System.exit(1);
  		}

  		if (bgType.toUpperCase().equals("MARKOV")) {
  			MarkovBackgroundModel mbg = BackgroundModelIO.parseMarkovBackgroundModel(bgName, bgFilename, gen);
  			BackgroundModelImport.insertMarkovBackgroundModel(mbg);
  		}
  	}
  	catch (NotFoundException nfex) {
  		logger.fatal(nfex);
  	}
  	catch (IOException ioex) {
  		logger.fatal(ioex);
  	}
  	catch (ParseException pex) {
  		logger.fatal(pex);  		
  	}
  	catch (SQLException sqlex) {
  		logger.fatal(sqlex);
  	}
  	System.exit(0);
  }
      

  private static void insertModelColumnsNoCounts(BackgroundModel model, int bggmID, Connection cxn) throws SQLException {
  	boolean isFreqModel;
  	if (model instanceof FrequencyBackgroundModel) {
  		isFreqModel = true;
  	}
  	else if (model instanceof MarkovBackgroundModel) {
  		isFreqModel = false;
  	}
  	else {
  		throw new IllegalArgumentException("Model must be a MarkovBackgroundModel or a FrequencyBackgroundModel. Use insertModelColumns() for a CountsBackgroundModel.");
  	}
  	
  	PreparedStatement insertCol = null;
  	for (int i = 1; i <= model.getMaxKmerLen(); i++) {
  		for (String kmer : model.getKmers(i)) {
  			double prob;
  			if (isFreqModel) {
  				prob = ((FrequencyBackgroundModel)model).getFrequency(kmer);  				
  			}
  			else {
  				prob = model.getMarkovProb(kmer);
  			}

				insertCol = cxn.prepareStatement("insert into background_model_cols(bggm_id,kmer,probability) values(?,?,?)");
  			insertCol.setInt(1, bggmID);
  			insertCol.setString(2, kmer);
  			insertCol.setDouble(3, prob);
  			insertCol.execute();
  		}
  	}
  	if (insertCol != null) {
  		insertCol.close();
  	}
  }
  
  
  public static int insertBackgroundModelTable(BackgroundModel model, String modelType, Connection cxn) throws SQLException {    
    int modelID = -1;
    PreparedStatement modelExists = cxn.prepareStatement("select id from background_model where name = ? and max_kmer_len = ? and model_type = " + modelType);
    modelExists.setString(1, model.name);
    modelExists.setInt(2, model.getMaxKmerLen());
    ResultSet rs = modelExists.executeQuery();
    
    if (rs.next()) {
      modelID = rs.getInt(1);
      rs.close();
      modelExists.close();
    }
    else {
      rs.close();
      modelExists.close();
      //Note: reuse the weightmatrix_id sequence
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
    return modelID;
  }
  
  
  public static int insertBackgroundGenomeMapTable(BackgroundModel model, int modelID, Connection cxn) throws SQLException {
  	int bggmID = -1;
  	PreparedStatement bgGenomeMapExists = cxn.prepareStatement("select id from background_genome_map where bg_model_id = " + modelID + " and genome_id = " + model.getGenome().getDBID());
    ResultSet rs = bgGenomeMapExists.executeQuery();
    if (rs.next()) {
    	rs.close();
    	bgGenomeMapExists.close();
      throw new IllegalArgumentException("Model already exists. Set a different name or use updateMarkovModel() to update the existing model.");
    }
    else {
      rs.close();
      bgGenomeMapExists.close();
      //Note: reuse the weightmatrix_id sequence
      PreparedStatement insertBGGM = 
        cxn.prepareStatement("insert into background_genome_map(id, genome_id, bg_model_id) values (weightmatrix_id.nextval," + model.getGenome().getDBID() + "," + modelID + ")");
      insertBGGM.execute();
      insertBGGM.close();
      PreparedStatement getBGGMID = cxn.prepareStatement("select weightmatrix_id.currval from dual");
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
    return bggmID;
  }
  
  
  public static int insertMarkovBackgroundModel(MarkovBackgroundModel model) throws SQLException {
  	java.sql.Connection cxn = DatabaseFactory.getConnection("annotations");
  	cxn.setAutoCommit(false);
  	int bggmID = -1;
  	try {
  	//make sure the model has a name and genome
  	if ((model.getName() == null) || (model.getName().isEmpty()) || (model.getGenome() == null)) {
  		throw new IllegalArgumentException("Model must have a name and genome specified to be imported to database.");
  	}
  	
    /**
     * Check whether there is already an entry for a model with this name, kmerlen, and type. If so, resuse the model ID, 
     * otherwise create one.
     */
  	int modelID = BackgroundModelImport.insertBackgroundModelTable(model, "'MARKOV'", cxn);
    
    /**
     * Check whether there is already an entry for this model's genome
     */
    bggmID = BackgroundModelImport.insertBackgroundGenomeMapTable(model, modelID, cxn);

    //now insert all the columns for this genome
    BackgroundModelImport.insertModelColumnsNoCounts(model, bggmID, cxn);
    
    cxn.commit();
  	}
  	catch (Exception ex) {
  		cxn.rollback();
  		logger.fatal(ex);
  	}
    finally {
    	if (cxn != null) {
    		DatabaseFactory.freeConnection(cxn);
    	}    	                  
    }
  	
    
    return bggmID;
  }
}
