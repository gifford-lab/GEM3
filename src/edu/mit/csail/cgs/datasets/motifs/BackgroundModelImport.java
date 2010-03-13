/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

import java.io.*;
import java.sql.*;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.io.motifs.BackgroundModelIO;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author rca
 * Code for database interaction involving background models
 */
public class BackgroundModelImport {

  private static Logger logger = Logger.getLogger(BackgroundModelImport.class);
  
  public static final String FREQUENCY_TYPE_STRING = "FREQUENCY";
  public static final String MARKOV_TYPE_STRING = "MARKOV";
  
  //TODO Make all the SQLs used by this class be constants
  private static final String SQL_GET_MODEL_ID = "select id from background_model where name = ? and max_kmer_len = ? and model_type = ?";
  
  private static final String SQL_GET_MAP_ID = "select id from background_genome_map where bg_model_id = ? and genome_id = ?";
  
  private static final String SQL_GET_ALL_MODELS = "select id, name, max_kmer_len, model_type from background_model";
  
  private static final String SQL_GET_GENOMES = "select genome_id from background_genome_map where bg_model_id = ?";
    
  private static final String SQL_GET_METADATA_BY_MODEL_ID = "select name, max_kmer_len, model_type from background_model where id = ?";
  
  private static final String SQL_GET_METADATA_CORE = 
    "select map.id, map.genome_id, map.bg_model_id, bm.name, bm.max_kmer_len, bm.model_type, map.has_counts"
    + " from background_model bm, background_genome_map map"
    + " where bm.id = map.bg_model_id";
  private static final String SQL_GET_METADATA_ORDER_BY = " order by bm.name, bm.max_kmer_len, map.genome_id";
  
  
  private static final String SQL_GET_METADATA_MAP_ID = " and map.id = ?";
  private static final String SQL_GET_METADATA_GENOME_ID = " and map.genome_id = ?";
  private static final String SQL_GET_METADATA_MODEL_ID = " and map.bg_model_id = ?";
  private static final String SQL_GET_METADATA_NAME = " and bm.name = ?";
  private static final String SQL_GET_METADATA_KMER_LEN = " and bm.max_kmer_len = ?";
  private static final String SQL_GET_METADATA_MODEL_TYPE = " and bm.model_type = ?";

  private static final int SQL_GET_METADATA_CORE_MAP_ID_INDEX = 1;
  private static final int SQL_GET_METADATA_CORE_GENOME_ID_INDEX = 2;
  private static final int SQL_GET_METADATA_CORE_MODEL_ID_INDEX = 3;
  private static final int SQL_GET_METADATA_CORE_NAME_INDEX = 4;
  private static final int SQL_GET_METADATA_CORE_KMER_LEN_INDEX = 5;
  private static final int SQL_GET_METADATA_CORE_MODEL_TYPE_INDEX = 6;
  private static final int SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX = 7;
  
  
  
  
  private static final String SQL_GET_MODEL_BY_ID = "select map_id, kmer, probability, count from background_model_cols where map_id = ? order by map_id, length(kmer), kmer";
  private static final int SQL_GET_MODEL_BY_ID_MAP_ID_INDEX = 1;
  private static final int SQL_GET_MODEL_BY_ID_KMER_INDEX = 2;
  private static final int SQL_GET_MODEL_BY_ID_PROB_INDEX = 3;
  private static final int SQL_GET_MODEL_BY_ID_COUNT_INDEX = 4;
  
  private static final String SQL_GET_MODEL_CORE = 
    "select map.id, map.genome_id, bm.name, bm.max_kmer_len, bm.id, bgmc.kmer, bgmc.probability, bgmc.count"
    + " from background_model bm, background_genome_map map, background_model_cols bgmc"
    + " where bm.id = map.bg_model_id and bgmc.map_id = map.id and bm.model_type = ?";
  private static final String SQL_GET_MODEL_ORDER_BY = " order by bm.name, bm.max_kmer_len, map.genome_id, length(bgmc.kmer), bgmc.kmer";
  private static final int SQL_GET_MODEL_CORE_MAP_ID_INDEX = 1;
  private static final int SQL_GET_MODEL_CORE_GENOME_ID_INDEX = 2;
  private static final int SQL_GET_MODEL_CORE_NAME_INDEX = 3;
  private static final int SQL_GET_MODEL_CORE_KMERLEN_INDEX = 4;
  private static final int SQL_GET_MODEL_CORE_MODEL_ID_INDEX = 5;
  private static final int SQL_GET_MODEL_CORE_KMER_INDEX = 6;
  private static final int SQL_GET_MODEL_CORE_PROB_INDEX = 7;
  private static final int SQL_GET_MODEL_CORE_COUNT_INDEX = 8;
  
  private static final String SQL_GET_MODEL_MAP_ID = " and map.map_id = ?";

  private static final String SQL_GET_MODEL_GENOME_ID = " and map.genome_id = ?";

  private static final String SQL_GET_MODEL_MODEL_ID = " and bm.id = ?";

  private static final String SQL_GET_MODEL_NAME = " and bm.name = ?";
  
  private static final String SQL_GET_MODEL_KMER_LEN = " and bm.max_kmer_len = ?";  
  
    
  /**
   * @see getBackgroundModelID(String, int, String, Connection)
   * Creates a database connection for the query
   */
  public static Integer getBackgroundModelID(String name, int kmerLen, String modelType) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModelID(name, kmerLen, modelType, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  /**
   * Look for a model with this name, maxKmerLen, and type and return its ID. 
   * 
   * @param name the name of the model
   * @param kmerLen the length of the longest kmer in the model
   * @param dbModelType the type of the model (FREQUENCY or MARKOV)
   * @param cxn an open database connection
   * @return the ID of the model, or null if there's no match
   * @throws SQLException
   */
  public static Integer getBackgroundModelID(String name, int kmerLen, String modelType, Connection cxn) throws SQLException {
    PreparedStatement getModelID = null;
    ResultSet rs = null;    
    try {
      getModelID = cxn.prepareStatement(SQL_GET_MODEL_ID);
      getModelID.setString(1, name);
      getModelID.setInt(2, kmerLen);
      getModelID.setString(3, modelType);
      rs = getModelID.executeQuery();

      if (rs.next()) {
        Integer modelID = rs.getInt(1);
        return modelID;
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModelID != null) {
        getModelID.close();
      }
    }
  }
  
  
  /**
   * @see getBackgroundGenomeMapID(int, int, Connection)
   * Creates a database connection for the query
   */
  public static Integer getBackgroundGenomeMapID(int bgModelID, int genomeID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundGenomeMapID(bgModelID, genomeID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  /**
   * Look for an entry in the background model genome map for the specified
   * background model and genome, and return its ID
   * @param bgModelID the id of the model to look up
   * @param genomeID a genome for which there may be an instance of the model
   * @param cxn an open database connection
   * @return the background genome map ID of the model, or null if there's no match
   * @throws SQLException
   */
  public static Integer getBackgroundGenomeMapID(int bgModelID, int genomeID, Connection cxn) throws SQLException {
    PreparedStatement getBGGenomeMapID = null;
    ResultSet rs = null;
    try {
      getBGGenomeMapID = cxn.prepareStatement(SQL_GET_MAP_ID);
      getBGGenomeMapID.setInt(1, bgModelID);
      getBGGenomeMapID.setInt(2, genomeID);
      rs = getBGGenomeMapID.executeQuery();
      if (rs.next()) {
        return rs.getInt(1);
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getBGGenomeMapID != null) {
        getBGGenomeMapID.close();
      }      
    }
  }
  
  
  /**************************************************************************
   * Methods for looking up which background models are in the database
   **************************************************************************/
  
  
  /**
   * 
   * @param modelID
   * @return
   * @throws SQLException
   */  
  public static BackgroundModelMetadata getBackgroundModelByModelID(int modelID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModelByModelID(modelID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModelByModelID(int modelID, Connection cxn) throws SQLException {
    PreparedStatement getModel = null;
    ResultSet rs = null;
    try {
      getModel = cxn.prepareStatement(SQL_GET_METADATA_BY_MODEL_ID);
      getModel.setInt(1, modelID);
      rs = getModel.executeQuery();
      if (rs.next()) {
        return new BackgroundModelMetadata(modelID, rs.getString(1), rs.getInt(2), rs.getString(3));
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModel != null) {
        getModel.close();
      }      
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModel(int modelID, int genomeID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModel(modelID, genomeID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModel(int modelID, int genomeID, Connection cxn) throws SQLException {
    PreparedStatement getModel = null;
    ResultSet rs = null;
    try {
      StringBuffer sql = new StringBuffer(SQL_GET_METADATA_CORE);
      sql.append(SQL_GET_METADATA_MODEL_ID);
      sql.append(SQL_GET_METADATA_GENOME_ID);
      sql.append(SQL_GET_METADATA_ORDER_BY);
      
      getModel = cxn.prepareStatement(sql.toString());
      getModel.setInt(1, modelID);
      getModel.setInt(2, genomeID);
      rs = getModel.executeQuery();
      if (rs.next()) {
        return new BackgroundModelMetadata(rs.getInt(SQL_GET_METADATA_CORE_MAP_ID_INDEX), genomeID, modelID, 
            rs.getString(SQL_GET_METADATA_CORE_NAME_INDEX), rs.getInt(SQL_GET_METADATA_CORE_KMER_LEN_INDEX), 
            rs.getString(SQL_GET_METADATA_CORE_MODEL_TYPE_INDEX), rs.getBoolean(SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX));
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModel != null) {
        getModel.close();
      }      
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModel(String name, int maxKmerLen, String modelType, int genomeID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModel(name, maxKmerLen, modelType, genomeID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModel(String name, int maxKmerLen, String modelType, int genomeID, Connection cxn) throws SQLException {
    PreparedStatement getModel = null;
    ResultSet rs = null;
    try {
      StringBuffer sql = new StringBuffer(SQL_GET_METADATA_CORE);
      sql.append(SQL_GET_METADATA_NAME);
      sql.append(SQL_GET_METADATA_KMER_LEN);
      sql.append(SQL_GET_METADATA_MODEL_TYPE);
      sql.append(SQL_GET_METADATA_GENOME_ID);
      sql.append(SQL_GET_METADATA_ORDER_BY);

      getModel = cxn.prepareStatement(sql.toString());
      getModel.setString(1, name);
      getModel.setInt(2, maxKmerLen);
      getModel.setString(3, modelType);
      getModel.setInt(4, genomeID);
      rs = getModel.executeQuery();
      if (rs.next()) {
        return new BackgroundModelMetadata(rs.getInt(SQL_GET_METADATA_CORE_MAP_ID_INDEX), genomeID, rs.getInt(SQL_GET_METADATA_CORE_MODEL_ID_INDEX), 
            name, maxKmerLen, modelType, rs.getBoolean(SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX));
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModel != null) {
        getModel.close();
      }      
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModelByMapID(int mapID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModelByMapID(mapID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static BackgroundModelMetadata getBackgroundModelByMapID(int mapID, Connection cxn) throws SQLException {
    PreparedStatement getModel = null;
    ResultSet rs = null;
    try {
      StringBuffer sql = new StringBuffer(SQL_GET_METADATA_CORE);
      sql.append(SQL_GET_METADATA_MAP_ID);
      sql.append(SQL_GET_METADATA_ORDER_BY);

      getModel = cxn.prepareStatement(sql.toString());
      getModel.setInt(1, mapID);
      rs = getModel.executeQuery();
      if (rs.next()) {
        return new BackgroundModelMetadata(mapID, rs.getInt(SQL_GET_METADATA_CORE_GENOME_ID_INDEX), rs.getInt(SQL_GET_METADATA_CORE_MODEL_ID_INDEX),
            rs.getString(SQL_GET_METADATA_CORE_NAME_INDEX), rs.getInt(SQL_GET_METADATA_CORE_KMER_LEN_INDEX), 
            rs.getString(SQL_GET_METADATA_CORE_MODEL_TYPE_INDEX), rs.getBoolean(SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX));
      }
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModel != null) {
        getModel.close();
      }      
    }
  }
  
  
  public static List<BackgroundModelMetadata> getAllBackgroundModels() throws SQLException {
    return BackgroundModelImport.getAllBackgroundModels(false);
  }
  
  
  public static List<BackgroundModelMetadata> getAllBackgroundModels(boolean ignoreGenome) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getAllBackgroundModels(ignoreGenome, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  

  public static List<BackgroundModelMetadata> getAllBackgroundModels(Connection cxn) throws SQLException {
    return BackgroundModelImport.getAllBackgroundModels(false, cxn);
  }

  
  public static List<BackgroundModelMetadata> getAllBackgroundModels(boolean ignoreGenome, Connection cxn) throws SQLException {
    PreparedStatement getAllModels = null;
    ResultSet rs = null;
    try {
      if (ignoreGenome) {
        getAllModels = cxn.prepareStatement(SQL_GET_ALL_MODELS);
        rs = getAllModels.executeQuery();
        List<BackgroundModelMetadata> results = new ArrayList<BackgroundModelMetadata>();
        while (rs.next()) {        
          results.add(new BackgroundModelMetadata(rs.getInt(1), rs.getString(2), rs.getInt(3), rs.getString(4)));
        }
        return results;
      }
      else {
        StringBuffer sql = new StringBuffer(SQL_GET_METADATA_CORE);
        sql.append(SQL_GET_METADATA_ORDER_BY);

        getAllModels = cxn.prepareStatement(sql.toString());
        rs = getAllModels.executeQuery();
        List<BackgroundModelMetadata> results = new ArrayList<BackgroundModelMetadata>();
        while (rs.next()) {    
          results.add(new BackgroundModelMetadata(rs.getInt(SQL_GET_METADATA_CORE_MAP_ID_INDEX), rs.getInt(SQL_GET_METADATA_CORE_GENOME_ID_INDEX), 
              rs.getInt(SQL_GET_METADATA_CORE_MODEL_ID_INDEX), rs.getString(SQL_GET_METADATA_CORE_NAME_INDEX), 
              rs.getInt(SQL_GET_METADATA_CORE_KMER_LEN_INDEX), rs.getString(SQL_GET_METADATA_CORE_MODEL_TYPE_INDEX), 
              rs.getBoolean(SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX)));
        }
        return results;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getAllModels != null) {
        getAllModels.close();
      }      
    }
  }
  
  
  public static List<BackgroundModelMetadata> getBackgroundModelsForGenome(int genomeID) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getBackgroundModelsForGenome(genomeID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static List<BackgroundModelMetadata> getBackgroundModelsForGenome(int genomeID, Connection cxn) throws SQLException {
    PreparedStatement getGenomeModels = null;
    ResultSet rs = null;    
    try {
      StringBuffer sql = new StringBuffer(SQL_GET_METADATA_CORE);
      sql.append(SQL_GET_METADATA_GENOME_ID);
      sql.append(SQL_GET_METADATA_ORDER_BY);

      getGenomeModels = cxn.prepareStatement(sql.toString());
      getGenomeModels.setInt(1, genomeID);
      rs = getGenomeModels.executeQuery();
      List<BackgroundModelMetadata> results = new ArrayList<BackgroundModelMetadata>();
      while (rs.next()) {
        results.add(new BackgroundModelMetadata(rs.getInt(SQL_GET_METADATA_CORE_MAP_ID_INDEX), genomeID, 
            rs.getInt(SQL_GET_METADATA_CORE_MODEL_ID_INDEX), rs.getString(SQL_GET_METADATA_CORE_NAME_INDEX), 
            rs.getInt(SQL_GET_METADATA_CORE_KMER_LEN_INDEX), rs.getString(SQL_GET_METADATA_CORE_MODEL_TYPE_INDEX), 
            rs.getBoolean(SQL_GET_METADATA_CORE_HAS_COUNTS_INDEX)));
      }
      return results;
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getGenomeModels != null) {
        getGenomeModels.close();
      }      
    }
  }
  
  
  public static List<Integer> getGenomesForBackgroundModel(int modelID) throws SQLException{
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getGenomesForBackgroundModel(modelID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static List<Integer> getGenomesForBackgroundModel(int modelID, Connection cxn) throws SQLException {
    PreparedStatement getModels = null;
    ResultSet rs = null;
    try {
      getModels = cxn.prepareStatement(SQL_GET_GENOMES);
      getModels.setInt(1, modelID);
      List<Integer> results = new ArrayList<Integer>();
      rs = getModels.executeQuery();
      while (rs.next()) {
        results.add(rs.getInt(1));
      }
      return results;
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModels != null) {
        getModels.close();
      }
    }
  }
  
  
  /**************************************************************************
   * Methods for loading frequency background models
   **************************************************************************/  
  
  /**
   * Given a FrequencyBackgroundModel object and a result set whose rows are 
   * kmer probabilities (and perhaps other data), parse those rows and 
   * set the model probabilities accordingly
   * @param fbm A FrequencyBackgroundModel object to initialize
   * @param rs a ResultSet of kmer probability data, already on the first row
   * for the specified background model
   * @param queryCore the query core used to obtain the result set. This is used
   * to determine which indices to check for the kmer probabilities
   * @return true if there are more rows in the result set for another background
   * model, false otherwise
   * @throws SQLException
   */
  private static boolean initFrequencyModelProbs(FrequencyBackgroundModel fbm, ResultSet rs, String queryCore) throws SQLException {
    int idIndex;
    int kmerIndex;
    int probIndex;
    
    if (queryCore.equals(SQL_GET_MODEL_BY_ID)) {
      idIndex = SQL_GET_MODEL_BY_ID_MAP_ID_INDEX;
      kmerIndex = SQL_GET_MODEL_BY_ID_KMER_INDEX;
      probIndex = SQL_GET_MODEL_BY_ID_PROB_INDEX;
    }
    else if (queryCore.equals(SQL_GET_MODEL_CORE)) {
      idIndex = SQL_GET_MODEL_CORE_MAP_ID_INDEX;
      kmerIndex = SQL_GET_MODEL_CORE_KMER_INDEX;
      probIndex = SQL_GET_MODEL_CORE_PROB_INDEX;
    }
    else {
      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
    }

    int currKmerLen = rs.getString(kmerIndex).length();
    HashMap<String, Double> probs = new HashMap<String, Double>();
    boolean hasNext = true;
    do {
      String kmer = rs.getString(kmerIndex);
      //if this row is a new kmer length then set the probabilities for the
      //previous length
      if (kmer.length() != currKmerLen) {
        if (!probs.isEmpty()) {
          fbm.setKmerFrequencies(probs);
          probs.clear();
        }
        currKmerLen = kmer.length();
      }
      probs.put(kmer, rs.getDouble(probIndex));
      hasNext = rs.next();
    } while (hasNext && (rs.getInt(idIndex) == fbm.getMapID()));

    //set the last batch of kmer probabilities
    fbm.setKmerFrequencies(probs);

    return hasNext;
  }


  /**
   * @see getFrequencyModel(BackgroundModelMetadata md, Connection cxn)
   */
  public static FrequencyBackgroundModel getFrequencyModel(BackgroundModelMetadata md) throws SQLException, NotFoundException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getFrequencyModel(md, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
    
  /**
   * Get a Frequency Background Model described by the specified metadata
   * @param md A complete metadata object specifying a background model. 
   * @param cxn an open db connection to the annotations schema
   * @return A FrequencyBackgroundModel, or null if none exists (and can't be generated)
   * @throws SQLException
   * @throws NotFoundException if the genomeID is invalid, most likely because
   * the metadata object wasn't loaded from the database or was modified afterwards
   */
  public static FrequencyBackgroundModel getFrequencyModel(BackgroundModelMetadata md, Connection cxn) throws SQLException, NotFoundException {
    PreparedStatement getProbs = null;
    ResultSet rs = null;
    try {
      //check that the metadata is complete
      if (!md.hasMapID() || !md.hasGenomeID() || !md.hasModelID() || (md.getName() == null) || !md.hasMaxKmerLen() || (md.getDBModelType() == null)) {
        throw new NullPointerException("Metadata object has one or more null fields");
      }
      
      int mapID = md.getMapID();
      //if the model type is frequency then load it directly
      if (md.getDBModelType().equals(FREQUENCY_TYPE_STRING)) {        
        FrequencyBackgroundModel fbm = new FrequencyBackgroundModel(md);
        
        getProbs = cxn.prepareStatement(SQL_GET_MODEL_BY_ID);
        getProbs.setInt(1, mapID);
        rs = getProbs.executeQuery();
        if (rs.next()) {
          BackgroundModelImport.initFrequencyModelProbs(fbm, rs, SQL_GET_MODEL_BY_ID);
        }
        return fbm;
      }
      //otherwise load the counts if available and convert from the counts model
      else if (md.getDBModelType().equals(MARKOV_TYPE_STRING)) {
        if (BackgroundModelImport.hasCounts(mapID, cxn)) {
          return new FrequencyBackgroundModel(BackgroundModelImport.getCountsModel(md, cxn));
        }
        else {
          return null;
        }
      }
      else {
        //otherwise the model type is invalid
        throw new IllegalArgumentException("Invalid model type: " + md.getDBModelType());
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getProbs != null) {
        getProbs.close();
      }
    }
  }
  
    
  /**
   * Returns a list of Frequency Background Models parsed out of a single result
   * set
   * @param rs the result set containing the data for the background model(s)
   * @param queryCore the SQL query core that was used to generate the specified
   * result set
   * @return A list of FrequencyBackgroundModels
   * @throws SQLException
   */  
  private static List<FrequencyBackgroundModel> createFrequencyModels(ResultSet rs, String queryCore) throws SQLException {
    int mapIDIndex;
    int nameIndex;
    int kmerLenIndex;
    int modelIDIndex;
    int genomeIDIndex;
    
    if (queryCore.equals(SQL_GET_MODEL_CORE)) {
      mapIDIndex = SQL_GET_MODEL_CORE_MAP_ID_INDEX;
      nameIndex = SQL_GET_MODEL_CORE_NAME_INDEX;
      kmerLenIndex = SQL_GET_MODEL_CORE_KMERLEN_INDEX;
      modelIDIndex = SQL_GET_MODEL_CORE_MODEL_ID_INDEX;
      genomeIDIndex = SQL_GET_MODEL_CORE_GENOME_ID_INDEX;
    }
    else {
      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
    }     
    
    try {
      List<FrequencyBackgroundModel> models = new ArrayList<FrequencyBackgroundModel>();
      if (rs.next()) {
        boolean hasNext = true;
        while (hasNext) {
          int mapID = rs.getInt(mapIDIndex);
          FrequencyBackgroundModel mbm = new FrequencyBackgroundModel(rs.getString(nameIndex), Organism.findGenome(rs.getInt(genomeIDIndex)), rs.getInt(kmerLenIndex));
          mbm.setMapID(mapID);
          mbm.setModelID(rs.getInt(modelIDIndex));
          models.add(mbm);
          //call the subroutine to parse all the probabilities from the result set
          hasNext = BackgroundModelImport.initFrequencyModelProbs(mbm, rs, queryCore);
        }
      }
      return models;
    }
    catch (NotFoundException nfex) {
      throw new DatabaseException("Error loading genome for model", nfex);
    }      
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  private static List<FrequencyBackgroundModel> getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getFrequencyModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }

  
  /**
   * Given the parameters that describe one or more background models construct
   * an appropriate SQL query and return the Background Models that are found 
   * @param mapID an ID from the background genome map (this will limit to a 
   * single model)
   * @param genomeID a genome ID for the models
   * @param modelID a model ID from the background_model table
   * @param name a model name
   * @param maxKmerLen a maximum kmer length
   * @param cxn an open DB connection to the annotations schema
   * @return A list of matching Frequency Background Models
   * @throws SQLException
   */
  private static List<FrequencyBackgroundModel> getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn) throws SQLException {
    PreparedStatement getModels = null;
    ResultSet rs = null;
    try {
      //build the appropriate sql query
      StringBuffer sql = new StringBuffer(SQL_GET_MODEL_CORE);
      if (mapID != null) {
        sql.append(SQL_GET_MODEL_MAP_ID);
      }
      if (genomeID != null) {
        sql.append(SQL_GET_MODEL_GENOME_ID);
      }
      if (modelID != null) {
        sql.append(SQL_GET_MODEL_MODEL_ID);
      }
      if (name != null) {
        sql.append(SQL_GET_MODEL_NAME);
      }
      if (maxKmerLen != null) {
        sql.append(SQL_GET_MODEL_KMER_LEN);
      }
      sql.append(SQL_GET_MODEL_ORDER_BY);
      //done building the query      
      
      getModels = cxn.prepareStatement(sql.toString());
      
      //set the sql params
      getModels.setString(1, FREQUENCY_TYPE_STRING);
      int argCount = 2;
      if (mapID != null) {
        getModels.setInt(argCount, mapID);
        argCount++;
      }
      if (genomeID != null) {
        getModels.setInt(argCount, genomeID);
        argCount++;
      }
      if (modelID != null) {
        getModels.setInt(argCount, modelID);
        argCount++;
      }
      if (name != null) {
        getModels.setString(argCount, name);
        argCount++;
      }
      if (maxKmerLen != null) {
        getModels.setInt(argCount, maxKmerLen);
        argCount++;
      }
      //done setting the params
      
      rs = getModels.executeQuery();
      return BackgroundModelImport.createFrequencyModels(rs, SQL_GET_MODEL_CORE);      
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModels != null) {
        getModels.close();
      }
    }
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(String name) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, name, null);
  }

  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(String name, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, name, null,cxn);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(String name, int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, name, maxKmerLen);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(String name, int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, name, maxKmerLen, cxn);
  }

  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(int modelID) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, modelID, null, null);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(int modelID, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, modelID, null, null, cxn);
  }
  

  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(BackgroundModelMetadata md) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getFrequencyModels(md, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }

  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModels(BackgroundModelMetadata md, Connection cxn) throws SQLException {
    Integer mapID = null;
    Integer genomeID = null;
    Integer modelID = null;
    String name = null;
    Integer maxKmerLen = null;
    
    //initialize any fields that are available 
    if (md.hasMapID()) {
      mapID = md.getMapID();
    }
    if (md.hasGenomeID()) {
     genomeID = md.getGenomeID();
    }
    if (md.hasModelID()) {
      modelID = md.getModelID();
    }
    name = md.getName();
    if (md.hasMaxKmerLen()) {
      maxKmerLen = md.getMaxKmerLen();
    }

    return BackgroundModelImport.getFrequencyModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByLength(int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, null, maxKmerLen);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByLength(int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, null, maxKmerLen, cxn);
  }
  
  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, null, null);
  }


  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, null, null, cxn);
  }


  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */  
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID, String name) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, name, null);
  }


  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID, String name, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, name, null, cxn);
  }

  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID, int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, null, maxKmerLen);
  }
  

  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getFrequencyModelsByGenome(int genomeID, int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, genomeID, null, null, maxKmerLen, cxn);
  }

  
  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getAllFrequencyModels() throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, null, null);
  }
  

  /**
   * @see getFrequencyModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<FrequencyBackgroundModel> getAllFrequencyModels(Connection cxn) throws SQLException {
    return BackgroundModelImport.getFrequencyModels(null, null, null, null, null, cxn);
  }

  
  /**************************************************************************
   * Methods for loading Markov background models
   **************************************************************************/

  /**
   * Given a MarkovBackgroundModel object and a result set whose rows are 
   * kmer probabilities (and perhaps other data), parse those rows and 
   * set the model probabilities accordingly
   * @param mbm A MarkovBackgroundModel object to initialize
   * @param rs a ResultSet of kmer probability data, already on the first row
   * for the specified background model
   * @param queryCore the query core used to obtain the result set. This is used
   * to determine which indices to check for the kmer probabilities
   * @return true if there are more rows in the result set for another background
   * model, false otherwise
   * @throws SQLException
   */
  private static boolean initMarkovModelProbs(MarkovBackgroundModel mbm, ResultSet rs, String queryCore) throws SQLException {
    int idIndex;
    int kmerIndex;
    int probIndex;
    
    if (queryCore.equals(SQL_GET_MODEL_BY_ID)) {
      idIndex = SQL_GET_MODEL_BY_ID_MAP_ID_INDEX;
      kmerIndex = SQL_GET_MODEL_BY_ID_KMER_INDEX;
      probIndex = SQL_GET_MODEL_BY_ID_PROB_INDEX;
    }
    else if (queryCore.equals(SQL_GET_MODEL_CORE)) {
      idIndex = SQL_GET_MODEL_CORE_MAP_ID_INDEX;
      kmerIndex = SQL_GET_MODEL_CORE_KMER_INDEX;
      probIndex = SQL_GET_MODEL_CORE_PROB_INDEX;
    }
    else {
      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
    }

    String firstKmer = rs.getString(kmerIndex);
    String currPrevBases = firstKmer.substring(0, firstKmer.length()- 1);
    double[] probs = new double[4];
    boolean hasNext = true;
    do {
      String kmer = rs.getString(kmerIndex);
      String prevBases = kmer.substring(0, kmer.length()-1);
      //if this row is a new value of prev bases then set the probabilities for
      //the old value
      if (!prevBases.equals(currPrevBases)) {
        mbm.setMarkovProb(currPrevBases, probs[0], probs[1], probs[2], probs[3]);
        Arrays.fill(probs, 0.0);
        currPrevBases = prevBases;
      }
      char currBase = kmer.charAt(kmer.length()-1);
      double prob = rs.getDouble(probIndex);
      switch (currBase) {
        case 'A': probs[0] = prob; break;
        case 'C': probs[1] = prob; break;
        case 'G': probs[2] = prob; break;
        case 'T': probs[3] = prob; break;
      }          
      hasNext = rs.next();
    } while (hasNext && (rs.getInt(idIndex) == mbm.getMapID()));
    //set the last batch of kmer probabilities
    mbm.setMarkovProb(currPrevBases, probs[0], probs[1], probs[2], probs[3]);
    return hasNext; 
  }
  
  
  /**
   * @see getMarkovModel(BackgroundModelMetadata md, Connection cxn)
   */
  public static MarkovBackgroundModel getMarkovModel(BackgroundModelMetadata md) throws SQLException, NotFoundException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getMarkovModel(md, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
    
  /**
   * Get a Markov Background Model described by the specified metadata
   * @param md A complete metadata object specifying a background model. 
   * @param cxn an open db connection to the annotations schema
   * @return A MarkovBackgroundModel, or null if none exists (and can't be generated)
   * @throws SQLException
   * @throws NotFoundException if the genomeID is invalid, most likely because
   * the metadata object wasn't loaded from the database or was modified afterwards
   */
  public static MarkovBackgroundModel getMarkovModel(BackgroundModelMetadata md, Connection cxn) throws SQLException, NotFoundException {
    PreparedStatement getProbs = null;
    ResultSet rs = null;
    try {
      //check that the metadata is complete
      if (!md.hasMapID() || !md.hasGenomeID() || !md.hasModelID() || (md.getName() == null) || !md.hasMaxKmerLen() || (md.getDBModelType() == null)) {
        throw new NullPointerException("Metadata object has one or more null fields");
      }
      
      int mapID = md.getMapID();
      //if the model's db model type is Markov then load it directly
      if (md.getDBModelType().equals(MARKOV_TYPE_STRING)) {        
        MarkovBackgroundModel mbm = new MarkovBackgroundModel(md);
        
        getProbs = cxn.prepareStatement(SQL_GET_MODEL_BY_ID);
        getProbs.setInt(1, mapID);
        rs = getProbs.executeQuery();
        if (rs.next()) {
          BackgroundModelImport.initMarkovModelProbs(mbm, rs, SQL_GET_MODEL_BY_ID);
        }
        return mbm;
      }
      //otherwise, if it's frequency, then load it as a frequency model and 
      //convert it to a markov model
      else if (md.getDBModelType().equals(FREQUENCY_TYPE_STRING)) {
        return new MarkovBackgroundModel(BackgroundModelImport.getFrequencyModel(md, cxn));
      }
      //otherwise there's some error with the model type
      else {
        throw new IllegalArgumentException("Invalid model type: " + md.getDBModelType());
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getProbs != null) {
        getProbs.close();
      }
    }
  }
  
  
  /**
   * Returns a list of Markov Background Models parsed out of a single result
   * set
   * @param rs the result set containing the data for the background model(s)
   * @param queryCore the SQL query core that was used to generate the specified
   * result set
   * @return A list of MarkovBackgroundModels
   * @throws SQLException
   */
  private static List<MarkovBackgroundModel> createMarkovModels(ResultSet rs, String queryCore) throws SQLException {
    int mapIDIndex;
    int nameIndex;
    int kmerLenIndex;
    int modelIDIndex;
    int genomeIDIndex;

    if (queryCore.equals(SQL_GET_MODEL_CORE)) {
      mapIDIndex = SQL_GET_MODEL_CORE_MAP_ID_INDEX;
      nameIndex = SQL_GET_MODEL_CORE_NAME_INDEX;
      kmerLenIndex = SQL_GET_MODEL_CORE_KMERLEN_INDEX;
      modelIDIndex = SQL_GET_MODEL_CORE_MODEL_ID_INDEX;
      genomeIDIndex = SQL_GET_MODEL_CORE_GENOME_ID_INDEX;
    }
    else {
      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
    }     

    try {
      List<MarkovBackgroundModel> models = new ArrayList<MarkovBackgroundModel>();
      if (rs.next()) {
        boolean hasNext = true;
        while (hasNext) {
          int mapID = rs.getInt(mapIDIndex);
          MarkovBackgroundModel mbm = new MarkovBackgroundModel(rs.getString(nameIndex), Organism.findGenome(rs.getInt(genomeIDIndex)), rs.getInt(kmerLenIndex));
          mbm.setMapID(mapID);
          mbm.setModelID(rs.getInt(modelIDIndex));
          models.add(mbm);
          //call the subroutine to parse all the probabilities from the result set
          hasNext = BackgroundModelImport.initMarkovModelProbs(mbm, rs, queryCore);
        }
      }
      return models;
    }
    catch (NotFoundException nfex) {
      throw new DatabaseException("Error loading genome for model", nfex);
    }      
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  private static List<MarkovBackgroundModel> getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getMarkovModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }

  
  /**
   * Given the parameters that describe one or more background models construct
   * an appropriate SQL query and return the Background Models that are found 
   * @param mapID an ID from the background genome map (this will limit to a 
   * single model)
   * @param genomeID a genome ID for the models
   * @param modelID a model ID from the background_model table
   * @param name a model name
   * @param maxKmerLen a maximum kmer length
   * @param cxn an open DB connection to the annotations schema
   * @return A list of matching Markov Background Models
   * @throws SQLException
   */
  private static List<MarkovBackgroundModel> getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn) throws SQLException {
    PreparedStatement getModels = null;
    ResultSet rs = null;
    try {
      //build the appropriate sql query
      StringBuffer sql = new StringBuffer(SQL_GET_MODEL_CORE);
      if (mapID != null) {
        sql.append(SQL_GET_MODEL_MAP_ID);
      }
      if (genomeID != null) {
        sql.append(SQL_GET_MODEL_GENOME_ID);
      }
      if (modelID != null) {
        sql.append(SQL_GET_MODEL_MODEL_ID);
      }
      if (name != null) {
        sql.append(SQL_GET_MODEL_NAME);
      }
      if (maxKmerLen != null) {
        sql.append(SQL_GET_MODEL_KMER_LEN);
      }
      sql.append(SQL_GET_MODEL_ORDER_BY);
      //done building the query      
      
      getModels = cxn.prepareStatement(sql.toString());
      
      //set the sql params
      getModels.setString(1, MARKOV_TYPE_STRING);
      int argCount = 2;
      if (mapID != null) {
        getModels.setInt(argCount, mapID);
        argCount++;
      }
      if (genomeID != null) {
        getModels.setInt(argCount, genomeID);
        argCount++;
      }
      if (modelID != null) {
        getModels.setInt(argCount, modelID);
        argCount++;
      }
      if (name != null) {
        getModels.setString(argCount, name);
        argCount++;
      }
      if (maxKmerLen != null) {
        getModels.setInt(argCount, maxKmerLen);
        argCount++;
      }
      //done setting the params
      
      rs = getModels.executeQuery();
      return BackgroundModelImport.createMarkovModels(rs, SQL_GET_MODEL_CORE);      
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getModels != null) {
        getModels.close();
      }
    }
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(String name) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, name, null);
  }

  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(String name, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, name, null,cxn);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(String name, int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, name, maxKmerLen);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(String name, int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, name, maxKmerLen, cxn);
  }

  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(int modelID) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, modelID, null, null);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(int modelID, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, modelID, null, null, cxn);
  }
  

  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(BackgroundModelMetadata md) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getMarkovModels(md, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }

  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModels(BackgroundModelMetadata md, Connection cxn) throws SQLException {
    Integer mapID = null;
    Integer genomeID = null;
    Integer modelID = null;
    String name = null;
    Integer maxKmerLen = null;
    
    //initialize any fields that are available 
    if (md.hasMapID()) {
      mapID = md.getMapID();
    }
    if (md.hasGenomeID()) {
     genomeID = md.getGenomeID();
    }
    if (md.hasModelID()) {
      modelID = md.getModelID();
    }
    name = md.getName();
    if (md.hasMaxKmerLen()) {
      maxKmerLen = md.getMaxKmerLen();
    }

    return BackgroundModelImport.getMarkovModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByLength(int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, null, maxKmerLen);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByLength(int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, null, maxKmerLen, cxn);
  }
  
  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, null, null);
  }


  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, null, null, cxn);
  }


  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */  
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID, String name) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, name, null);
  }


  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID, String name, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, name, null, cxn);
  }

  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID, int maxKmerLen) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, null, maxKmerLen);
  }
  

  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getMarkovModelsByGenome(int genomeID, int maxKmerLen, Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, genomeID, null, null, maxKmerLen, cxn);
  }

  
  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getAllMarkovModels() throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, null, null);
  }
  

  /**
   * @see getMarkovModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
   */
  public static List<MarkovBackgroundModel> getAllMarkovModels(Connection cxn) throws SQLException {
    return BackgroundModelImport.getMarkovModels(null, null, null, null, null, cxn);
  }

  
  /**************************************************************************
   * Methods for loading Counts background models
   **************************************************************************/
  
  /**
   * 
   * @param mapID
   * @return
   * @throws SQLException
   */

  public static boolean hasCounts(int mapID) throws SQLException{
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.hasCounts(mapID, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
  
  public static boolean hasCounts(int mapID, Connection cxn) throws SQLException {
    PreparedStatement checkCounts = null;
    ResultSet rs = null;
    try {
      //check if the model has any kmers with a null count. 
      checkCounts = cxn.prepareStatement("select bgm.has_counts from background_model bgm, background_genome_map map where map.model_id = bgm.id and map.id = ?");
      checkCounts.setInt(1, mapID);
      rs = checkCounts.executeQuery();      
      if (rs.next()) {
        return false;
      }
      else {
        return true;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (checkCounts != null) {
        checkCounts.close();
      }
    }
  }
  
  
  /**
   * Given a CountsBackgroundModel object and a result set whose rows are 
   * kmer counts (and perhaps other data), parse those rows and 
   * set the model counts accordingly
   * @param fbm A CountsBackgroundModel object to initialize
   * @param rs a ResultSet of kmer probability data, already on the first row
   * for the specified background model
   * @param queryCore the query core used to obtain the result set. This is used
   * to determine which indices to check for the kmer probabilities
   * @return true if there are more rows in the result set for another background
   * model, false otherwise
   * @throws SQLException
   */
  private static boolean initCountsModel(CountsBackgroundModel cbm, ResultSet rs, String queryCore) throws SQLException {
    int idIndex;
    int kmerIndex;
    int countIndex;
    
    if (queryCore.equals(SQL_GET_MODEL_BY_ID)) {
      idIndex = SQL_GET_MODEL_BY_ID_MAP_ID_INDEX;
      kmerIndex = SQL_GET_MODEL_BY_ID_KMER_INDEX;
      countIndex = SQL_GET_MODEL_BY_ID_COUNT_INDEX;
    }
    else {
      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
    }

    boolean hasNext = true;
    do {
      String kmer = rs.getString(kmerIndex);
      cbm.setKmerCount(kmer, rs.getLong(countIndex));
      hasNext = rs.next();
    } while (hasNext && (rs.getInt(idIndex) == cbm.getMapID()));

    return hasNext;
  }


  /**
   * @see getCountsModel(BackgroundModelMetadata md, Connection cxn)
   */
  public static CountsBackgroundModel getCountsModel(BackgroundModelMetadata md) throws SQLException, NotFoundException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      return BackgroundModelImport.getCountsModel(md, cxn);
    }
    finally {
      DatabaseFactory.freeConnection(cxn);
    }
  }
  
    
  /**
   * Get a Counts Background Model described by the specified metadata
   * @param md A complete metadata object specifying a background model.  
   * @param cxn an open db connection to the annotations schema
   * @return A CountsBackgroundModel, or null if none exists (and can't be generated)
   * @throws SQLException
   * @throws NotFoundException if the genomeID is invalid, most likely because
   * the metadata object wasn't loaded from the database or was modified afterwards
   */
  public static CountsBackgroundModel getCountsModel(BackgroundModelMetadata md, Connection cxn) throws SQLException, NotFoundException {
    PreparedStatement getCounts = null;
    ResultSet rs = null;
    try {
      //check that the metadata is complete (except for the model type)
      if (!md.hasMapID() || !md.hasGenomeID() || !md.hasModelID() || (md.getName() == null) || !md.hasMaxKmerLen() || (!md.hasDBModelType())) {
        throw new NullPointerException("Metadata object has one or more null fields");
      }
      
      int mapID = md.getMapID();
      if (BackgroundModelImport.hasCounts(mapID)) {        
        CountsBackgroundModel cbm = new CountsBackgroundModel(md);
        
        getCounts = cxn.prepareStatement(SQL_GET_MODEL_BY_ID);
        getCounts.setInt(1, mapID);
        rs = getCounts.executeQuery();
        if (rs.next()) {
          BackgroundModelImport.initCountsModel(cbm, rs, SQL_GET_MODEL_BY_ID);
        }
        return cbm;
      }      
      //otherwise there's some error with the model type
      else {
        return null;
      }
    }
    finally {
      if (rs != null) {
        rs.close();
      }
      if (getCounts != null) {
        getCounts.close();
      }
    }
  }
  
  
//  /**
//   * Returns a list of Counts Background Models parsed out of a single result
//   * set
//   * @param rs the result set containing the data for the background model(s)
//   * @param queryCore the SQL query core that was used to generate the specified
//   * result set
//   * @return A list of CountsBackgroundModels
//   * @throws SQLException
//   */
//  private static List<CountsBackgroundModel> createCountsModels(ResultSet rs, String queryCore) throws SQLException {
//    int mapIDIndex;
//    int nameIndex;
//    int kmerLenIndex;
//    int modelIDIndex;
//    int genomeIDIndex;
//
//    if (queryCore.equals(SQL_GET_MODEL_CORE)) {
//      mapIDIndex = SQL_GET_MODEL_CORE_MAP_ID_INDEX;
//      nameIndex = SQL_GET_MODEL_CORE_NAME_INDEX;
//      kmerLenIndex = SQL_GET_MODEL_CORE_KMERLEN_INDEX;
//      modelIDIndex = SQL_GET_MODEL_CORE_MODEL_ID_INDEX;
//      genomeIDIndex = SQL_GET_MODEL_CORE_GENOME_ID_INDEX;
//    }
//    else {
//      throw new IllegalArgumentException("Unrecognized query core: " + queryCore);
//    }     
//
//    try {
//      List<CountsBackgroundModel> models = new ArrayList<CountsBackgroundModel>();
//      if (rs.next()) {
//        boolean hasNext = true;
//        while (hasNext) {
//          int mapID = rs.getInt(mapIDIndex);
//          CountsBackgroundModel cbm = new CountsBackgroundModel(rs.getString(nameIndex), Organism.findGenome(rs.getInt(genomeIDIndex)), rs.getInt(kmerLenIndex));
//          cbm.setMapID(mapID);
//          cbm.setModelID(rs.getInt(modelIDIndex));
//          models.add(cbm);
//          //call the subroutine to parse all the probabilities from the result set
//          hasNext = BackgroundModelImport.initCountsModelProbs(cbm, rs, queryCore);
//        }
//      }
//      return models;
//    }
//    catch (NotFoundException nfex) {
//      throw new DatabaseException("Error loading genome for model", nfex);
//    }      
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  private static List<CountsBackgroundModel> getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen) throws SQLException {
//    java.sql.Connection cxn = null;
//    try {
//      cxn = DatabaseFactory.getConnection("annotations");
//      return BackgroundModelImport.getCountsModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
//    }
//    finally {
//      DatabaseFactory.freeConnection(cxn);
//    }
//  }
//
//  
//  /**
//   * Given the parameters that describe one or more background models construct
//   * an appropriate SQL query and return the Background Models that are found 
//   * @param mapID an ID from the background genome map (this will limit to a 
//   * single model)
//   * @param genomeID a genome ID for the models
//   * @param modelID a model ID from the background_model table
//   * @param name a model name
//   * @param maxKmerLen a maximum kmer length
//   * @param cxn an open DB connection to the annotations schema
//   * @return A list of matching Counts Background Models
//   * @throws SQLException
//   */
//  private static List<CountsBackgroundModel> getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn) throws SQLException {
//    PreparedStatement getModels = null;
//    ResultSet rs = null;
//    try {
//      //build the appropriate sql query
//      StringBuffer sql = new StringBuffer(SQL_GET_MODEL_CORE);
//      if (mapID != null) {
//        sql.append(SQL_GET_MODEL_MAP_ID);
//      }
//      if (genomeID != null) {
//        sql.append(SQL_GET_MODEL_GENOME_ID);
//      }
//      if (modelID != null) {
//        sql.append(SQL_GET_MODEL_MODEL_ID);
//      }
//      if (name != null) {
//        sql.append(SQL_GET_MODEL_NAME);
//      }
//      if (maxKmerLen != null) {
//        sql.append(SQL_GET_MODEL_KMER_LEN);
//      }
//      sql.append(SQL_GET_MODEL_ORDER_BY);
//      //done building the query      
//      
//      getModels = cxn.prepareStatement(sql.toString());
//      
//      //set the sql params
//      getModels.setString(1, Counts_TYPE_STRING);
//      int argCount = 2;
//      if (mapID != null) {
//        getModels.setInt(argCount, mapID);
//        argCount++;
//      }
//      if (genomeID != null) {
//        getModels.setInt(argCount, genomeID);
//        argCount++;
//      }
//      if (modelID != null) {
//        getModels.setInt(argCount, modelID);
//        argCount++;
//      }
//      if (name != null) {
//        getModels.setString(argCount, name);
//        argCount++;
//      }
//      if (maxKmerLen != null) {
//        getModels.setInt(argCount, maxKmerLen);
//        argCount++;
//      }
//      //done setting the params
//      
//      rs = getModels.executeQuery();
//      return BackgroundModelImport.createCountsModels(rs, SQL_GET_MODEL_CORE);      
//    }
//    finally {
//      if (rs != null) {
//        rs.close();
//      }
//      if (getModels != null) {
//        getModels.close();
//      }
//    }
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(String name) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, name, null);
//  }
//
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(String name, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, name, null,cxn);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(String name, int maxKmerLen) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, name, maxKmerLen);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(String name, int maxKmerLen, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, name, maxKmerLen, cxn);
//  }
//
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(int modelID) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, modelID, null, null);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(int modelID, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, modelID, null, null, cxn);
//  }
//  
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(BackgroundModelMetadata md) throws SQLException {
//    java.sql.Connection cxn = null;
//    try {
//      cxn = DatabaseFactory.getConnection("annotations");
//      return BackgroundModelImport.getCountsModels(md, cxn);
//    }
//    finally {
//      DatabaseFactory.freeConnection(cxn);
//    }
//  }
//
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModels(BackgroundModelMetadata md, Connection cxn) throws SQLException {
//    Integer mapID = null;
//    Integer genomeID = null;
//    Integer modelID = null;
//    String name = null;
//    Integer maxKmerLen = null;
//    
//    //initialize any fields that are available 
//    if (md.hasMapID()) {
//      mapID = md.getMapID();
//    }
//    if (md.hasGenomeID()) {
//     genomeID = md.getGenomeID();
//    }
//    if (md.hasModelID()) {
//      modelID = md.getModelID();
//    }
//    name = md.getName();
//    if (md.hasMaxKmerLen()) {
//      maxKmerLen = md.getMaxKmerLen();
//    }
//
//    return BackgroundModelImport.getCountsModels(mapID, genomeID, modelID, name, maxKmerLen, cxn);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByLength(int maxKmerLen) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, null, maxKmerLen);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByLength(int maxKmerLen, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, null, maxKmerLen, cxn);
//  }
//  
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, null, null);
//  }
//
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, null, null, cxn);
//  }
//
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */  
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID, String name) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, name, null);
//  }
//
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID, String name, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, name, null, cxn);
//  }
//
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID, int maxKmerLen) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, null, maxKmerLen);
//  }
//  
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getCountsModelsByGenome(int genomeID, int maxKmerLen, Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, genomeID, null, null, maxKmerLen, cxn);
//  }
//
//  
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getAllCountsModels() throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, null, null);
//  }
//  
//
//  /**
//   * @see getCountsModels(Integer mapID, Integer genomeID, Integer modelID, String name, Integer maxKmerLen, Connection cxn)
//   */
//  public static List<CountsBackgroundModel> getAllCountsModels(Connection cxn) throws SQLException {
//    return BackgroundModelImport.getCountsModels(null, null, null, null, null, cxn);
//  }

  
  /**************************************************************************
   * Inserting and Updating Background Models
   **************************************************************************/
  
  /**
   * Insert a markov background model into the database
   * @param model the model to insert
   * @return the DBID of the model
   * @throws SQLException
   */
  public static int insertMarkovModel(MarkovBackgroundModel model) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if ((model.getName() == null) || (model.getName().isEmpty()) || (model.getGenome() == null)) {
      throw new IllegalArgumentException("Model must have a name and genome specified to be imported to database.");
    }
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);

      //insert into the background model and background
      int mapID = BackgroundModelImport.insertBackgroundModelAndMap(model.getName(), model.getMaxKmerLen(), 
          MARKOV_TYPE_STRING, model.getGenome().getDBID(), false, cxn);

      //Finally, insert all the "columns" of the background model
      BackgroundModelImport.insertMarkovModelColumns(model, mapID, cxn);
      
      //If everything has worked then commit
      cxn.commit();

      //update the model with its new ID
      model.setMapID(mapID);

      return mapID;
    }
    catch (SQLException sqlex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw sqlex;      
    }
    catch (CGSException cgsex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw cgsex;
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }


  /**
   * insert a frequency background model into the database
   * @param model the model to insert
   * @return the DBID of the model
   * @throws SQLException
   */
  public static int insertFrequencyModel(FrequencyBackgroundModel model) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if ((model.getName() == null) || (model.getName().isEmpty()) || (model.getGenome() == null)) {
      throw new IllegalArgumentException("Model must have a name and genome specified to be imported to database.");
    }
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);

      //insert into the background model and background
      int mapID = BackgroundModelImport.insertBackgroundModelAndMap(model.getName(), model.getMaxKmerLen(), 
          FREQUENCY_TYPE_STRING, model.getGenome().getDBID(), false, cxn);

      //insert all the "columns" of the background model
      BackgroundModelImport.insertFrequencyModelColumns(model, mapID, cxn);

      //If everything has worked then commit
      cxn.commit();

      //update the model with its new ID
      model.setMapID(mapID);

      return mapID;
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }
  

  /**
   * Insert a counts background model into the database
   * @param model the model to insert
   * @param insertAsMarkov indicates whether to insert the markov probabilities
   * or the frequencies corresponding to the counts 
   * @return the DBID of the model
   * @throws SQLException
   */
  public static int insertCountsModel(CountsBackgroundModel model, boolean insertAsMarkov) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if ((model.getName() == null) || (model.getName().isEmpty()) || (model.getGenome() == null)) {
      throw new IllegalArgumentException("Model must have a name and genome specified to be imported to database.");
    }
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);

      String modelType;
      if (insertAsMarkov) {
        modelType = MARKOV_TYPE_STRING;
      }
      else {
        modelType = FREQUENCY_TYPE_STRING;
      }
      
      //insert into the background model and background
      int mapID = BackgroundModelImport.insertBackgroundModelAndMap(model.getName(), model.getMaxKmerLen(), 
          modelType, model.getGenome().getDBID(), true, cxn);

      //insert all the "columns" of the background model
      BackgroundModelImport.insertCountsModelColumns(model, mapID, insertAsMarkov, cxn);

      //If everything has worked then commit
      cxn.commit();

      //update the model with its new ID
      model.setMapID(mapID);

      return mapID;
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }
  
  
  /**
   * Update a markov background model already in the database
   * @param model the model to insert
   * @throws SQLException
   */
  public static void updateMarkovModel(MarkovBackgroundModel model) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if (!model.hasMapID()) {
      throw new IllegalArgumentException("Model must already have a database ID to be updated in the database.");
    }
    
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);

      int mapID = model.getMapID();
      //remove from the database all the existing entries for the model columns
      BackgroundModelImport.removeModelColumns(mapID, cxn);
      
      //Insert all the "columns" of the background model
      BackgroundModelImport.insertMarkovModelColumns(model, mapID, cxn);
      
      //If everything has worked then commit
      cxn.commit();
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }

  
  /**
   * Update a frequency background model already in the database
   * @param model the model to insert
   * @throws SQLException
   */
  public static void updateFrequencyModel(FrequencyBackgroundModel model) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if (!model.hasMapID()) {
      throw new IllegalArgumentException("Model must already have a database ID to be updated in the database.");
    }
    
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);

      int mapID = model.getMapID();
      //remove from the database all the existing entries for the model columns
      BackgroundModelImport.removeModelColumns(mapID, cxn);
      
      //Insert all the "columns" of the background model
      BackgroundModelImport.insertFrequencyModelColumns(model, mapID, cxn);
      
      //If everything has worked then commit
      cxn.commit();
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }

  
  /**
   * Update a counts background model already in the database
   * @param model the model to insert
   * @throws SQLException
   */
  public static void updateCountsModel(CountsBackgroundModel model) throws SQLException, CGSException {
    //make sure the model has a name and genome
    if (!model.hasMapID()) {
      throw new IllegalArgumentException("Model must already have a database ID to be updated in the database.");
    }
    
    java.sql.Connection cxn = null;
    PreparedStatement getModelType = null;
    ResultSet rs = null;
    try {
      cxn = DatabaseFactory.getConnection("annotations");
      cxn.setAutoCommit(false);      
      
      int mapID = model.getMapID();
      
      //determine whether this model exists in the database as a markov model or
      //a frequency model, so that it can be updated in the same format
      boolean isMarkov;
      getModelType = 
      	cxn.prepareStatement("select model_type from background_model bm, background_genome_map map where map.id = ? and map.bg_model_id = bm.id");
      getModelType.setInt(1, mapID);
      rs = getModelType.executeQuery();
      if (rs.next()) {
      	isMarkov = rs.getString(1).equals(MARKOV_TYPE_STRING);
      }
      else {
      	throw new DatabaseException("Unable to find Background Model in database.");
      }
      
      //remove from the database all the existing entries for the model columns
      BackgroundModelImport.removeModelColumns(mapID, cxn);
      
      //Insert all the "columns" of the background model
      BackgroundModelImport.insertCountsModelColumns(model, mapID, isMarkov, cxn);
      
      //If everything has worked then commit
      cxn.commit();
    }
    catch (RuntimeException ex) {
      //If any runtime exceptions come up rollback the transaction and then
      //rethrow the exception
      cxn.rollback();
      throw ex;
    }
    finally {
    	if (rs != null) {
    		rs.close();
    	}
    	if (getModelType != null) {
    		getModelType.close();
    	}
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }

  
  /**
   * insert entries for the model in the background model table and the
   * background model genome map table
   * @param name the name of the model
   * @param kmerLen the length of the longest kmer in the model
   * @param dbModelType the type of the model ("MARKOV" or "FREQUENCY")
   * @param genomeID the DBID of the genome the model is for
   * @param cxn an open database connection
   * @return the background genome map ID of the model
   * @throws SQLException
   */
  private static Integer insertBackgroundModelAndMap(String name, int kmerLen, String modelType, int genomeID, boolean hasCounts, Connection cxn) throws SQLException, CGSException {
    //FIXME make this return a pair for the modelID and the mapID
    /**
     * Check whether there is already an entry for a model with this name, maxKmerLen, and type. If so, reuse the model ID, 
     * otherwise create one.
     */
    Integer modelID = BackgroundModelImport.getBackgroundModelID(name, kmerLen, modelType, cxn);
    if (modelID == null) {
      modelID = BackgroundModelImport.insertBackgroundModel(name, kmerLen, modelType, cxn);
    }

    /**
     * Check whether there is already an entry for this model's genome in the
     * background model genome map. If not then create one, otherwise throw
     * an exception indicating that the update method should be called to
     * update an existing model.
     */
    Integer mapID = BackgroundModelImport.getBackgroundGenomeMapID(modelID, genomeID, cxn);
    if (mapID != null) {
      throw new CGSException("Model already exists. Select a different name or use updateMarkovModel() to update the existing model.");
    }
    else {
      mapID = BackgroundModelImport.insertBackgroundGenomeMap(modelID, genomeID, hasCounts, cxn);
    }

    return mapID;
  }


  /**
   * insert the model in the the background model table
   * @param name the name of the model
   * @param kmerLen the length of the longest kmer in the model
   * @param dbModelType the type of the model ("MARKOV" or "FREQUENCY")
   * @param cxn an open database connection
   * @return the background model ID 
   * @throws SQLException
   */
  private static Integer insertBackgroundModel(String name, int kmerLen, String modelType, Connection cxn) throws SQLException {
    PreparedStatement insertBG = cxn.prepareStatement("insert into background_model(id,name,max_kmer_len,model_type) values (weightmatrix_id.nextval,?,?,?)");
    insertBG.setString(1, name);
    insertBG.setInt(2, kmerLen);
    insertBG.setString(3, modelType);
    insertBG.execute();
    insertBG.close();
    PreparedStatement getModelID = cxn.prepareStatement("select weightmatrix_id.currval from dual");
    ResultSet rs = getModelID.executeQuery();
    try {
      if (rs.next()) {
        return rs.getInt(1);
      }
      else {
        throw new SQLException("Failed to get Model ID following insert into Background_Model table.");
      }    
    }
    finally {
      rs.close();
      getModelID.close();
    }
  }
  
  
  /**
   * insert the model in the background model genome map table
   * @param modelID the background model table ID of the model
   * @param genomeID the DBID of the genome the model is for
   * @param cxn an open database connection
   * @return the
   * @throws SQLException
   */
  private static Integer insertBackgroundGenomeMap(int modelID, int genomeID, boolean hasCounts, Connection cxn) throws SQLException {
    PreparedStatement insertMap = 
      cxn.prepareStatement("insert into background_genome_map(id, genome_id, bg_model_id, has_counts) values (weightmatrix_id.nextval, ?, ?, ?)");
    insertMap.setInt(1, genomeID);
    insertMap.setInt(2, modelID);
    if (hasCounts) {
      insertMap.setInt(3, 1);      
    }
    else {
      insertMap.setInt(3, 0);
    }
    
    insertMap.execute();
    insertMap.close();
    PreparedStatement getMapID = cxn.prepareStatement("select weightmatrix_id.currval from dual");
    ResultSet rs = getMapID.executeQuery();
    try {
      if (rs.next()) {
        return rs.getInt(1);
      }
      else {
        throw new SQLException("Failed to get Background Model Genome Map ID following insert into Background_Genome_Map table.");
      }
    }
    finally {
      rs.close();
      getMapID.close();
    }
  }
  
  
  /**
   * remove the model's kmer probabilities 
   * @param mapID the model's ID in the background genome map table
   * @param cxn an open database connection
   * @throws SQLException
   */
  private static void removeModelColumns(int mapID, Connection cxn) throws SQLException {
  	PreparedStatement deleteOld = cxn.prepareStatement("delete from background_model_cols where map_id = ?");
    try {
      deleteOld.setInt(1,mapID);
      deleteOld.execute();
    }
    finally {
      deleteOld.close();
    }
  }
  
  
  /**
   * insert the model's kmer probabilities 
   * @param model the background model
   * @param mapID the model's ID in the background genome map table
   * @param cxn an open database connection
   * @throws SQLException
   */
  private static void insertMarkovModelColumns(MarkovBackgroundModel model, int mapID, Connection cxn) throws SQLException {
    PreparedStatement insertCol = cxn.prepareStatement("insert into background_model_cols(map_id,kmer,probability) values(?,?,?)");
    try {
      for (int i = 1; i <= model.getMaxKmerLen(); i++) {
        for (String kmer : model.getKmers(i)) {
          double prob = model.getMarkovProb(kmer);

          insertCol.setInt(1, mapID);
          insertCol.setString(2, kmer);
          insertCol.setDouble(3, prob);
          insertCol.execute();
        }
      }
    }
    finally {
      insertCol.close();
    }
  }
  
  
  /**
   * insert the model's kmer probabilities 
   * @param model the background model
   * @param mapID the model's ID in the background genome map table
   * @param cxn an open database connection
   * @throws SQLException
   */
  private static void insertFrequencyModelColumns(FrequencyBackgroundModel model, int mapID, Connection cxn) throws SQLException {
    PreparedStatement insertCol = cxn.prepareStatement("insert into background_model_cols(map_id,kmer,probability) values(?,?,?)");
    try {
      for (int i = 1; i <= model.getMaxKmerLen(); i++) {
        for (String kmer : model.getKmers(i)) {
          double prob = model.getFrequency(kmer);

          insertCol.setInt(1, mapID);
          insertCol.setString(2, kmer);
          insertCol.setDouble(3, prob);
          insertCol.execute();
        }
      }
    }
    finally {
      insertCol.close();
    }
  }
  

  /**
   * insert the model's kmer probabilities 
   * @param model the background model
   * @param mapID the model's ID in the background genome map table
   * @param insertAsMarkov if true then insert the markov probabilities, if
   * false then insert the kmer frequencies
   * @param cxn an open database connection
   * @throws SQLException
   */
  private static void insertCountsModelColumns(CountsBackgroundModel model, int mapID, boolean insertAsMarkov, Connection cxn) throws SQLException {
    PreparedStatement insertCol = cxn.prepareStatement("insert into background_model_cols(map_id,kmer,probability,count) values(?,?,?,?)");
    try {
      for (int i = 1; i <= model.getMaxKmerLen(); i++) {
        for (String kmer : model.getKmers(i)) {
          double prob;
          if (insertAsMarkov) {
            prob = model.getMarkovProb(kmer);
          }
          else {
            prob = model.getFrequency(kmer);
          }          

          insertCol.setInt(1, mapID);
          insertCol.setString(2, kmer);
          insertCol.setDouble(3, prob);
          insertCol.setLong(4, model.getKmerCount(kmer));
          insertCol.execute();
        }
      }
    }
    finally {
      insertCol.close();
    }
  }
  
  
  public static void main(String[] args) {
    BackgroundModelImport.testDBLoading();
    //BackgroundModelImport.importModel(args);
  }
  
  
  /** 
   * Imports a weightmatrix from a file into the DB 
   *
   * Usage:
   * java edu.mit.csail.cgs.datasets.motifs.BackgroundModelImport --genome "Mus musculus;mm8" --bgname "whole genome" --bgtype MARKOV --bgfile foo.back
   */
  public static void importModel(String args[]) {
    try {
      Genome gen = null;
      String bgName = null;
      String bgType = null;
      String bgFilename = null;
      
      //args = new String[] {"--species", "Mus musculus;mm5", "--bgname", "test1", "--bgtype", MARKOV_TYPE_STRING, "--bgfile", "mm8_2.back"};
//      args = new String[] {"--species", "Saccharomyces cerevisiae;sacCer1", "--bgname", "test2", "--bgtype", MARKOV_TYPE_STRING, "--bgfile", "yeast1.back"};
      args = new String[] {"--species", "Homo sapiens;hg17", "--bgname", "test1", "--bgtype", MARKOV_TYPE_STRING, "--bgfile", "human_1.back"};
      
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

      if (bgType.toUpperCase().equals(MARKOV_TYPE_STRING)) {
        MarkovBackgroundModel mbg = BackgroundModelIO.parseMarkovBackgroundModel(bgName, bgFilename, gen);
        BackgroundModelImport.insertMarkovModel(mbg);
      }
      else if (bgType.toUpperCase().equals(FREQUENCY_TYPE_STRING)) {
        FrequencyBackgroundModel fbg = BackgroundModelIO.parseFreqBackgroundModel(bgName, bgFilename, gen);
        BackgroundModelImport.insertFrequencyModel(fbg);
      }
      else if (bgType.toUpperCase().equals("COUNTS;MARKOV")) {
        CountsBackgroundModel cbg = BackgroundModelIO.parseCountsBackgroundModel(bgName, bgFilename, gen);
        BackgroundModelImport.insertCountsModel(cbg, true);
      }
      else if (bgType.toUpperCase().equals("COUNTS;FREQUENCY")) {
        CountsBackgroundModel cbg = BackgroundModelIO.parseCountsBackgroundModel(bgName, bgFilename, gen);
        BackgroundModelImport.insertCountsModel(cbg, false);
      }
      else {
        logger.fatal("Background type must be one of: Markov, Frequency, Counts;Markov, Counts;Frequency.");
        System.exit(1);
      }
    }
    catch (NotFoundException nfex) {
      logger.fatal(nfex);
      System.exit(1);
    }
    catch (IOException ioex) {
      logger.fatal(ioex);
      System.exit(1);
    }
    catch (ParseException pex) {
      logger.fatal(pex);      
      System.exit(1);
    }
    catch (SQLException sqlex) {
      logger.fatal(sqlex);
      System.exit(1);
    }
    catch (CGSException cgsex) {
      logger.fatal(cgsex);
      System.exit(1);
    }
    System.exit(0);
  }

  
  private static void testDBLoading() {
    ClassLoader loader = BackgroundModelImport.class.getClassLoader();
    PropertyConfigurator.configure(loader.getResource("edu/mit/csail/cgs/utils/config/log4j.properties"));      
    
    MarkovBackgroundModel testValuesMM8_3 = null;
    MarkovBackgroundModel testValuesMM8_2 = null;
    MarkovBackgroundModel testValuesYeast_3 = null;
    MarkovBackgroundModel testValuesYeast_2 = null;
    MarkovBackgroundModel testValuesHuman_2 = null;
    try {
      testValuesMM8_3 = BackgroundModelIO.parseMarkovBackgroundModel("mm8.back", Organism.findGenome("mm8"));
      testValuesMM8_2 = BackgroundModelIO.parseMarkovBackgroundModel("mm8_1.back", Organism.findGenome("mm8"));
      testValuesYeast_3 = BackgroundModelIO.parseMarkovBackgroundModel("yeast2.back", Organism.findGenome("sacCer1"));
      testValuesYeast_2 = BackgroundModelIO.parseMarkovBackgroundModel("yeast1.back", Organism.findGenome("sacCer1"));
      testValuesHuman_2 = BackgroundModelIO.parseMarkovBackgroundModel("human_1.back", Organism.findGenome("hg17"));
    }
    catch (IOException ex) {
      // TODO Auto-generated catch block
      ex.printStackTrace();
    }
    catch (ParseException ex) {
      // TODO Auto-generated catch block
      ex.printStackTrace();
    }
    catch (NotFoundException nfex) {
      // TODO Auto-generated catch block
      nfex.printStackTrace();
    }

    
    try {
      //test1: 5370, 22, 5369, test1, 3, MARKOV
      //test1: 5372, 22, 5371, test1, 2, MARKOV
      //test2: 5374, 23, 5373, test2, 3, MARKOV
      //test2: 5376, 23, 5375, test2, 2, MARKOV
      
      int expectedMapID_test1_2 = 5372;
      int expectedGenomeID_test1_2 = 22;
      int expectedModelID_test1_2 = 5371;
      String expectedName_test1_2 = "test1";
      int expectedKmerLen_test1_2 = 2;
      String expectedType_test1_2 = MARKOV_TYPE_STRING;

      logger.debug("Testing getAllBackgroundModels()");
      List<BackgroundModelMetadata> foo = BackgroundModelImport.getAllBackgroundModels();
      for (BackgroundModelMetadata md : foo) {
        System.out.println(md.toString());        
      }
      
      logger.debug("Testing getAllBackgroundModels(ignoreGenome)");
      foo = BackgroundModelImport.getAllBackgroundModels(true);
      for (BackgroundModelMetadata md : foo) {
        System.out.println(md.toString());        
      }
      
      logger.debug("Testing getBackgroundModelsForGenome(genomeID)");
      foo = BackgroundModelImport.getBackgroundModelsForGenome(expectedGenomeID_test1_2); 
      for (BackgroundModelMetadata md : foo) {
        System.out.println(md.toString());        
      }
      
      logger.debug("Testing getGenomesForBackgroundModel(modelID)");
      List<Integer> bar = BackgroundModelImport.getGenomesForBackgroundModel(expectedModelID_test1_2); //mouse
      for (int genomeID : bar) {
        System.out.println(genomeID);        
      }
      
      
      BackgroundModelMetadata expectedMDModelOnly_test1_2 = new BackgroundModelMetadata(expectedModelID_test1_2, expectedName_test1_2, expectedKmerLen_test1_2, expectedType_test1_2);
      BackgroundModelMetadata expectedMD_test1_2 = new BackgroundModelMetadata(expectedMapID_test1_2, expectedGenomeID_test1_2, 
          expectedModelID_test1_2, expectedName_test1_2, 
          expectedKmerLen_test1_2, expectedType_test1_2, false);
      
      
      int expectedMapID_test1_3 = 5370;
      int expectedGenomeID_test1_3 = 22;
      int expectedModelID_test1_3 = 5369;
      String expectedName_test1_3 = "test1";
      int expectedKmerLen_test1_3 = 3;
      String expectedType_test1_3 = MARKOV_TYPE_STRING;
      
      BackgroundModelMetadata expectedMDModelOnly_test1_3 = new BackgroundModelMetadata(expectedModelID_test1_3, expectedName_test1_3, expectedKmerLen_test1_3, expectedType_test1_3);
      BackgroundModelMetadata expectedMD_test1_3 = new BackgroundModelMetadata(expectedMapID_test1_3, expectedGenomeID_test1_3, 
          expectedModelID_test1_3, expectedName_test1_3, 
          expectedKmerLen_test1_3, expectedType_test1_3, false);
      
      int expectedMapID_test2_3 = 5374;
      int expectedGenomeID_test2_3 = 23;
      int expectedModelID_test2_3 = 5373;
      String expectedName_test2_3 = "test2";
      int expectedKmerLen_test2_3 = 3;
      String expectedType_test2_3 = MARKOV_TYPE_STRING;
      
      BackgroundModelMetadata expectedMDModelOnly_test2_3 = new BackgroundModelMetadata(expectedModelID_test2_3, expectedName_test2_3, expectedKmerLen_test2_3, expectedType_test2_3);
      BackgroundModelMetadata expectedMD_test2_3 = new BackgroundModelMetadata(expectedMapID_test2_3, expectedGenomeID_test2_3, 
          expectedModelID_test2_3, expectedName_test2_3, 
          expectedKmerLen_test2_3, expectedType_test2_3, false);
      
      
      
      //test getBackgroundModelID(name, kmerlen, type)
      BackgroundModelImport.getBackgroundModelID(expectedName_test1_3, expectedKmerLen_test1_3, expectedType_test1_3);
      if (BackgroundModelImport.getBackgroundModelID(expectedName_test1_3, expectedKmerLen_test1_3, expectedType_test1_3) == expectedModelID_test1_3) {
        logger.debug("OK: getBackgroundModelID(name, kmerlen, type)");        
      }
      else {
        logger.debug("FAIL: getBackgroundModelID(name, kmerlen, type)");
      }

      //test getBackgroundModelID(name, kmerlen, type)
      if (BackgroundModelImport.getBackgroundGenomeMapID(expectedModelID_test1_3, expectedGenomeID_test1_3) == expectedMapID_test1_3) {
        logger.debug("OK: getBackgroundGenomeMapID(modelID, genomeID)");        
      }
      else {
        logger.debug("FAIL: getBackgroundGenomeMapID(modelID, genomeID)");
      }

      //test getBackgroundModelByModelID(modelID)
      if (expectedMDModelOnly_test1_3.equals(BackgroundModelImport.getBackgroundModelByModelID(expectedModelID_test1_3))) {
        logger.debug("OK: getBackgroundModelByModelID(modelID)");        
      }
      else {
        logger.debug("FAIL: getBackgroundModelByModelID(modelID)");
      }

      //test getBackgroundModel(modelID, genomeID)
      if (expectedMD_test1_3.equals(BackgroundModelImport.getBackgroundModel(expectedModelID_test1_3, expectedGenomeID_test1_3))) {
        logger.debug("OK: getBackgroundModel(modelID, genomeID)");        
      }
      else {
        logger.debug("FAIL: getBackgroundModel(modelID, genomeID)");
      }

      //test getBackgroundModel(name, kmerlen, type, genomeID)
      if (expectedMD_test1_3.equals(BackgroundModelImport.getBackgroundModel(expectedName_test1_3, expectedKmerLen_test1_3, expectedType_test1_3, expectedGenomeID_test1_3))) {
        logger.debug("OK: getBackgroundModel(name, kmerlen, type, genomeID)");        
      }
      else {
        logger.debug("FAIL: getBackgroundModel(name, kmerlen, type, genomeID)");
      }

      //test getBackgroundModelByMapID(mapID)
      if (expectedMD_test1_3.equals(BackgroundModelImport.getBackgroundModelByMapID(expectedMapID_test1_3))) {
        logger.debug("OK: getBackgroundModelByMapID(mapID)");        
      }
      else {
        logger.debug("FAIL: getBackgroundModelByMapID(mapID)");
      }

      /**********************************************************************
       * Done with basic tests for getting ids and metadata
       **********************************************************************/

      /**********************************************************************
       * Test methods for getting Markov Models
       **********************************************************************/
      
      
      //test getMarkovModel(metadata)
      MarkovBackgroundModel mbm = BackgroundModelImport.getMarkovModel(expectedMD_test1_3);
      if (mbm.equalValues(testValuesMM8_3)) {
        logger.debug("OK: getMarkovModel(metadata)");        
      }
      else {
        logger.debug("FAIL: getmarkovModel(metadata)");
      }
      MarkovBackgroundModel mbm2 = BackgroundModelImport.getMarkovModel(expectedMD_test2_3);
      if (mbm2.equalValues(testValuesYeast_3)) {
        logger.debug("OK: getMarkovModel(metadata)");        
      }
      else {
        logger.debug("FAIL: getmarkovModel(metadata)");
      }
      
      //test getMarkovModels methods
      List<MarkovBackgroundModel> mbmList;      
      
      mbmList = BackgroundModelImport.getMarkovModels(expectedMDModelOnly_test1_2);
      if ((mbmList.size() == 2) && mbmList.get(0).equalValues(testValuesHuman_2) && mbmList.get(1).equalValues(testValuesMM8_2)) {          
        logger.debug("OK: getMarkovModels(metadata)");        
      }
      else {
        logger.debug("FAIL: getmarkovModels(metadata)");
      }

      mbmList = BackgroundModelImport.getMarkovModels("test1");
      if ((mbmList.size() == 3) && mbmList.get(1).equalValues(testValuesMM8_2) && mbmList.get(2).equalValues(testValuesMM8_3)
          && mbmList.get(0).equalValues(testValuesHuman_2)) { 
        logger.debug("OK: getMarkovModels(name)");        
      }
      else {
        logger.debug("FAIL: getmarkovModels(name)");
      }
      
      mbmList = BackgroundModelImport.getMarkovModels("test1", 2);
      if ((mbmList.size() == 2) && mbmList.get(1).equalValues(testValuesMM8_2)) { 
        logger.debug("OK: getMarkovModels(name, kmerlen)");        
      }
      else {
        logger.debug("FAIL: getmarkovModels(name, kmerlen)");
      }
      
      mbmList = BackgroundModelImport.getMarkovModels(expectedModelID_test1_2);
      if ((mbmList.size() == 2) && mbmList.get(0).equalValues(testValuesHuman_2) && mbmList.get(1).equalValues(testValuesMM8_2)) { 
        logger.debug("OK: getMarkovModels(modelID)");        
      }
      else {
        logger.debug("FAIL: getmarkovModels(modelID)");
      }

      mbmList = BackgroundModelImport.getMarkovModelsByLength(2);
      if ((mbmList.size() == 3) && mbmList.get(1).equalValues(testValuesMM8_2) && mbmList.get(2).equalValues(testValuesYeast_2)
          && mbmList.get(0).equalValues(testValuesHuman_2)) { 
        logger.debug("OK: getMarkovModelsByLength(kmerlen)");        
      }
      else {
        logger.debug("FAIL: getmarkovModelsByLength(kmerlen)");
      }
      
      mbmList = BackgroundModelImport.getMarkovModelsByGenome(expectedGenomeID_test1_2);
      if ((mbmList.size() == 2) && mbmList.get(0).equalValues(testValuesMM8_2) && mbmList.get(1).equalValues(testValuesMM8_3)) {          
        logger.debug("OK: getMarkovModelsByGenome(genomeID)");        
      }
      else {
        logger.debug("FAIL: getmarkovModelsByGenome(genomeID)");
      }
      
      mbmList = BackgroundModelImport.getMarkovModelsByGenome(expectedGenomeID_test1_2, expectedName_test1_2);
      if ((mbmList.size() == 2) && mbmList.get(0).equalValues(testValuesMM8_2) && mbmList.get(1).equalValues(testValuesMM8_3)) {          
        logger.debug("OK: getMarkovModelsByGenome(genomeID, name)");        
      }
      else {
        logger.debug("FAIL: getmarkovModelsByGenome(genomeID, name)");
      }

      mbmList = BackgroundModelImport.getMarkovModelsByGenome(expectedGenomeID_test1_2, expectedKmerLen_test1_2);
      if ((mbmList.size() == 1) && mbmList.get(0).equalValues(testValuesMM8_2)) {          
        logger.debug("OK: getMarkovModelsByGenome(genomeID, kmerlen)");        
      }
      else {
        logger.debug("FAIL: getmarkovModelsByGenome(genomeID, kmerlen)");
      }
    }
    catch (NotFoundException nfex) {
      // TODO Auto-generated catch block
      nfex.printStackTrace();
    }
    catch (SQLException ex) {
      // TODO Auto-generated catch block
      ex.printStackTrace();
    }
  }
}
