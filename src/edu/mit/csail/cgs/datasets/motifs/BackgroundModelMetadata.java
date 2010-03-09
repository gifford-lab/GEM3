/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

/**
 * @author rca
 * This class is used by the BackgroundModelLoader for returning the metadata
 * for background models that are in the database.  
 */
public class BackgroundModelMetadata {

  /**
   * fields from the background model table also used when constructing parsed
   * background models
   */
  protected String name;
  protected final int maxKmerLen;

  
  /**
   * fields from the background model table that may be null/-1 for parsed 
   * background models
   */
  protected int modelID;
  protected String dbModelType; //only null for parsed Count-based models
  
  /**
   * these fields may be null/-1 for parsed background models or metadata 
   * instances that just represent a row from the background model table
   */
  protected int mapID;
  protected int genomeID; 


  public BackgroundModelMetadata(String name, int kmerlen, int genomeID) {
    this.name = name;
    this.maxKmerLen = kmerlen;
    this.genomeID = genomeID;
    this.modelID = -1;
    this.dbModelType = null;
    this.mapID = -1;    
  }
  
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType) {
    this(modelID, name, kmerlen, modelType, -1, -1);
  }
  
  
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType, int mapID, int genomeID) {
    this.modelID = modelID;
    this.name = name;
    this.maxKmerLen = kmerlen;
    this.dbModelType = modelType;
    this.mapID = mapID;
    this.genomeID = genomeID;
  }
  
  
  /**
   * Returns the name of this model
   * @return
   */
  public String getName() {
    return name;
  }
  
  
  /**
   * Set the name of this model
   * @param name
   */
  public void setName(String name) {
    this.name = name;
  }
  
  
  /**
   * Returns the max kmer length in this model
   * @return
   */
  public int getMaxKmerLen() {
    return maxKmerLen;
  }

  
  /**
   * Returns the db model type ("FREQUENCY" or "MARKOV"). This may be null for
   * parsed Counts models
   * @return
   */
  public String getDbModelType() {
    return dbModelType;
  }

  
  /**
   * Sets the db model type for this model ("FREQUENCY" or "MARKOV").
   * @param dbModelType
   */
  public void setDbModelType(String dbModelType) {
    if (dbModelType.equals("MARKOV") || dbModelType.equals("FREQUENCY")) {
      this.dbModelType = dbModelType;
    }
    else {
      throw new IllegalArgumentException("Invalid model type: " + dbModelType + 
          ". DB Model type must be either FREQUENCY or MARKOV");
    }    
  }


  /**
   * Returns whether this object has a model ID in the background model table
   * @return
   */
  public boolean hasModelID() {
    return (modelID != -1);
  }
  
  
  /**
   * Returns the id from the background model table 
   * @return
   */
  public int getModelID() {
    return modelID;
  }
  
  
  /**
   * Sets the id for the background model table
   * @param modelID
   */
  public void setModelID(Integer modelID) {
    this.modelID = modelID;
  }

  
  /**
   * Returns true if this model has a database ID
   * @return
   */
  public boolean hasMapID() {
    return (mapID != -1);
  }

  
  /**
   * Returns this model's background genome map ID, which is -1 if it doesn't 
   * have one
   * @return
   */
  public int getMapID() {
    return mapID;
  }

  /**
   * Set the background genome map ID
   * @param mapID
   */
  public void setMapID(int mapID) {    
    this.mapID = mapID;
  }

  
  public Integer getGenomeID() {
    return genomeID;
  }

  public void setGenomeID(Integer genomeID) {
    this.genomeID = genomeID;
  }  
}
