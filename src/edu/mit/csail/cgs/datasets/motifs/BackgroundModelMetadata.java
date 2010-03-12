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
  protected int maxKmerLen;

  
  /**
   * fields from the background model table that may be null for parsed 
   * background models
   */
  protected int modelID;
  protected String dbModelType; //only null for parsed Count-based models
  
  /**
   * these fields may be null for parsed background models or metadata 
   * instances that just represent a row from the background model table
   */
  protected int mapID;
  protected int genomeID; 


  public BackgroundModelMetadata(BackgroundModelMetadata md) {
    this(md.getMapID(), md.getGenomeID(), md.getModelID(), md.getName(), md.getMaxKmerLen(), md.getDBModelType());
  }
  
  public BackgroundModelMetadata(String name, int kmerlen, int genomeID) {
    this(-1, genomeID, -1, name, kmerlen, null);
  }
  
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType) {
    this(-1, -1, modelID, name, kmerlen, modelType);
  }
  
  
  public BackgroundModelMetadata(int mapID, int genomeID, int modelID, String name, int kmerlen, String modelType) {
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
    if (name == null) {
      throw new NullPointerException("Name can not be set to null");
    }
    else {
      this.name = name;
    }
  }
  
  
  public boolean hasMaxKmerLen() {
    return (maxKmerLen != -1);
  }
  
  
  /**
   * Returns the max kmer length in this model
   * @return
   */
  public int getMaxKmerLen() {
    return maxKmerLen;
  }

  
  public void setMaxKmerLen(int maxKmerLen) {
    this.maxKmerLen = maxKmerLen;
  }
  
  
  public boolean hasDBModelType() {
    return (dbModelType != null);
  }
  
  /**
   * Returns the db model type ("FREQUENCY" or "MARKOV"). This may be null for
   * parsed Counts models
   * @return
   */
  public String getDBModelType() {
    return dbModelType;
  }

  
  /**
   * Sets the db model type for this model ("FREQUENCY" or "MARKOV").
   * @param dbModelType
   */
  public void setDBModelType(String dbModelType) {
    if (dbModelType.equals(BackgroundModelImport.MARKOV_TYPE_STRING) || dbModelType.equals(BackgroundModelImport.FREQUENCY_TYPE_STRING)) {
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
  public void setModelID(int modelID) {
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

  
  public boolean hasGenomeID() {
    return (genomeID != -1);
  }
  
  public int getGenomeID() {
    return genomeID;
  }

  public void setGenomeID(int genomeID) {
    this.genomeID = genomeID;
  }
  
  public boolean equals(BackgroundModelMetadata other) {
    if ((other != null) && this.name.equals(other.getName()) && (this.maxKmerLen == other.getMaxKmerLen()) 
        && (((this.dbModelType == null) && (other.getDBModelType() == null)) || this.dbModelType.equals(other.dbModelType)) 
        && (this.modelID == other.getModelID() && this.mapID == other.getMapID())
        && (this.genomeID == other.getGenomeID())) {
      return true;
    }
    else {
      return false;
    }
  }
  
  public String toString() {
    return mapID + "\t" + genomeID + "\t" + modelID + "\t" + name  + "\t" + dbModelType + "\t" + maxKmerLen; 
  }
}
