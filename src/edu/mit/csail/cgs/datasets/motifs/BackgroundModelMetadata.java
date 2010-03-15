/**
 * 
 */
package edu.mit.csail.cgs.datasets.motifs;

/**
 * @author rca
 * This class is a lightweight object that contains the basic fields describing 
 * a background model. It is primarily used by the BackgroundModelLoader for 
 * returning the metadata for background models that are in the database. It 
 * is also the superclass for the BackgroundModel class. It intentially doesn't
 * constrain what values can be set for any of the fields except dbModelType - GIGO.
 * 
 * Note: -1 is used as a value for numeric fields that are unknown or not set.
 */
public class BackgroundModelMetadata {

  /**
   * fields from the background model table also used when constructing parsed
   * background models
   */
  protected String name;
  protected int maxKmerLen;

  
  /**
   * fields from the background model table that may be null/-1 for parsed 
   * background models
   */
  protected int modelID;
  protected String dbModelType; //only null for parsed Count-based models
  
  /**
   * these fields may be -1 for parsed background models or metadata 
   * instances that just represent a row from the background model table
   */
  protected int mapID;
  protected int genomeID; 
  
  protected boolean hasCounts;

  public BackgroundModelMetadata(BackgroundModelMetadata md) {
    this(md.getMapID(), md.getGenomeID(), md.getModelID(), md.getName(), md.getMaxKmerLen(), md.getDBModelType(), md.hasCounts());
  }
  
  //FIXME when is this used?
  public BackgroundModelMetadata(int genomeID, String name, int kmerlen) {
    this(-1, genomeID, -1, name, kmerlen, null, false);
  }
  
  //FIXME when is this used?
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType) {
    this(-1, -1, modelID, name, kmerlen, modelType, false);
  }
  
  public BackgroundModelMetadata(int mapID, int genomeID, int modelID, String name, int kmerlen, String modelType, boolean hasCounts) {
    this.modelID = modelID;
    this.name = name;
    this.maxKmerLen = kmerlen;
    this.dbModelType = modelType;
    this.mapID = mapID;
    this.genomeID = genomeID;
    this.hasCounts = hasCounts;
  }
  
  
  /**
   * Returns true if this object has an ID for the background genome map table
   * @return
   */
  public boolean hasMapID() {
    return (mapID != -1);
  }

  
  /**
   * Returns this object's background genome map ID
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

  
  /**
   * Checks whether this object has a database ID for a genome
   * @return
   */
  public boolean hasGenomeID() {
    return (genomeID != -1);
  }
  
  /**
   * Returns this object's genome ID
   * @return
   */
  public int getGenomeID() {
    return genomeID;
  }

  
  /**
   * Sets this object's genome ID
   * @param genomeID
   */
  public void setGenomeID(int genomeID) {
    this.genomeID = genomeID;
  }
  
  
  /**
   * Checks whether this object has an ID in the background model table
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
   * Returns the name of this object
   * @return
   */
  public String getName() {
    return name;
  }
  
  
  /**
   * Set the name of this object
   * @param name
   */
  public void setName(String name) {
    this.name = name;
  }
  
  
  /**
   * Check whether this object has a max kmer length set
   * @return true if the max kmer length is not -1
   */
  public boolean hasMaxKmerLen() {
    return (maxKmerLen != -1);
  }
  
  
  /**
   * Returns the max kmer length in this object
   * @return
   */
  public int getMaxKmerLen() {
    return maxKmerLen;
  }

  
  /**
   * Set the max kmer length for this object
   * @param maxKmerLen
   */
  public void setMaxKmerLen(int maxKmerLen) {
    this.maxKmerLen = maxKmerLen;
  }
  
  
  /**
   * Check whether this object has a DB model type
   * @return true if the db model type is not null
   */
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
   * Sets the db model type for this object ("FREQUENCY" or "MARKOV").
   * @param dbModelType
   */
  public void setDBModelType(String dbModelType) {
    if ((dbModelType == null) || dbModelType.equals(BackgroundModelLoader.MARKOV_TYPE_STRING) || dbModelType.equals(BackgroundModelLoader.FREQUENCY_TYPE_STRING)) {
      this.dbModelType = dbModelType;
    }
    else {
      throw new IllegalArgumentException("Invalid model type: " + dbModelType + 
          ". DB Model type must be either FREQUENCY or MARKOV or null");
    }    
  }


  /**
   * Checks whether this object has counts
   * @return
   */
  public boolean hasCounts() {
    return hasCounts;
  }

  
  /**
   * Set whether this object has counts
   * @param hasCounts
   */
  public void setHasCounts(boolean hasCounts) {
    this.hasCounts = hasCounts;
  }

  
  /**
   * Checks whether this object is equivalent to the specified other object
   * @param other
   * @return
   */
  public boolean equals(BackgroundModelMetadata other) {
    if ((other != null) && this.name.equals(other.getName()) && (this.maxKmerLen == other.getMaxKmerLen()) 
        && (((this.dbModelType == null) && (other.getDBModelType() == null)) || this.dbModelType.equals(other.dbModelType)) 
        && (this.modelID == other.getModelID() && this.mapID == other.getMapID())
        && (this.genomeID == other.getGenomeID())
        && (this.hasCounts == other.hasCounts())) {
      return true;
    }
    else {
      return false;
    }
  }
  
  
  /**
   * Formats this object as a String
   */
  public String toString() {
    return mapID + "\t" + genomeID + "\t" + modelID + "\t" + name + "\t" + maxKmerLen  + "\t" + dbModelType + "\t" + hasCounts; 
  }
}
