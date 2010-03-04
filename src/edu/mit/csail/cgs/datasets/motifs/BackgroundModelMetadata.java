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

  //fields from the background model table
  public final int modelID;
  public final String name;
  public final int kmerlen;
  public final String modelType;
  
  /**
   * these fields may be empty if the instance justs represents a row from the
   * background model table
   */
  public final Integer mapID;
  public final Integer genomeID;

  
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType) {
    this(modelID, name, kmerlen, modelType, null, null);
  }
  
  
  public BackgroundModelMetadata(int modelID, String name, int kmerlen, String modelType, Integer mapID, Integer genomeID) {
    this.modelID = modelID;
    this.name = name;
    this.kmerlen = kmerlen;
    this.modelType = modelType;
    this.mapID = mapID;
    this.genomeID = genomeID;
  }
}
