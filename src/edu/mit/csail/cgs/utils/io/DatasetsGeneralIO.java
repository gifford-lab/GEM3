/**
 * 
 */
package edu.mit.csail.cgs.utils.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author rca
 *
 */
public class DatasetsGeneralIO {

  /**
   * Read regions from a flat file with one region listed per line
   * 
   * @param filename
   * @return
   */
  public static Vector<Region> readRegionsFromFile(Genome genome, String filename) throws IOException {
    
    Vector<Region> regions = new Vector<Region>();
        
    String[] lines = LineByLineFileReader.readFile(filename, LineByLineFileReader.DEFAULT_COMMENT_PREFIXES);
    for (String currLine : lines) {
      Region region = Region.fromString(genome, currLine);
      regions.add(region);
    }
    
    return regions;
  }
  
  /**
   * Written to read regions from a flat file with one region listed per line
   * 
   * @param filename
   * @return
   */
  public static Vector<Point> readPointsFromFile(Genome genome, String filename) throws IOException {
    Vector<Point> points = new Vector<Point>();
    
    String[] lines = LineByLineFileReader.readFile(filename, LineByLineFileReader.DEFAULT_COMMENT_PREFIXES);
    
    for (String currLine : lines) {
      Point point = Point.fromString(genome, currLine);
      points.add(point);
    }
    
    return points;
  }
}
