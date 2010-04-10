package edu.mit.csail.cgs.datasets.general;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * A <code>Point</code> represents a single base position along some chromosome
 * in a genome <br>
 * <i>Note</i>: We assume 0-based, inclusive coordinate.
 */
public class Point implements Comparable<Point> {

  private static final Pattern POINT_PATTERN = Pattern.compile(Region.POINT_REGION_REG_EX);
  
  /**
   * The genome that this point corresponds to
   */
  private Genome g;

  /**
   * The chromosome that this point corresponds to
   */
  private String chrom;

  /**
   * The location that this point corresponds to
   */
  private int location;


  public Point(Genome g, String c, int position) {
    this.g = g;
    chrom = c;
    this.location = position;
  }


  public Point(Point copied) {
    this.g = copied.getGenome();
    this.chrom = copied.getChrom();
    this.location = copied.getLocation();
  }


  public Genome getGenome() {
    return g;
  }


  public String getChrom() {
    return chrom;
  }


  public int getLocation() {
    return location;
  }


  public Point clone() {
    return new Point(g, chrom, location);
  }


  /**
   * The <code>expand</code> method returns the region which is around
   * <code>distance</code> genomic coordinates (in bp) from this point.<br>
   * <i>Note</i>: The expansion happens in both directions.<br>
   * However, if <code>distance</code> exceeds any of chromosome ends, the
   * expansion takes place until this end.<br>
   * <b>In this method, indexing starts from 0</b>
   * 
   * @param distance
   *          The length of bidirectional expansion
   * @return The region (as a Region) that corresponds to this expansion
   */
  public Region expand(int distance) {
	  int ns = Math.max(0, location - distance);
	  int ne = Math.min(location + distance, g.getChromLength(chrom)-1);
	  return new Region(g, chrom, ns, ne);
  }


  public String toString() {
    return getLocationString();
  }


  public String getLocationString() {
    return chrom + ":" + location;
  }


  public int hashCode() {
    int code = 17;
    code += g.hashCode();
    code *= 37;
    code += chrom.hashCode();
    code *= 37;
    code += location;
    code *= 37;
    return code;
  }


  public boolean equals(Object o) {
    if (o instanceof Point) {
      Point r = (Point) o;
      return g.equals(r.getGenome()) && chrom.equals(r.getChrom()) && (location == r.getLocation());
    }
    else {
      return false;
    }
  }


  /**
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  public int compareTo(Point p) {
    if (!chrom.equals(p.chrom)) {
      return chrom.compareTo(p.chrom);
    }
    if (location < p.location) {
      return -1;
    }
    if (location > p.location) {
      return 1;
    }
    return 0;
  }


  // Distance to another point
  public int distance(Point p) {
    if (!chrom.equals(p.getChrom())) {
      throw new IllegalArgumentException(p.getChrom());
    }

    return Math.abs(location - p.getLocation());
  }


  // Offset relative to another point (offset = this-p)
  public int offset(Point p) {
    if (!chrom.equals(p.getChrom())) {
      throw new IllegalArgumentException(p.getChrom());
    }

    return location - p.getLocation();
  }
  
  /**
   * parses the input String into a Point. Understands abbreviates in the
   * coordinates such as k and m. <br>
   * The method accepts the form: <blockquote>
   * 
   * <pre>
   * chromosome:start
   * </pre>
   * 
   * </blockquote>
   */
  public static Point fromString(Genome genome, String input) throws NumberFormatException {
    String pieces[] = input.split(":");
    char strand = ' ';
    if (pieces.length == 3 && pieces[2].length() == 1) {
      strand = pieces[2].charAt(0);
      input = pieces[0] + ":" + pieces[1];
    }

    String trimmed = input.trim();
    
    Matcher pointmatcher = POINT_PATTERN.matcher(trimmed);

    Point output = null;
    
    if (pointmatcher.find()) {
      if (pointmatcher.groupCount() != 2) {
        return null;
      }
      String chromStr = pointmatcher.group(1);
      String locStr = pointmatcher.group(2);
      if (chromStr.startsWith("chr")) {
        chromStr = chromStr.substring(3, chromStr.length());
      }
      locStr = locStr.replaceAll(",", "");
      output = new Point(genome, chromStr, Region.stringToNum(locStr));
    }

    return output;
  }

}
