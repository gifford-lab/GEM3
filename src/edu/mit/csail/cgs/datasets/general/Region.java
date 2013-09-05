package edu.mit.csail.cgs.datasets.general;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Genome.ChromosomeInfo;
import edu.mit.csail.cgs.utils.Saveable;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.strings.StringUtils;

/**
 * A <code>Region</code> represents an interval along some chromosome in a
 * genome <br>
 * <i>Note</i>: We assume 0-based, inclusive coordinate.
 */

public class Region implements Comparable<Region>, Saveable {

  /**
   * regex for matching a region: newline-whitespace-one or more word chars or
   * digits-:-whitespace- one or more digits or commas-an optional m, M, k, or
   * K-whitespace- hyphen-whitespace-one or more digits or commas-an optional m,
   * M, k, or K- whitespace
   */
  public static final String COMPLETE_REGION_REG_EX = "^\\s*([\\w\\d]+):\\s*([,\\d]+[mMkK]?)\\s*-\\s*([,\\d]+[mMkK]?)\\s*";

  private static final Pattern REGION_PATTERN = Pattern.compile(COMPLETE_REGION_REG_EX);

  // private static final Matcher REGION_MATCHER = REGION_PATTERN.matcher("");

  /**
   * regex for matching a point: newline-whitespace-one or more word chars or
   * digits-:-whitespace- one or more digits or commas-an optional m, M, k, or
   * K-whitespace
   */
  public static final String POINT_REGION_REG_EX = "^\\s*([\\w\\d]+):\\s*([,\\d]+[mMkK]?)\\s*";

  private static final Pattern POINT_PATTERN = Pattern.compile(POINT_REGION_REG_EX);

  // private static final Matcher POINT_MATCHER = POINT_PATTERN.matcher("");

  /**
   * regex for matching a chromosome: newline-whitespace-one or more word chars
   * or digits-whitespace
   */
  public static final String CHROM_REGION_REG_EX = "^\\s*([\\w\\d]+)\\s*";

  private static final Pattern CHROM_PATTERN = Pattern.compile(CHROM_REGION_REG_EX);

  // private static final Matcher CHROM_MATCHER = CHROM_PATTERN.matcher("");


  /**
   * The genome that this region corresponds to
   */
  private Genome g;

  /**
   * The chromosome that this region corresponds to
   */
  private String chrom;

  /**
   * The start of this region (in bp)<br>
   * <i>Note</i>: The chromosome starts at 1
   */
  private int start;

  /**
   * The end of this region (in bp)<br>
   */
  private int end;


  public Region(Region copied) {
    g = copied.g;
    chrom = copied.chrom;
    start = copied.start;
    end = copied.end;
  }


  /**
   * Creates a new Region that encompasses first and second. Throws
   * IllegalArgumentException if <code>first</code> and <code>second</code> 
   * aren't on the same chromosome or they don't correspond to the same genome
   */
  public Region(Region first, Region second) {
    if (!first.g.equals(second.g)) {
      throw new IllegalArgumentException();
    }
    if (!first.chrom.equals(second.chrom)) {
      throw new IllegalArgumentException();
    }
    g = first.g;
    chrom = first.chrom;
    start = Math.min(first.start, second.start);
    end = Math.max(first.end, second.end);
  }


  /**
   * Creates a Region. Our codebase is probably not fully consistent about 0 or
   * 1 based coordinates and whether the start and end are inclusive or
   * exclusive.
   * 
   * @param g
   *          : the genome
   * @param c
   *          : the name of the chromosome
   * @param start
   *          : the start coordinate
   * @param end
   *          : the end coordinate
   */
  public Region(Genome g, String c, int start, int end) {
    if (start > end) {
      throw new IllegalArgumentException(String.format("Start > End for this region (chr %s) : %d > %d", c, start, end));
    }
    if (g == null) {
      throw new NullPointerException("Can't have a null genome");
    }
    if (c == null) {
      throw new NullPointerException("Can't have a null chromosome");
    }

    this.g = g;
    chrom = c;
    this.start = start;
    this.end = end;
  }


  public Region(Genome g, DataInputStream dis) throws IOException {
    this.g = g;
    chrom = dis.readUTF();
    start = dis.readInt();
    end = dis.readInt();
  }


  public void save(DataOutputStream dos) throws IOException {
    dos.writeUTF(chrom);
    dos.writeInt(start);
    dos.writeInt(end);
  }


  public Genome getGenome() {
    return g;
  }


  public String getChrom() {
    return chrom;
  }


  public int getStart() {
    return start;
  }


  public int getEnd() {
    return end;
  }


  public Region clone() {
    return new Region(g, chrom, start, end);
  }


  /**
   * @deprecated Some objects rely on Region being immutable. Change the start
   *             and end of the region.
   */
  protected void setStartEnd(int start, int end) {
    this.start = start;
    this.end = end;
    if (this.start < 1) {
      this.start = 0;
    }
    if (this.start > this.end) {
      this.start = this.end;
    }
  }


  /**
   * Subtracts the input Region from this Region and returns the remainder.<br>
   * In case, the two regions do not overlap, the region between them is
   * returned.
   * 
   * @throws IllegalArgumentException
   *           if <code>r</code> is not on the same chromosome as this or the
   *           two regions correspond to different genomes
   */
  public Collection<Region> getSubtractionFragments(Region r) {
    if (!g.equals(r.g)) {
      throw new IllegalArgumentException();
    }
    if (!chrom.equals(r.chrom)) {
      throw new IllegalArgumentException();
    }

    LinkedList<Region> lst = new LinkedList<Region>();

    if (start != r.start || end != r.end) {
      if (overlaps(r)) {
        lst.addLast(new Region(g, chrom, Math.min(start, r.start), Math.max(start, r.start) - 1));
        lst.addLast(new Region(g, chrom, Math.min(end, r.end) + 1, Math.max(end, r.end)));
      }
      else {
        if (before(r)) {
          lst.addLast(new Region(g, chrom, end + 1, r.start - 1));
        }
        else {
          lst.addLast(new Region(g, chrom, r.end + 1, start - 1));
        }
      }
    }

    return lst;
  }


  /**
   * Returns the distance between <code>r</code> and <code>this</code>.
   * Overlapping Regions have distance zero.
   * 
   * @throws IllegalArgumentException
   *           if <code>r</code> is not on the same chromosome as
   *           <code>this</code>
   */
  public int distance(Region r) {
    if (!chrom.equals(r.chrom)) {
      throw new IllegalArgumentException(r.chrom + " is not " + chrom);
    }
    if (overlaps(r)) {
      return 0;
    }
    if (before(r)) {
      return r.start - end;
    }
    if (after(r)) {
      return start - r.end;
    }
    return 0;
  }


  /**
   * Returns the distance between <code>p</code> and <code>this</code>
   * 
   * @throws IllegalArgumentException
   *           if <code>p</code> is not on the same chromosome as
   *           <code>this</code>
   */
  public int distance(Point p) {
    if (!chrom.equals(p.getChrom())) {
      throw new IllegalArgumentException(p.getChrom());
    }
    if (p.getLocation() < start) {
      return start - p.getLocation();
    }
    if (p.getLocation() > end) {
      return p.getLocation() - end;
    }
    return 0;
  }


  public Region getOverlap(Region r) {
    if (!overlaps(r)) {
      return null;
    }
    int ns = Math.max(start, r.start);
    int ne = Math.min(end, r.end);
    return new Region(g, chrom, ns, ne);
  }


  /**
   * Returns the number of bp overlap between <code>r</code> and
   * <code>this</code>.
   */
  public int getOverlapSize(Region r) {
    if (!overlaps(r)) {
      return 0;
    }
    int ns = Math.max(start, r.start), ne = Math.min(end, r.end);
    return ne - ns + 1;
  }


  /**
   * Returns the width (or size) of this region. This method assumes that both
   * <code>start</code> and <code>end</code> are inclusive
   */
  public int getWidth() {
    return end - start + 1;
  }


  /**
   * Returns true iff <code>this</code> comes before <code>r</code> along the
   * chromosome.
   * 
   * @throws IllegalArgumentException
   *           if <code>r</code> is not on the same chromosome as
   *           <code>this</code> or corresponds to a different genome
   */
  public boolean before(Region r) {
    if (!g.equals(r.g)) {
      throw new IllegalArgumentException();
    }
    if (!chrom.equals(r.chrom)) {
      throw new IllegalArgumentException();
    }
    return end < r.start;
  }


  /**
   * Returns true iff <code>this</code> comes after <code>r</code> along the
   * chromosome.
   * 
   * @throws IllegalArgumentException
   *           if <code>r</code> is not on the same chromosome as
   *           <code>this</code> or corresponds to a different genome
   */
  public boolean after(Region r) {
    if (!g.equals(r.g)) {
      throw new IllegalArgumentException();
    }
    if (!chrom.equals(r.chrom)) {
      throw new IllegalArgumentException();
    }
    return start > r.end;
  }


  /**
   * Returns a new Region that is <code>this</code> expanded by
   * <code>dstart</code> bases to the left and <code>dend</code> bases to the
   * right.<br>
   * <i>Note</i>: If any of the arguments exceeds the corresponding ends of the
   * chromosome, the expansion is limited to that end.
   */
  public Region expand(int dstart, int dend) {
    int chromLength = g.getChromLength(chrom);
    int ns = start - dstart;
    int ne = end + dend;
    if (ns < 0) {
      ns = 0;
    }
    if (ne > chromLength - 1) {
      ne = chromLength - 1;
      if (ns > ne) {
        ns = ne;
      }
    }
    return new Region(g, chrom, ns, ne);
  }
  
  /**
   * Returns a new Region with the same midpoint as <code>this</code>
   * but with width <code>width</code>
   * @param width
   * @return
   */
  public Region resize(int width) {
	  return this.getMidpoint().expand(width/2);
  }


  /**
   * Checks for full or partial overlapping between <code>this</code> and
   * <code>r</code> region.
   * 
   * @param r
   * @return
   */
  public boolean overlaps(Region r) {
    if (!chrom.equals(r.chrom)) {
      return false;
    }
    if (start <= r.start && end >= r.start) {
      return true;
    }
    if (r.start <= start && r.end >= start) {
      return true;
    }
    return false;
  }


  /**
   * Checks for full or partial overlapping between <code>this</code> region and
   * the region that corresponds between <code>otherstart</code> and
   * <code>otherend</code>.
   * 
   * @param otherstart
   * @param otherend
   * @return
   */
  public boolean overlaps(int otherstart, int otherend) {
    if (start <= otherstart && end >= otherstart) {
      return true;
    }
    if (otherstart <= start && otherend >= start) {
      return true;
    }
    return false;
  }


  /**
   * Checks whether <tt>r</tt> region is fully contained in <tt>this</tt>
   * region.
   * 
   * @param r
   * @return
   */
  public boolean contains(Region r) {
    if (!chrom.equals(r.chrom)) {
      return false;
    }
    return start <= r.start && end >= r.end;
  }


  /**
   * Checks whether point <tt>p</tt> is contained in <tt>this</tt> region.
   * 
   * @param p
   * @return
   */
  public boolean contains(Point p) {
    if (!chrom.equals(p.getChrom())) {
      return false;
    }
    return start <= p.getLocation() && end >= p.getLocation();
  }


  public Point startPoint() {
    return new Point(getGenome(), getChrom(), getStart());
  }


  public Point endPoint() {
    return new Point(getGenome(), getChrom(), getEnd());
  }


  /**
   * Checks whether region <tt>r</tt> is the same as <tt>this</tt> region.
   * 
   * @param r
   * @return
   */
  public boolean matchesRegion(Region r) {
    return g.equals(r.getGenome()) && chrom.equals(r.getChrom()) && (start == r.getStart()) && (end == r.getEnd());
  }


  /**
   * Returns the <i>super-region</i> that contains both regions <tt>this</tt>
   * and <tt>o</tt>.
   * 
   * @param o
   * @return
   */
  public Region combine(Region o) {
    if (!getChrom().equals(o.getChrom())) {
      throw new IllegalArgumentException(getChrom() + " != " + o.getChrom());
    }
    int ns = Math.min(getStart(), o.getStart());
    int ne = Math.max(getEnd(), o.getEnd());
    return new Region(getGenome(), getChrom(), ns, ne);
  }


  /**
   * Returns the mid-point (as Point) that corresponds to this region.
   * 
   * @return
   */
  public Point getMidpoint() {
    int middle = (start + end) / 2;
    return new Point(g, chrom, middle);
  }


  /**
   * regionString() should <b>not</b> be overridden in any subclasses of Region
   * -- it is meant to always return a string with the chromosome, start, and
   * stop (and only that information).
   * 
   * @return A string with the format, "chrom:start-end"
   */
  public String regionString() {
    return String.format("%s:%d-%d", chrom, start, end);
  }


  /**
   * 
   * @param input
   * @return
   */
  public static boolean isValidCompleteRegionString(String input) {
    String trimmed = input.trim();
    Matcher regionMatcher = REGION_PATTERN.matcher(trimmed);
    return regionMatcher.matches();
  }


  /**
   * 
   * @param input
   * @return
   */
  public static boolean isValidPointRegionString(String input) {
    String trimmed = input.trim();
    Matcher pointMatcher = POINT_PATTERN.matcher(trimmed);
    return pointMatcher.matches();
  }


  /**
   * 
   * @param input
   * @return
   */
  public static boolean isValidChromRegionString(String input) {
    String trimmed = input.trim();
    Matcher chromMatcher = CHROM_PATTERN.matcher(trimmed);
    return chromMatcher.matches();
  }


  /**
   * 
   * @param input
   * @return
   */
  public static String chromFromLocationString(String input) {
    String trimmed = input.trim();
    Matcher chromMatcher = CHROM_PATTERN.matcher(trimmed);
    String chromStr = null;
    if (chromMatcher.find()) {
      chromStr = chromMatcher.group(1);
      if (chromStr.startsWith("chr")) {
        chromStr = chromStr.substring(3, chromStr.length());
      }
    }
    return chromStr;
  }


  /**
   * 
   * @param input
   * @return
   */
  public static int startFromLocationString(String input) {
    String trimmed = input.trim();
    Matcher regionMatcher = REGION_PATTERN.matcher(trimmed);
    int start = -1;
    if (regionMatcher.matches()) {
      String startStr = regionMatcher.group(2);
      String endStr = regionMatcher.group(3);
      startStr = startStr.replaceAll(",", "");
      endStr = endStr.replaceAll(",", "");
      start = Math.min(stringToNum(startStr), stringToNum(endStr));
    }
    return start;
  }


  /**
   * 
   * @param input
   * @return
   */
  public static synchronized int stopFromLocationString(String input) {   
    String trimmed = input.trim();
    Matcher regionMatcher = REGION_PATTERN.matcher(trimmed);
    int end = -1;
    if (regionMatcher.matches()) {
      String startStr = regionMatcher.group(2);
      String endStr = regionMatcher.group(3);
      startStr = startStr.replaceAll(",", "");
      endStr = endStr.replaceAll(",", "");
      end = Math.max(stringToNum(startStr), stringToNum(endStr));
    }
    return end;    
  }


  /**
   * getLocationString() returns, by default, the same as regionString().
   * However, this method can be overridden by subclasses to add additional
   * information to the returned String. Any additional information added by a
   * subclass should be separated by a ":". For instance, a StrandedRegion might
   * return "chrom:start-stop:strand".
   * 
   * @return A string representation of the Region's genomic location
   */
  public String getLocationString() {
    return regionString();
  }


  /**
   * By default, returns the result of getLocationString(). Can be modified in
   * any appropriate way by a subclass.
   * 
   * @return A string represention of the object.
   */
  public String toString() {
    return getLocationString();
  }

  /**
   * parses the input String into a Region. The Region can either be strict
   * coordinates (eg, chrom/start/stop). Understands abbreviates in the
   * coordinates such as k and m. <br>
   * The method accepts either the form: <blockquote>
   * 
   * <pre>
   * chromosome:start-end:[strand] ([] stands for optional) - <i>strictly specified region</i>  
   * <i>or</i>
   * chromosome:start:[strand] - <i>region from start to start+1</i>
   * <i>or</i>
   * chromosome:[strand] - <i>whole chromosome</i>
   * </pre>
   * 
   * </blockquote>
   */
  public static Region fromString(Genome genome, String input) throws NumberFormatException {
    String pieces[] = input.split(":");
    char strand = ' ';
    if (pieces.length == 3 && pieces[2].length() == 1) {
      strand = pieces[2].charAt(0);
      input = pieces[0] + ":" + pieces[1];
    }

    String trimmed = input.trim();

    Matcher regionmatcher = REGION_PATTERN.matcher(trimmed);
    Matcher pointmatcher = POINT_PATTERN.matcher(trimmed);
    Matcher chrommatcher = CHROM_PATTERN.matcher(trimmed);

    Region output = null;
    if (regionmatcher.find()) {
      if (regionmatcher.groupCount() != 3) {
        return null;
      }
      String chromStr = regionmatcher.group(1);
      String startStr = regionmatcher.group(2);
      String endStr = regionmatcher.group(3);
      if (chromStr.startsWith("chr")) {
        chromStr = chromStr.substring(3, chromStr.length());
      }
      startStr = startStr.replaceAll(",", "");
      endStr = endStr.replaceAll(",", "");
      int start = Math.min(stringToNum(startStr), stringToNum(endStr));
      int end = Math.max(stringToNum(startStr), stringToNum(endStr));
      output = new Region(genome, chromStr, start, end);
    }
    else if (pointmatcher.find()) {
      if (pointmatcher.groupCount() != 2) {
        return null;
      }
      String chromStr = pointmatcher.group(1);
      String startStr = pointmatcher.group(2);
      if (chromStr.startsWith("chr")) {
        chromStr = chromStr.substring(3, chromStr.length());
      }
      startStr = startStr.replaceAll(",", "");
      output = new Region(genome, chromStr, stringToNum(startStr), stringToNum(startStr) + 1);
    }
    else if (chrommatcher.matches()) {
      String chromStr = chrommatcher.group(1);
      if (chromStr.startsWith("chr")) {
        chromStr = chromStr.substring(3, chromStr.length());
      }
      ChromosomeInfo info = genome.getChrom(chromStr);
      if (info != null) {
        int length = info.getLength();
        output = new Region(genome, chromStr, 0, length - 1);
      }
    }
    if (output != null && strand != ' ') {
      output = new StrandedRegion(output, strand);
    }

    return output;
  }


  /**
   * Parses integers but understand k and m abbreviations for kilo and mega.
   */
  public static int stringToNum(String s) {
    if (s.matches(".*[kK]$")) {
      return 1000 * Integer.parseInt(s.substring(0, s.length() - 1));
    }
    if (s.matches(".*[mM]$")) {
      return 1000000 * Integer.parseInt(s.substring(0, s.length() - 1));
    }
    return Integer.parseInt(s);
  }
  
  public static boolean[] overlap(List<Region> sourceRegions, List<Region> targetRegions) {
	  return overlap(sourceRegions.toArray(new Region[0]), targetRegions.toArray(new Region[0]));
  }//end of overlap method
  
  /**
   * This method finds all possible overlaps between <tt>N source</tt> and <tt>M target</tt>
   * regions. 																				<br>
   * It runs in <tt>O(N log(M)) + O(M log(M))</tt> time compared to the naive (nested loop) search 
   * between all the elements of <tt>source</tt> and <tt>target</tt> regions which runs in O(N M) time.
   * @param sourceRegions <tt>N</tt> source regions where we want to find overlaps between these and the <tt>target</tt> regions 
   * @param targetRegions <tt>M</tt> target regions
   * @return <tt>true</tt> in all positions where there is an overlap between a <tt>source</tt> and
   * a <tt>target</tt> region. <tt>false</tt>, otherwise.
   */
  public static boolean[] overlap(Region[] sourceRegions, Region[] targetRegions) {
	  boolean[] overlaps = new boolean[sourceRegions.length];
	  
	  if(sourceRegions.length == 0 || targetRegions.length == 0) { return overlaps; }
	  
	  Arrays.sort(targetRegions);
	  int[] target_start = new int[targetRegions.length];
	  for(int i = 0; i < target_start.length; i++) { target_start[i] = targetRegions[i].getStart(); }
	  int[] target_end = new int[targetRegions.length];
	  for(int i = 0; i < target_end.length; i++) { target_end[i] = targetRegions[i].getEnd(); }
	  int[] idx = StatUtil.findSort(target_end);
	  
	  for(Region sourceReg:sourceRegions) {
		  boolean sourceRegOverlaps = false;
		  int start_idx = Arrays.binarySearch(target_end, sourceReg.getStart());
		  int end_idx   = Arrays.binarySearch(target_start, sourceReg.getEnd());
		  if( start_idx < 0 ) { start_idx = Math.min(-start_idx - 1, target_start.length-1); }
		  if( end_idx < 0 )   { end_idx   = Math.min(-end_idx - 1  , target_start.length-1); }
		  
		  for(int i = idx[start_idx]; i <= end_idx; i++) {
			  if(targetRegions[i].overlaps(sourceReg)) {
				  sourceRegOverlaps = true;
				  break;	  
			  }
		  }
		  
		  if(!sourceRegOverlaps && (end_idx - idx[start_idx] + 1 != target_start.length)) {
			  for(int i = Math.min(idx[start_idx]+1, target_start.length-1); i > 0; i--) {
				  if(targetRegions[i].overlaps(targetRegions[i-1])) {
					  if(targetRegions[i-1].overlaps(sourceReg)){
						  sourceRegOverlaps = true;
						  break;
					  }
				  }
					else
						break;
			  }
		  }
	  }
	  
	  return overlaps;
  }//end of overlap method


  public boolean equals(Object o) {
    if (o instanceof Region) {
      Region r = (Region) o;
      return matchesRegion(r);
    }
    else {
      return false;
    }
  }


  public int hashCode() {
    int code = 17;
    code += g.hashCode();
    code *= 37;
    code += chrom.hashCode();
    code *= 37;
    code += start;
    code *= 37;
    code += end;
    code *= 37;
    return code;
  }


  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  public int compareTo(Region r) {
    if (!chrom.equals(r.chrom)) {
      return chrom.compareTo(r.chrom);
    }
    if (start < r.start) {
      return -1;
    }
    if (start > r.start) {
      return 1;
    }
    if (end < r.end) {
      return -1;
    }
    if (end > r.end) {
      return 1;
    }
    return 0;
  }
  
  /**
   * Given a list of regions, filter out all the regions that overlaps with any reference regions, return the non-overlap regions
   * @param rs
   * @param refs
   * @return
   */
  public static ArrayList<Region> filterOverlapRegions(ArrayList<Region> rs, ArrayList<Region> refs){
	    ArrayList<Region> results = new ArrayList<Region>();
		Map<String, ArrayList<Region>> chr2rs = new HashMap<String, ArrayList<Region>>();
		// group regions by chrom
		for(Region r:rs) {
			String chrom = r.getChrom();
			if(!chr2rs.containsKey(chrom))
				chr2rs.put(chrom, new ArrayList<Region>());
			chr2rs.get(chrom).add(r);
		}
		Map<String, ArrayList<Region>> chr2refs = new HashMap<String, ArrayList<Region>>();
		for(Region r:refs) {
			String chrom = r.getChrom();
			if(!chr2rs.containsKey(chrom))
				continue;
			if(!chr2refs.containsKey(chrom))
				chr2refs.put(chrom, new ArrayList<Region>());
			chr2refs.get(chrom).add(r);
		}	
		// for each chrom
		for (String chr : chr2rs.keySet()){
			ArrayList<Region> rsc = chr2rs.get(chr);
			if (!chr2refs.containsKey(chr)){
				results.addAll(rsc);
				continue;
			}
			ArrayList<Region> refsc = chr2refs.get(chr);
			Collections.sort(rsc);
			Collections.sort(refsc);
			int rId=0;
			int refId=0;
			while(rId<rsc.size()){
				Region r = rsc.get(rId);
				Region ref = refsc.get(refId);
				if (r.overlaps(ref)){
					rId++;
					continue;
				}
				else{ 
					if(r.before(ref)){
						results.add(r);
						rId++;
						continue;
					}
					else if(r.after(ref)){
						refId++;
						if (refId==refsc.size()){		// end of ref regions
							for (int i=rId+1;i<rsc.size();i++)
								results.add(rsc.get(i));
							break;
						}
					}
				}
			}
		}
		
	    return results;
  }
	/**
	 * Merge the overlapped regions<br>
	 * The regions will be sorted
	 * @param regions
	 * @return
	 */
	public static ArrayList<Region> mergeRegions(ArrayList<Region> regions){
		if (regions.isEmpty())
			return regions;
		ArrayList<Region> mergedRegions = new ArrayList<Region>();
		Collections.sort(regions);
		Region previous = regions.get(0);

		for (int i = 1; i < regions.size(); i++) {
          Region region = regions.get(i);
			// if overlaps with previous region, combine the regions
			if (previous.overlaps(region)){
				previous = previous.combine(region);
			} else{
				mergedRegions.add(previous);
				previous = region;
			}
		}
		mergedRegions.add(previous);
		mergedRegions.trimToSize();
		return mergedRegions;
	}//end of mergeRegions method

}
