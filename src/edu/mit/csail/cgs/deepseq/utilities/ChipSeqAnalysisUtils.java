/**
 * 
 */
package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.ewok.verbs.RegionParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.utils.Pair;

/**
 * @author rca
 * 
 */
public class ChipSeqAnalysisUtils {


  private static final Logger logger = Logger.getLogger(ChipSeqAnalysisUtils.class);


  public static double matrixMax(DoubleMatrix1D mat) {
    return mat.aggregate(Functions.max, Functions.identity);
  }
  
  /**
   * Loads a set of regions from a MACS peak file. <br>
   * For the proper function of the method, regions contained in the MACS file
   * should be sorted, something done inside the method's body.
   * 
   * @param peaks
   *          The <tt>List</tt> created by a MACS file
   * @return
   */
  public static ArrayList<Region> selectEnrichedMACSRegions(int modelRange, List<MACSPeakRegion> peaks) {
    ArrayList<Region> regions = new ArrayList<Region>();
    for (MACSPeakRegion p : peaks) {
      // filter towers
      if (p.getTags() > 1000 && p.getFold_enrichment() < 5) {
        // TODO: should compare with WCE to decide it is really a tower
        continue;
      }
      regions.add(p);
    }
    return mergeRegions(modelRange, regions);
  }


  /**
   * Loads a set of regions (Format: Chrom:start-end). <br>
   * For the proper function of the method, regions contained in the file should
   * be sorted something done inside the method's body.
   * 
   * @param fname
   *          the file containing the regions
   * @return
   */
  public static ArrayList<Region> loadRegionsFromFile(Genome gen, String fname,
      boolean expandedRegion, int modelRange) {
    logger.info("loading regions ... ");
    ArrayList<Region> rset = new ArrayList<Region>();
    try {
      File rFile = new File(fname);
      if (!rFile.isFile()) {
        logger.fatal("Invalid file name");
        System.exit(1);
      }
      BufferedReader reader = new BufferedReader(new FileReader(rFile));
      String line;
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        if (!line.equals("")) {
          String[] words = line.split("\\s+");
          RegionParser parser = new RegionParser(gen);
          Region r = parser.execute(words[0]);
          rset.add(r);
        }
      }
      reader.close();
    }
    catch (FileNotFoundException fnfex) {
      logger.error(fnfex.getMessage(), fnfex);
    }
    catch (IOException ioex) {
      logger.error(ioex.getMessage(), ioex);
    }

    // if regions are from previous round of analysis
    if (expandedRegion) {
      return rset;
    }
    else { // regions are from other sources (i.e. StatPeakFinder)
      return ChipSeqAnalysisUtils.mergeRegions(modelRange, rset);
    }
  }


  /**
   * Expand selected regions to include areas with reads that could be accounted
   * for by binding events at the edges of the regions.
   * 
   * @param model
   *          BindingModel object specifying the fragment length distribution
   * @param regions
   *          the regions to expand and merge
   * @return
   */
  public static ArrayList<Region> mergeRegions(int modelRange, ArrayList<Region> regions) {
    ArrayList<Region> mergedRegions = new ArrayList<Region>();
    Region previous = regions.get(0).expand(modelRange, modelRange);
    mergedRegions.add(previous);

    Collections.sort(regions);
    for (Region region : regions) {
      // expand on both side to leave enough space
      // to include every potential read
      Region r = region.expand(modelRange, modelRange);
      // if overlaps with previous region, combine the regions
      // note that region has non-expanded original width
      // overlaps if region and previous region are within getRange()
      if (previous.overlaps(region)) {
        mergedRegions.remove(previous);
        mergedRegions.add(previous.combine(r));
      }
      else {
        mergedRegions.add(r);
      }
      previous = r;
    }
    return mergedRegions;
  } // end of mergeRegions method


  /*
   * filter out duplicate reads (potential tower, needles, but could be real
   * reads) assuming the reads are sorted
   */
  public static List<ReadHit> filterDuplicateReads(List<ReadHit> reads, int maxDuplicates) {
    ReadHit currentHit = reads.get(0);
    int count = 0;
    List<ReadHit> filteredReads = new ArrayList<ReadHit>();
    for (ReadHit hit : reads) {
      // if read from a new position
      if (!(currentHit.getStart() == hit.getStart())) {
        currentHit = hit;
        count = 0;
        filteredReads.add(hit);
      }
      else {// if duplicate
        count++;
        if (count <= maxDuplicates) {
          filteredReads.add(hit);
        }
      }
    }
    return filteredReads;
  }

  
  /**
   * Tests if a region <tt>r</tt> being represented by the <tt>ReadHits signals</tt>
   * is a tower or not. 
   * @param r  Region to be tested if it is a tower
   * @param signals Readhits corresponding to region <tt>r</tt>
   * @return <tt>true</tt> if the regions is a tower, <tt>false</tt> otherwise
   */
  public static boolean isTower(Region r, ArrayList<List<ReadHit>> signals, double totalSigCount, int towerThreshold) {
    boolean isTower = false;
    if (totalSigCount > towerThreshold) {
      // HashMap storing
      HashMap<Integer, Integer> readCoverage = new HashMap<Integer, Integer>();
      for (List<ReadHit> condSignal : signals) {
        for (ReadHit h : condSignal) {
          int read = h.getFivePrime();
          if (readCoverage.containsKey(read)) {
            readCoverage.put(read, readCoverage.get(read) + 1);
          }
          else {
            readCoverage.put(read, 1);
          }
        }
      }

      int maxCoverage = 0;
      for (int read : readCoverage.keySet()) {
        maxCoverage = Math.max(maxCoverage, readCoverage.get(read));
      }

      if (maxCoverage >= 10) {
        logger.debug("Suspicious region " + r.getLocationString() + "\tMaxCoverage=" + maxCoverage
            + "\tTotalReads=" + totalSigCount + "\tSkipped ...");

        // TODO: store the tower regions
        isTower = true;
      }
    }
    return isTower;
  }// end of isTower method

  
  /**
   * Alternate signature for single condition analysis
   * @param r
   * @param signalReads
   * @param totalSigCount
   * @param towerThreshold
   * @return
   */
  public static boolean isTower(Region r, List<ReadHit> signalReads, double totalSigCount, int towerThreshold) {
    ArrayList<List<ReadHit>> signals = new ArrayList<List<ReadHit>>();
    signals.add(signalReads);
    return isTower(r, signals, totalSigCount, towerThreshold);
  }
  
  /**
   * Write out the list of regions to a file
   * 
   * @param restrictRegions
   *          the list of regions
   * @param conditionNames
   *          names of the experiment conditions
   * @param isInitial
   *          @todo ???
   */
  public static void saveRestrictRegions(List<Region> restrictRegions, List<String> conditionNames, boolean isInitial) {

    // save the list of regions to file
    try {
      StringBuilder fileName = new StringBuilder("PEAK_");
      if (conditionNames != null) {
        for (String cond : conditionNames) {
          fileName.append(cond).append("_");
        }
      }
      if (isInitial) {
        fileName.append("Init_");
      }
      else {
        fileName.append("Refined_");
      }

      fileName.append("Regions.txt");
      FileWriter fw = new FileWriter(fileName.toString());

      StringBuilder txt = new StringBuilder();
      for (Region r : restrictRegions) {
        txt.append(r.toString()).append("\n");
      }
      fw.write(txt.toString());
      fw.close();
    }
    catch (IOException ioex) {
      logger.error(ioex.getMessage(), ioex);
    }
  }


  /**
   * Alternate signature for single condition experiments
   * 
   * @param restrictRegions
   * @param isInitial
   */
  public static void saveRestrictRegions(List<Region> restrictRegions, boolean isInitial) {
    ChipSeqAnalysisUtils.saveRestrictRegions(restrictRegions, null, isInitial);
  }
}
