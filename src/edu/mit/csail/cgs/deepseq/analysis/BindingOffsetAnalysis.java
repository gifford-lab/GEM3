package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.io.LineByLineFileWriter;

public class BindingOffsetAnalysis {

  private static final int MAX_DIST = 16;
  private static final int MOTIF_DISTANCE = 16;  
  private static final int MIN_DIST = 100;


  public BindingOffsetAnalysis(String inputFilenameOne, String inputFilenameTwo, String outFile, 
      String tfOneName, String tfTwoName, Genome genome, Organism org, String motifString, String motifVersion, double motifThreshold, 
      boolean printZeroDist, boolean requireMotifContainment) throws IOException, NotFoundException {


    List<GPSPeak> gpsPeaksOne;
    List<GPSPeak> gpsPeaksTwo;

    ArrayList<Point> gpsPeakLocsOne;
    TreeSet<Point> gpsPeakLocsTwo;

    TreeMap<Point, GPSPeak> peaksByPointOne = new TreeMap<Point, GPSPeak>();
    TreeMap<Point, GPSPeak> peaksByPointTwo = new TreeMap<Point, GPSPeak>();

    TreeMap<Point, Pair<Point, Point>> pairsByPoint = null; //= new TreeMap<Point, Pair<Point, Point>>();
    TreeMap<Integer, List<Pair<Point, Point>>> pairsByDist = null; //= new TreeMap<Integer, List<Pair<Point, Point>>>();
    TreeMap<Integer, List<GPSPeak>> pointsBySize = new TreeMap<Integer, List<GPSPeak>>();


    WeightMatrixScorer scorer;
    WeightMatrix motif = null;	  

    /**
     * 
     */

    int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), motifString, motifVersion);
    motif = WeightMatrix.getWeightMatrix(wmid);
    scorer = new WeightMatrixScorer(motif);


    File gpsFileOne = new File(inputFilenameOne);
    gpsPeaksOne = GPSParser.parseGPSOutput(gpsFileOne.getAbsolutePath(), genome);

    File gpsFileTwo = new File(inputFilenameTwo);
    gpsPeaksTwo = GPSParser.parseGPSOutput(gpsFileTwo.getAbsolutePath(), genome);

    gpsPeakLocsOne = new ArrayList<Point>();
    gpsPeakLocsTwo = new TreeSet<Point>();

    for (GPSPeak peak : gpsPeaksOne) {
      Point peakPoint = new Point(genome, peak.getChrom(), peak.getLocation());
      gpsPeakLocsOne.add(peakPoint);
      peaksByPointOne.put(peakPoint, peak);

      int size = (int)Math.round(peak.getStrength());
      if (!pointsBySize.containsKey(size)) {
        pointsBySize.put(size, new ArrayList<GPSPeak>());
      }
      pointsBySize.get(size).add(peak);
    }

    for (GPSPeak peak : gpsPeaksTwo) {
      Point peakPoint = new Point(genome, peak.getChrom(), peak.getLocation());
      gpsPeakLocsTwo.add(peakPoint);
      peaksByPointTwo.put(peakPoint, peak);

      int size = (int)Math.round(peak.getStrength());
      if (!pointsBySize.containsKey(size)) {
        pointsBySize.put(size, new ArrayList<GPSPeak>());
      }
      pointsBySize.get(size).add(peak);
    }

    //		Pair<TreeMap<Point, Pair<Point, Point>>, TreeMap<Integer, List<Pair<Point, Point>>>> matchedPeaks = this.findMatchedPeaks(gpsPeakLocsOne, gpsPeakLocsTwo);
    //		pairsByPoint = matchedPeaks.car();
    //		pairsByDist = matchedPeaks.cdr();
    //		
    //		this.writeMatchedPeaks(outFile, tfOneName, tfTwoName, printZeroDist, requireMotifContainment, scorer, motif, motifThreshold, pairsByPoint, pairsByDist, peaksByPointOne, peaksByPointTwo);
    //		

    Pair<TreeMap<Point, Pair<Integer, Pair<Character, Double>>>, 
    TreeMap<Point, Pair<Integer, Pair<Character, Double>>>> unmatchedMotifPeaks = 
      this.findUnmatchedMotifPeaks(gpsPeakLocsOne, gpsPeakLocsTwo, scorer, motifThreshold);
    this.writeUnmatchedPeaks(outFile, tfOneName, tfTwoName, unmatchedMotifPeaks, peaksByPointOne, peaksByPointTwo);

    System.out.println("DONE!");
  }


  /**
   * Find pairs of peaks that are closest to each other
   * @param gpsPeakLocsOne
   * @param gpsPeakLocsTwo
   * @return
   */
  public Pair<TreeMap<Point, Pair<Point, Point>>, TreeMap<Integer, List<Pair<Point, Point>>>> findMatchedPeaks(ArrayList<Point> gpsPeakLocsOne, TreeSet<Point> gpsPeakLocsTwo) {
    TreeMap<Point, Pair<Point, Point>> pairsByPoint = new TreeMap<Point, Pair<Point, Point>>();
    TreeMap<Integer, List<Pair<Point, Point>>> pairsByDist = new TreeMap<Integer, List<Pair<Point, Point>>>();


    ArrayList<Point> currList = gpsPeakLocsOne;
    ArrayList<Point> nextList = new ArrayList<Point>();

    while (!currList.isEmpty()) {
      for (Point peakLocOne : currList) {
        //find the points from gpsPeakLocsTwo closest to the current point from
        //gpsPeakLocsOne (and not closer to any other point from gpsPeakLocsOne)
        Point floor = findFloor(peakLocOne, gpsPeakLocsTwo, pairsByPoint, pairsByDist, nextList);
        Point ceil = findCeil(peakLocOne, gpsPeakLocsTwo, pairsByPoint, pairsByDist, nextList);
        int floorDist = Integer.MAX_VALUE;
        int ceilDist = Integer.MAX_VALUE;
        if (floor != null) {
          floorDist = (int)Math.abs(peakLocOne.offset(floor));
        }
        if (ceil != null) {
          ceilDist = (int)Math.abs(peakLocOne.offset(ceil));
        }
        //pick the closer of the two points and add the mapping
        if ((floorDist < MAX_DIST) || (ceilDist < MAX_DIST)) {
          if (floorDist < ceilDist) {
            Pair<Point, Point> pointPair = new Pair<Point, Point>(peakLocOne, floor);
            //point2PairMap.put(peakLocOne, pointPair);
            pairsByPoint.put(floor, pointPair);
            if (!pairsByDist.containsKey(floorDist)) {
              pairsByDist.put(floorDist, new ArrayList<Pair<Point, Point>>());
            }
            pairsByDist.get(floorDist).add(pointPair);
          }
          else {
            Pair<Point, Point> pointPair = new Pair<Point, Point>(peakLocOne, ceil);
            //point2PairMap.put(peakLocOne, pointPair);
            pairsByPoint.put(ceil, pointPair);
            if (!pairsByDist.containsKey(ceilDist)) {
              pairsByDist.put(ceilDist, new ArrayList<Pair<Point, Point>>());
            }
            pairsByDist.get(ceilDist).add(pointPair);       
          }
        }
      }
      currList = nextList;
      nextList = new ArrayList<Point>();
    }

    return new Pair<TreeMap<Point, Pair<Point, Point>>, TreeMap<Integer, List<Pair<Point, Point>>>>(pairsByPoint, pairsByDist);
  }


  /**
   * Find peaks which occur near a motif but not near a peak from the other factor 
   */
  public Pair<TreeMap<Point, Pair<Integer, Pair<Character, Double>>>, TreeMap<Point, Pair<Integer, Pair<Character, Double>>>> findUnmatchedMotifPeaks(ArrayList<Point> gpsPeakLocsOne, TreeSet<Point> gpsPeakLocsTwo, 
      WeightMatrixScorer scorer, double motifThreshold) {

    TreeMap<Point, Pair<Point, Point>> pairsByPoint = new TreeMap<Point, Pair<Point, Point>>();
    TreeMap<Integer, List<Pair<Point, Point>>> pairsByDist = new TreeMap<Integer, List<Pair<Point, Point>>>();
    ArrayList<Point> nextList = new ArrayList<Point>();

    ArrayList<Point> peakLocsTwoToCheck = new ArrayList<Point>(gpsPeakLocsTwo);

    ArrayList<Point> lonePeaksOne = new ArrayList<Point>();
    ArrayList<Point> lonePeaksTwo = new ArrayList<Point>();

    TreeSet<Point> nearPeaksOne = new TreeSet<Point>();

    for (Point peakLocOne : gpsPeakLocsOne) {
      Point floor = findFloor(peakLocOne, gpsPeakLocsTwo, pairsByPoint, pairsByDist, nextList);
      Point ceil = findCeil(peakLocOne, gpsPeakLocsTwo, pairsByPoint, pairsByDist, nextList);

      int floorDist = Integer.MAX_VALUE;
      int ceilDist = Integer.MAX_VALUE;
      if (floor != null) {
        floorDist = (int)Math.abs(peakLocOne.offset(floor));
      }
      if (ceil != null) {
        ceilDist = (int)Math.abs(peakLocOne.offset(ceil));
      }
      //pick the closer of the two points and add the mapping
      if ((floorDist >= MIN_DIST) && (ceilDist >= MIN_DIST)) {
        lonePeaksOne.add(peakLocOne);
      }
      else {
        nearPeaksOne.add(peakLocOne);
        if (floorDist < MIN_DIST) {
          peakLocsTwoToCheck.remove(floor);
        }
        if (ceilDist < MIN_DIST) {
          peakLocsTwoToCheck.remove(ceil);
        }
      }
    }

    for (Point peakLocTwo : peakLocsTwoToCheck) {
      Point floor = findFloor(peakLocTwo, nearPeaksOne, pairsByPoint, pairsByDist, nextList);
      Point ceil = findCeil(peakLocTwo, nearPeaksOne, pairsByPoint, pairsByDist, nextList);

      int floorDist = Integer.MAX_VALUE;
      int ceilDist = Integer.MAX_VALUE;
      if (floor != null) {
        floorDist = (int)Math.abs(peakLocTwo.offset(floor));
      }
      if (ceil != null) {
        ceilDist = (int)Math.abs(peakLocTwo.offset(ceil));
      }
      //pick the closer of the two points and add the mapping
      if ((floorDist >= MIN_DIST) && (ceilDist >= MIN_DIST)) {
        lonePeaksTwo.add(peakLocTwo);
      }
    }

    TreeMap<Point, Pair<Integer, Pair<Character, Double>>> motifPeaksOne = new TreeMap<Point, Pair<Integer, Pair<Character, Double>>>();
    TreeMap<Point, Pair<Integer, Pair<Character, Double>>> motifPeaksTwo = new TreeMap<Point, Pair<Integer, Pair<Character, Double>>>();

    int i = 0;
    System.out.println("Checking for motifs, Points to check: " + (lonePeaksOne.size() + lonePeaksTwo.size()));
    for (Point peakLocOne : lonePeaksOne) {      
      Pair<Integer, Pair<Character, Double>> motifData = this.findNearestMotif(peakLocOne, scorer, motifThreshold);
      if (motifData.car() < MOTIF_DISTANCE) {
        motifPeaksOne.put(peakLocOne, motifData);
      }
      if (i % 1000 == 0) {
        if (i % 10000 == 0) {
          System.out.println(i);
        }
        else {
          System.out.print(i);
        }
      }
      else if (i % 100==0) {
        System.out.print(".");
      }
      i++;
    }

    for (Point peakLocTwo : lonePeaksTwo) {
      Pair<Integer, Pair<Character, Double>> motifData = this.findNearestMotif(peakLocTwo, scorer, motifThreshold);
      if (motifData.car() < MOTIF_DISTANCE) {
        motifPeaksTwo.put(peakLocTwo, motifData);
      }
      if (i % 1000 == 0) {
        if (i % 10000 == 0) {
          System.out.println(i);
        }
        else {
          System.out.print(i);
        }
      }
      else if (i % 100==0) {
        System.out.print(".");
      }
      i++;
    }

    return new Pair<TreeMap<Point, Pair<Integer, Pair<Character, Double>>>, TreeMap<Point, Pair<Integer, Pair<Character, Double>>>>(motifPeaksOne, motifPeaksTwo);
  }

  /**
   * Find the point from the testSet that is closest to pOne and closer to 
   * pOne than any other point that it may already be paired with. 
   * @param pOne
   * @param nextList
   * @return
   */
  public Point findFloor(Point pOne, TreeSet<Point> testSet, TreeMap<Point, Pair<Point, Point>> pairByPointMap, 
      TreeMap<Integer, List<Pair<Point, Point>>> pairsByDistMap, List<Point> nextList) {

    Point floor = null;

    //figure out the maximum floor point given the specified Max Dist
    int maxFloorLoc = Math.max(0, pOne.getLocation()-MAX_DIST);
    Point maxFloor = new Point(pOne.getGenome(), pOne.getChrom(), maxFloorLoc);
    //get the set of all points with locations "before" pOne and within range MAX DIST
    NavigableSet<Point> floorSet = testSet.headSet(pOne, true).tailSet(maxFloor, true);


    //find the closest point to pOne that is closer to pOne than to any other
    //point that it may already be paired with. 
    Iterator<Point> iter = floorSet.descendingIterator();
    while (iter.hasNext()) {
      floor = iter.next();
      if (!floor.getChrom().equals(pOne.getChrom())) {
        //there are no more points to check on the same chromosome as pOne
        return null;
      }
      if (!pairByPointMap.containsKey(floor)) {
        return floor;
      }
      else {
        //figure out if pOne is closer to floor than the other point it's
        //paired with
        int dist = (int)Math.abs(pOne.offset(floor));
        Pair<Point, Point> otherPair = pairByPointMap.get(floor);
        int otherDist = (int)Math.abs(otherPair.car().offset(otherPair.cdr()));
        if (otherDist > dist) {
          Point bumpedPoint;
          if (otherPair.car() == floor) {
            bumpedPoint = otherPair.car();
          }
          else {
            bumpedPoint = otherPair.cdr();
          }
          nextList.add(bumpedPoint);
          pairsByDistMap.get(otherDist).remove(otherPair);
          pairByPointMap.remove(floor);
          return floor;
        }
      }  		
    }
    //if no more points than return null
    return null;
  }

  public Point findCeil(Point pOne, TreeSet<Point> testSet, TreeMap<Point, Pair<Point, Point>> pairByPointMap, 
      TreeMap<Integer, List<Pair<Point, Point>>> pairsByDistMap, List<Point> nextList) {
    Point ceil = null;

    //figure out the maximum floor point given the specified Max Dist
    int maxCeilLoc = Math.min(pOne.getGenome().getChromLength(pOne.getChrom()), pOne.getLocation()+MAX_DIST);
    Point maxCeil = new Point(pOne.getGenome(), pOne.getChrom(), maxCeilLoc);
    //get the set of all points with locations "after" pOne and within range MAX DIST     
    NavigableSet<Point> ceilSet = testSet.tailSet(pOne, true).headSet(maxCeil, true);

    Iterator<Point> iter = ceilSet.iterator();
    while (iter.hasNext()) {
      ceil = iter.next();
      if (!ceil.getChrom().equals(pOne.getChrom())) {
        return null;
      }

      if (!pairByPointMap.containsKey(ceil)) {
        return ceil;
      }
      else {
        //figure out if pOne is closer to floor than the other point it's
        //paired with
        int dist = (int)Math.abs(pOne.offset(ceil));
        Pair<Point, Point> otherPair = pairByPointMap.get(ceil);
        int otherDist = (int)Math.abs(otherPair.car().offset(otherPair.cdr()));
        if (otherDist > dist) {
          Point bumpedPoint;
          if (otherPair.car() == ceil) {
            bumpedPoint = otherPair.car();
          }
          else {
            bumpedPoint = otherPair.cdr();
          }
          nextList.add(bumpedPoint);
          pairsByDistMap.get(otherDist).remove(otherPair);
          pairByPointMap.remove(ceil);
          return ceil;
        }
      }  		
    }
    //if no more points than return null
    return null;  	
  }



  public Pair<Integer, Pair<Character, Double>> findNearestMotif(Point peak, WeightMatrixScorer scorer, double motifThreshold) {
    int motif_5prime_offset = Integer.MAX_VALUE;
    double motif_score = Double.NEGATIVE_INFINITY;
    char motif_strand = 'x';

    Region r= peak.expand(MOTIF_DISTANCE);
    WeightMatrixScoreProfile profiler = scorer.execute(r);
    int halfWidth = profiler.getMatrix().length()/2;
    //search from BS outwards
    for(int z=0; z<=MOTIF_DISTANCE; z++){
      double leftScore = Double.NEGATIVE_INFINITY;
      double rightScore = Double.NEGATIVE_INFINITY;
      if ((MOTIF_DISTANCE+z) < profiler.length()) {
        rightScore= profiler.getMaxScore(MOTIF_DISTANCE+z);        
        if(rightScore>=motifThreshold){
          motif_score = rightScore;
          motif_strand = profiler.getMaxStrand(MOTIF_DISTANCE+z);
          if (motif_strand == '+') {
            motif_5prime_offset = z;
          }
          else {
            motif_5prime_offset = z + (profiler.getMatrix().length()-1);
          }					
        }
      }
      if ((MOTIF_DISTANCE-z) >= 0) {
        leftScore= profiler.getMaxScore(MOTIF_DISTANCE-z);
        if ((leftScore>=motifThreshold) && (leftScore > motif_score)) {					
          motif_score = leftScore;
          motif_strand = profiler.getMaxStrand(MOTIF_DISTANCE-z);
          if (motif_strand == '+') {
            motif_5prime_offset = -z;
          }
          else {
            motif_5prime_offset = -z + (profiler.getMatrix().length()-1);
          }					
        }
      }
      //if either leftScore or rightScore was above the thresh then stop search
      if (motif_score > motifThreshold) {
        break;
      }
      // if motif score at this position is too small, and reach end of region
      if (z==r.getWidth()/2){
        motif_5prime_offset = Integer.MAX_VALUE;
        motif_score = Math.max(leftScore, rightScore);
      }
    }
    return new Pair<Integer, Pair<Character, Double>>(motif_5prime_offset, new Pair<Character, Double>(motif_strand, motif_score));
  }


  public void writeMatchedPeaks(String outFile, String tfOneName, String tfTwoName, boolean printZeroDist, boolean requireMotifContainment,
      WeightMatrixScorer scorer, WeightMatrix motif, double motifThreshold, TreeMap<Point, Pair<Point, Point>> pairsByPoint,
      TreeMap<Integer, List<Pair<Point, Point>>> pairsByDist, TreeMap<Point, GPSPeak> peaksByPointOne, TreeMap<Point, GPSPeak> peaksByPointTwo) {
    System.out.println("Checking for motifs and Writing file, Point Pairs to check: " + pairsByPoint.size());
    LineByLineFileWriter lblfw = null;
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(outFile);
      lblfw.writeLine("TF1 = " + tfOneName + "\nTF2 = " + tfTwoName);
      lblfw.writeLine("TF1 Peak" + "\t" + "TF2 Peak" + "\t" + "Interpeak dist" + "\t" + "stranded dist" + "\t" + "motif 5' start" 
          + "\t" + "motif strand" + "\t" + "motif score" + "\t" + "TF1 pos on motif" + "\t" + "TF2 pos on motif" + "\t" 
          + "TF1 peak size" + "\t" + "TF2 peak size");

      int i = 0;
      int startDist = 0;
      if (!printZeroDist) {
        startDist = 1;
      }
      for (int k = startDist; k < MAX_DIST; k++) {
        List<Pair<Point, Point>>pairList = pairsByDist.get(k);
        if (pairList == null) {
          continue;
        }
        //for (List<Pair<Point, Point>> pairList : pointPairs.values()) {        
        for (Pair <Point, Point> pointPair : pairList) {
          if (i % 1000 == 0) {
            if (i % 10000 == 0) {
              System.out.println(i);
            }
            else {
              System.out.print(i);
            }
          }
          else if (i % 100==0) {
            System.out.print(".");
          }

          Point pOne = pointPair.car();
          Point pTwo = pointPair.cdr();
          Pair<Integer, Pair<Character, Double>> motifData = this.findNearestMotif(pOne, scorer, motifThreshold);

          int dist = pOne.offset(pTwo);
          int strandedDist = Integer.MAX_VALUE;
          int motif5primeStartOffset = Integer.MAX_VALUE;
          int fOnePosOnMotif = Integer.MAX_VALUE;
          int fTwoPosOnMotif = Integer.MAX_VALUE;  

          if (motifData.car() <= MOTIF_DISTANCE) {           
            strandedDist = dist;
            motif5primeStartOffset = motifData.car();
            fOnePosOnMotif = -motif5primeStartOffset;
            fTwoPosOnMotif = fOnePosOnMotif - dist;  
            if (motifData.cdr().car() == '-') {
              strandedDist = -dist;
              fOnePosOnMotif = motif5primeStartOffset;
              fTwoPosOnMotif = fOnePosOnMotif + dist;
            }
          }
          else {
            //try the other factor
            motifData = this.findNearestMotif(pTwo, scorer, motifThreshold);
            if (motifData.car() <= MOTIF_DISTANCE) {           
              strandedDist = dist;
              motif5primeStartOffset = motifData.car();
              fTwoPosOnMotif = -motif5primeStartOffset;
              fOnePosOnMotif = fTwoPosOnMotif + dist;               
              if (motifData.cdr().car() == '-') {
                strandedDist = -dist;
                fTwoPosOnMotif = motif5primeStartOffset;
                fOnePosOnMotif = fTwoPosOnMotif - dist;                 
              }
            }
          }

          if (motifData.car() <= MOTIF_DISTANCE) {
            if (!requireMotifContainment || ((fOnePosOnMotif >= 0) && (fTwoPosOnMotif >= 0) 
                && (fOnePosOnMotif < motif.length()) && (fTwoPosOnMotif < motif.length()))) {
              lblfw.writeLine(pOne.toString() + "\t" + pTwo.toString() + "\t" + dist + "\t" + strandedDist + "\t" + motif5primeStartOffset 
                  + "\t" + motifData.cdr().car() + "\t" + motifData.cdr().cdr() + "\t" + fOnePosOnMotif + "\t" + fTwoPosOnMotif
                  + "\t" + peaksByPointOne.get(pOne).getStrength() + "\t" + peaksByPointTwo.get(pTwo).getStrength());
            }
          }
          i++;
        }
      }
    }
    catch (IOException ioex) {
      ioex.printStackTrace();
    }
    finally {
      if (lblfw != null) {
        try {
          lblfw.closeFile();
        }
        catch (IOException ioex) {
          ioex.printStackTrace();
        }
      }     
    }
  }


  public void writeUnmatchedPeaks(String outFile, String tfOneName, String tfTwoName, 
      Pair<TreeMap<Point, Pair<Integer, Pair<Character, Double>>>, TreeMap<Point, Pair<Integer, Pair<Character, Double>>>> unmatchedMotifPeaks,
      TreeMap<Point, GPSPeak> peaksByPointOne, TreeMap<Point, GPSPeak> peaksByPointTwo) {
    System.out.println("Writing file");
    LineByLineFileWriter lblfw = null;
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(outFile + "." + tfOneName);
      lblfw.writeLine("TF1 = " + tfOneName + "\nTF2 = " + tfTwoName);
      lblfw.writeLine("TF1 Peak" + "\t" + "motif 5' start" + "\t" + "motif strand" + "\t" + "motif score" + "\t" 
          + "TF1 pos on motif" + "\t" + "TF1 peak size");

      for (Point peakLocOne : unmatchedMotifPeaks.car().keySet()) {
        Pair<Integer, Pair<Character, Double>> motifData = unmatchedMotifPeaks.car().get(peakLocOne);
        int motif5primeStartOffset = motifData.car();
        int fOnePosOnMotif = -motif5primeStartOffset; 
        if (motifData.cdr().car() == '-') {            
          fOnePosOnMotif = motif5primeStartOffset;
        }
        lblfw.writeLine(peakLocOne.toString() + "\t" + motif5primeStartOffset + "\t" + motifData.cdr().car() + "\t" + motifData.cdr().cdr() + 
            "\t" + fOnePosOnMotif + "\t" + peaksByPointOne.get(peakLocOne).getStrength());
      }	      
    }
    catch (IOException ioex) {
      ioex.printStackTrace();
    }
    finally {
      if (lblfw != null) {
        try {
          lblfw.closeFile();
        }
        catch (IOException ioex) {
          ioex.printStackTrace();
        }
      }     
    }

    //write peak two
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(outFile + "." + tfTwoName);
      lblfw.writeLine("TF1 = " + tfOneName + "\nTF2 = " + tfTwoName);
      lblfw.writeLine("TF2 Peak" + "\t" + "motif 5' start" + "\t" + "motif strand" + "\t" + "motif score" + "\t" 
          + "TF2 pos on motif" + "\t" + "TF2 peak size");

      for (Point peakLocTwo : unmatchedMotifPeaks.cdr().keySet()) {
        Pair<Integer, Pair<Character, Double>> motifData = unmatchedMotifPeaks.cdr().get(peakLocTwo);
        int motif5primeStartOffset = motifData.car();
        int fTwoPosOnMotif = -motif5primeStartOffset; 
        if (motifData.cdr().car() == '-') {            
          fTwoPosOnMotif = motif5primeStartOffset;
        }
        lblfw.writeLine(peakLocTwo.toString() + "\t" + motif5primeStartOffset + "\t" + motifData.cdr().car() + "\t" + motifData.cdr().cdr() + 
            "\t" + fTwoPosOnMotif + "\t" + peaksByPointTwo.get(peakLocTwo).getStrength());
      }       
    }
    catch (IOException ioex) {
      ioex.printStackTrace();
    }
    finally {
      if (lblfw != null) {
        try {
          lblfw.closeFile();
        }
        catch (IOException ioex) {
          ioex.printStackTrace();
        }
      }     
    }
  }


  /**
   * @param args
   */
  public static void main(String[] args) {
    try {
      String tfOneName = "Oct4";
      String tfTwoName = "Sox2";


      String GPSfilenameOne = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/YL_Oct4_ES_1_signal.peaks.txt";
      String GPSfilenameTwo = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/YL_Sox2_ES_2_signal.peaks.txt";
      String outFile = "yl_2_oct4_sox2_unmatched_motif_offsets_01.txt";
      //    String GPSfilenameOne = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/Oct4_YL_0_signal.peaks.txt";
      //    String GPSfilenameTwo = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/Sox2_rep1_YL_0_signal.peaks.txt";
      //	    String GPSfilenameOne = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/Sing_Oct4_ES_2_signal.peaks.txt";
      //	    String GPSfilenameTwo = "/afs/csail.mit.edu/group/psrg/projects/GPS/GPS_results/Sing_Sox2_ES_2_signal.peaks.txt";
      //      String outFile = "sing_oct4_sox2_motif_offsets_lax_01.txt";

      Genome genome = Organism.findGenome("mm8");
      Organism org = Organism.getOrganism("Mus musculus");     
      String motifString = "OSNT";
      String motifVersion = "YoungLab";

      double motifThreshold = 12.72;

      boolean printZeroDist = false;
      boolean requireInMotif = true;

      new BindingOffsetAnalysis(GPSfilenameOne, GPSfilenameTwo, outFile, tfOneName, tfTwoName, genome, org, 
          motifString, motifVersion, motifThreshold, printZeroDist, requireInMotif);
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
  }

}
