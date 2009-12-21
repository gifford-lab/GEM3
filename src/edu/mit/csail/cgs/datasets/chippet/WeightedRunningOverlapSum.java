/*
 * Created on Aug 22, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.chippet;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Interval;

public class WeightedRunningOverlapSum {

  private Genome genome;

  private String chrom;

  private TreeMap<Integer, Double> changes;


  public WeightedRunningOverlapSum(Genome g, String c) {
    genome = g;
    chrom = c;
    changes = new TreeMap<Integer, Double>();
  }


  public void clear() {
    changes.clear();
  }


  public void addInterval(int start, int end) {
    this.addWeightedInterval(start, end, 1.0);
  }


  public void addWeightedInterval(int start, int end, double weight) {
    if (changes.containsKey(start)) {
      changes.put(start, changes.get(start) + weight);
    }
    else {
      changes.put(start, weight);
    }

    /**
     * I'm assuming that intervals, like regions, are inclusive of their end
     * points... This needs to put the changepoint at one past the position
     * where the interval ends to account for intervals being inclusive. -Bob
     */
    if (changes.containsKey(end + 1)) {
      changes.put(end + 1, changes.get(end + 1) - weight);
    }
    else {
      changes.put(end + 1, -weight);
    }
  }


  public void addInterval(Interval intv) {
    this.addWeightedInterval(intv, 1.0);
  }


  public void addWeightedInterval(Interval intv, double weight) {
    int start = intv.start;
    int end = intv.end;
    if (changes.containsKey(start)) {
      changes.put(start, changes.get(start) + weight);
    }
    else {
      changes.put(start, weight);
    }

    /**
     * I'm assuming that intervals, like regions, are inclusive of their end
     * points... This needs to put the changepoint at one past the position
     * where the interval ends to account for intervals being inclusive. -Bob
     */
    if (changes.containsKey(end + 1)) {
      changes.put(end + 1, changes.get(end + 1) - weight);
    }
    else {
      changes.put(end + 1, -weight);
    }

  }


  public double[][] getChangePoints() {
    double[][] array = new double[changes.size()][];
    int i = 0;
    for (int change : changes.keySet()) {
      double[] pair = new double[2];
      pair[0] = change;
      pair[1] = changes.get(change);
      array[i++] = pair;
    }
    return array;
  }


  public void addRegion(Region r) {
    this.addWeightedRegion(r, 1.0);
  }


  public void addWeightedRegion(Region r, double weight) {
    if (!genome.equals(r.getGenome())) {
      throw new IllegalArgumentException(r.getGenome().toString());
    }
    if (!chrom.equals(r.getChrom())) {
      throw new IllegalArgumentException(r.getChrom());
    }

    int start = r.getStart();
    int end = r.getEnd();
    if (changes.containsKey(start)) {
      changes.put(start, changes.get(start) + weight);
    }
    else {
      changes.put(start, weight);
    }

    /**
     * This needs to put the changepoint at one past the position where the
     * region ends to account for regions being inclusive. -Bob
     */
    if (changes.containsKey(end + 1)) {
      changes.put(end + 1, changes.get(end + 1) - weight);
    }
    else {
      changes.put(end + 1, -weight);
    }
  }


  public double getMaxOverlap() {
    double max = 0;
    double running = 0;
    for (int pc : changes.keySet()) {
      double delta = changes.get(pc);
      running += delta;
      max = Math.max(max, running);
    }
    return max;
  }


  public double getMaxOverlap(int start, int end) {
    double max = 0;
    double running = 0;
    Map<Integer, Double> map = changes.headMap(end);
    for (int pc : map.keySet()) {
      double delta = map.get(pc);
      running += delta;
      if (pc >= start && pc <= end) {
        max = Math.max(max, running);
      }
    }
    return max;
  }


  public double countOverlapping(int start, int end) {
    double count = 0;
    Map<Integer, Double> map = changes.headMap(end);
    double running = 0;
    for (int pc : map.keySet()) {
      double delta = map.get(pc);
      running += delta;
      if (pc >= start && delta == -1) {
        count += 1;
      }
    }
    count += running;
    return count;
  }


  public Collection<Region> collectRegions(double threshold) {
    if (threshold <= 0) {
      throw new IllegalArgumentException(String.valueOf(threshold));
    }

    LinkedList<Region> regions = new LinkedList<Region>();
    double runningSum = 0;

    int rstart = -1, rend = -1;

    for (int pc : changes.keySet()) {
      double delta = changes.get(pc);
      if (runningSum < threshold && runningSum + delta >= threshold) {
        rstart = pc;
      }

      if (runningSum >= threshold && runningSum + delta < threshold) {
        // subtract 1 from the endpoint because regions are inclusive
        rend = pc - 1;
        Region r = new Region(genome, chrom, rstart, rend);
        regions.addLast(r);
        rstart = rend = -1;
      }

      runningSum += delta;
      if (runningSum < 0) {
        throw new IllegalStateException(String.valueOf(runningSum));
      }
    }

    if (runningSum > 0) {
      throw new IllegalStateException(String.valueOf(runningSum));
    }
    return regions;
  }


  public TreeMap<Integer, Double> getChangeMap() {
    return changes;
  }
}
