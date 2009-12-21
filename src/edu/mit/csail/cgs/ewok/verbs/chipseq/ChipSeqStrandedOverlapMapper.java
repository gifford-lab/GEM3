package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.sql.*;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.chippet.WeightedRunningOverlapSum;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class ChipSeqStrandedOverlapMapper implements Closeable, Mapper<Region, WeightedRunningOverlapSum[]> {

  public static final int TOTAL_SUM_INDEX = 0;
  public static final int POS_SUM_INDEX = 1;
  public static final int NEG_SUM_INDEX = 2;
  
  private ChipSeqLoader loader;

  private LinkedList<ChipSeqAlignment> alignments;

  private int extension;

  private int shift = 0;

  private boolean shifting = false;


  public ChipSeqStrandedOverlapMapper(ChipSeqLocator loc, int extension) throws SQLException, IOException {
    this.extension = extension;
    loader = new ChipSeqLoader();
    alignments = new LinkedList<ChipSeqAlignment>();

    try {
      alignments.addAll(loc.loadAlignments(loader));
    }
    catch (SQLException e) {
      e.printStackTrace(System.err);
    }
    catch (NotFoundException e) {
      e.printStackTrace();
    }
  }


  public WeightedRunningOverlapSum[] execute(Region a) {
    try {
      Genome g = a.getGenome();
      String chrom = a.getChrom();
      WeightedRunningOverlapSum[] sums = { new WeightedRunningOverlapSum(g, chrom), 
          new WeightedRunningOverlapSum(g, chrom), 
          new WeightedRunningOverlapSum(g, chrom) };
      
      Collection<ChipSeqHit> hits = loader.loadByRegion(alignments, a);
      for (ChipSeqHit hit : hits) {
        if (hit.getStrand() == '+') {
          if (shifting) {
            int intervalStart = hit.getStart() + shift - (extension / 2);
            int intervalEnd = hit.getEnd() + shift + (extension / 2); 
            sums[POS_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());            
          }
          else {
            sums[POS_SUM_INDEX].addWeightedInterval(hit.getStart(), hit.getEnd() + extension, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(hit.getStart(), hit.getEnd() + extension, hit.getWeight());
          }
        }
        else {
          if (shifting) {
            int intervalStart = hit.getStart() - shift - (extension / 2);
            int intervalEnd = hit.getEnd() - shift + (extension / 2);
            sums[NEG_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
          }
          else {
            sums[NEG_SUM_INDEX].addWeightedInterval(hit.getStart() - extension, hit.getEnd(), hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(hit.getStart() - extension, hit.getEnd(), hit.getWeight());
          }
        }
      }
      return sums;
    }
    catch (Exception sqlex) {
      throw new DatabaseException(sqlex.toString(), sqlex);
    }
  }


  public void setExtension(int e) {
    extension = e;
  }


  public void setShift(int s) {
    shift = s;
    if (s > 0)
      shifting = true;
    else shifting = false;
  }


  public void close() {
    if (loader != null) {
      loader.close();
      loader = null;
      alignments.clear();
    }
  }


  public boolean isClosed() {
    return loader == null;
  }
}