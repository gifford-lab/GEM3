/**
 * 
 */
package edu.mit.csail.cgs.utils.io.motifs;

import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mit.csail.cgs.datasets.motifs.BackgroundModel;
import edu.mit.csail.cgs.datasets.motifs.CountsBackgroundModel;
import edu.mit.csail.cgs.datasets.motifs.MarkovBackgroundModel;
import edu.mit.csail.cgs.datasets.motifs.FrequencyBackgroundModel;
import edu.mit.csail.cgs.utils.io.LineByLineFileReader;
import edu.mit.csail.cgs.utils.io.LineByLineFileWriter;

/**
 * @author rca
 * 
 */
public class BackgroundModelIO {

  public static final String BG_LINE_COUNTS_REG_EX = "^([\\d]+)\\s*([ACGTacgt]+)\\s*(\\d+)\\s*";
  public static final String BG_LINE_PROBS_REG_EX = "^([\\d]+)\\s*([ACGTacgt]+)\\s*([01][\\.\\d+]?)\\s*";

  public static final Pattern BG_LINE_COUNTS_PATTERN = Pattern.compile(BG_LINE_COUNTS_REG_EX);
  public static final Pattern BG_LINE_PROBS_PATTERN = Pattern.compile(BG_LINE_PROBS_REG_EX);

  public static void main(String[] args) throws ParseException, IOException {
    BackgroundModelIO.parseMarkovBackgroundModel("/home/rca/projects/cgs_repos/cgs/mm8.back");
    
//    String testLine = "12      TA      0.21842155511047961";
//    Matcher bgLineMatcher = BG_LINE_PROBS_PATTERN.matcher(testLine);
//    boolean found = bgLineMatcher.matches();
//    System.out.println(found);
//    if (found) {
//      System.out.println("start: " + bgLineMatcher.start() + "  end: " + bgLineMatcher.end() + "  group: " + bgLineMatcher.group());
//    }    
  }
  
  private static void initFreqBackgroundModel(FrequencyBackgroundModel model, String[] lines) throws ParseException {
    // initialize background model object

    Matcher bgLineMatcher = BG_LINE_PROBS_PATTERN.matcher(""); 
    int lineIndex = 0;
    for (int i = 0; i < model.getMaxKmerLen(); i++) {
    	HashMap<String, Double> freqs = new HashMap<String, Double>();
      for (int j = 0; j < Math.pow(4, i + 1); j++) {
        bgLineMatcher.reset(lines[lineIndex]);
        if (bgLineMatcher.matches()) {
          int intVal = Integer.valueOf(bgLineMatcher.group(1));
          String mer = bgLineMatcher.group(2).toUpperCase();
          if ((intVal != j) || !mer.equals(BackgroundModel.int2seq(j, (i + 1)))) {
            throw new ParseException("Expected index " + j + " and kmer " + BackgroundModel.int2seq(j, (i + 1))
                + ", but got " + lines[lineIndex], 0);
          }
          double prob = Double.valueOf(bgLineMatcher.group(3));
          freqs.put(mer, prob);
          lineIndex++;
        }
        else {
          throw new ParseException("Incorrectly formatted line: " + lines[lineIndex], 0);
        }
      }
      model.setKmerFrequencies(freqs);
    }
  }

  
  private static void initMarkovBackgroundModel(MarkovBackgroundModel model, String[] lines) throws ParseException {
    // initialize background model object

    Matcher bgLineMatcher = BG_LINE_PROBS_PATTERN.matcher(""); 
    int lineIndex = 0;
    for (int i = 0; i < model.getMaxKmerLen(); i++) {
      for (int j = 0; j < Math.pow(4, i + 1); j += 4) {
      	double[] probs = new double[4];
      	for (int k = 0; k < 4; k++) {
          bgLineMatcher.reset(lines[lineIndex]);
          if (bgLineMatcher.matches()) {
            int intVal = Integer.valueOf(bgLineMatcher.group(1));
            String mer = bgLineMatcher.group(2).toUpperCase();
            if ((intVal != (j+k)) || !mer.equals(BackgroundModel.int2seq(j, (i + 1)))) {
              throw new ParseException("Expected index " + j + " and kmer " + BackgroundModel.int2seq(j, (i + 1))
                  + ", but got " + lines[lineIndex], 0);
            }
            probs[k] = Double.valueOf(bgLineMatcher.group(3));
            lineIndex++;
          }
          else {
            throw new ParseException("Incorrectly formatted line: " + lines[lineIndex], 0);
          }      		
      	}
      	model.setMarkovProb(BackgroundModel.int2seq(j, i).substring(0, i-1), probs[0], probs[1], probs[2], probs[3]);
      }
    }
  }

  private static void initCountsBackgroundModel(CountsBackgroundModel model, String[] lines) throws ParseException {
    // initialize background model object
    Matcher bgLineMatcher = BG_LINE_COUNTS_PATTERN.matcher("");
    int lineIndex = 0;
    for (int i = 0; i < model.getMaxKmerLen(); i++) {
      for (int j = 0; j < Math.pow(4, i + 1); j++) {
        bgLineMatcher.reset(lines[lineIndex]);
        if (bgLineMatcher.matches()) {
          int intVal = Integer.valueOf(bgLineMatcher.group(1));
          String mer = bgLineMatcher.group(2).toUpperCase();
          if ((intVal != j) || !mer.equals(BackgroundModel.int2seq(j, (i + 1)))) {
            throw new ParseException("Expected index " + j + " and kmer " + BackgroundModel.int2seq(j, (i + 1))
                + ", but got " + lines[lineIndex], 0);
          }
          int count = Integer.valueOf(bgLineMatcher.group(3));
          model.setKmerCount(mer, count);
          lineIndex++;
        }
        else {
          throw new ParseException("Incorrectly formatted line: " + lines[lineIndex], 0);
        }
      }
    }
  }
  
  
  /**
   * Checks the number of lines to make sure there are a proper number for all 
   * kmers up to the determined maximum length
   * @param lines
   * @return
   * @throws ParseException
   */
  private static int checkModelMaxKmerLen(String[] lines) throws ParseException {
    // determine the max kmer length for model and check that the number of lines is valid
    int maxKmerLen = 1;
    int numLines = 4;
    while (numLines < lines.length) {
      maxKmerLen++;
      numLines = numLines + (int) Math.pow(4, maxKmerLen);
    }
    if (numLines != lines.length) {
      throw new ParseException("Parsed " + lines.length
          + " non-comment lines, but number of lines must be 4 + 4^2 + 4^3 + ... + 4^(k+1)", 0);
    }
    return maxKmerLen;
  }

  
  /**
   * 
   * @param filename
   * @return
   * @throws IOException
   * @throws ParseException
   */
  public static CountsBackgroundModel parseCountsBackgroundModel(String filename) throws IOException, ParseException {
    LineByLineFileReader lblfr = null;

    try {
      String[] lines = LineByLineFileReader.readFile(filename, LineByLineFileReader.DEFAULT_COMMENT_PREFIXES, true);

      int maxKmerLen = BackgroundModelIO.checkModelMaxKmerLen(lines);

      // construct the Counts Background Model object of that order
      CountsBackgroundModel bgModel = new CountsBackgroundModel(filename, null, maxKmerLen);

      // initialize the model
      BackgroundModelIO.initCountsBackgroundModel(bgModel, lines);

      return bgModel;
    }
    finally {
      if (lblfr != null) {
        lblfr.closeFile();
      }
    }
  }


  /**
   * 
   * @param filename
   * @return
   * @throws IOException
   * @throws ParseException
   */
  public static MarkovBackgroundModel parseMarkovBackgroundModel(String filename) throws IOException, ParseException {
    LineByLineFileReader lblfr = null;

    try {
      String[] lines = LineByLineFileReader.readFile(filename, LineByLineFileReader.DEFAULT_COMMENT_PREFIXES, true);

      int maxKmerLen = BackgroundModelIO.checkModelMaxKmerLen(lines);

      /**
       * construct the Markov Background Model object for that kmer length
       * Note: the markov order will be maxKmerLen - 1
       */      
      MarkovBackgroundModel bgModel = new MarkovBackgroundModel(filename, null, maxKmerLen);

      // initialize the model
      BackgroundModelIO.initMarkovBackgroundModel(bgModel, lines);
      
      String[] normKmers = bgModel.verifyNormalization();
      if (normKmers != null) {
        throw new ParseException("Model is not normalized for kmers: " + normKmers[0] + "," + normKmers[1] + ","
            + normKmers[2] + "," + normKmers[3], -1);
      }

      return bgModel;
    }
    finally {
      if (lblfr != null) {
        lblfr.closeFile();
      }
    }
  }


  /**
   * 
   * @param filename
   * @return
   * @throws IOException
   * @throws ParseException
   */
  public static FrequencyBackgroundModel parseFreqBackgroundModel(String filename) throws IOException, ParseException {
    LineByLineFileReader lblfr = null;

    try {
      String[] lines = LineByLineFileReader.readFile(filename, LineByLineFileReader.DEFAULT_COMMENT_PREFIXES, true);

      int maxKmerLen = BackgroundModelIO.checkModelMaxKmerLen(lines);

      // construct the Frequency Background Model object for that kmer length
      FrequencyBackgroundModel bgModel = new FrequencyBackgroundModel(filename, null, maxKmerLen);

      // initialize the model
      BackgroundModelIO.initFreqBackgroundModel(bgModel, lines);

      int normKmers = bgModel.verifyNormalization();
      if (normKmers != -1) {
        throw new ParseException("Model is not normalized for kmers of length " + normKmers, -1);
      }
      
      return bgModel;
    }
    finally {
      if (lblfr != null) {
        lblfr.closeFile();
      }
    }
  }
  

  /**
   * Write out the background model to a file
   * @param bgModel
   * @param filename
   * @throws IOException
   */
  public static void printProbsToFile(BackgroundModel bgModel, String filename) throws IOException {   
    LineByLineFileWriter lblfw = null;
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(filename);      
      
      for (int i = 1; i <= bgModel.getMaxKmerLen(); i++) {
        for (int j = 0; j < Math.pow(4, i); j++) {
          String currKmer = BackgroundModel.int2seq(j, i);
          if (bgModel instanceof MarkovBackgroundModel) {
          	lblfw.writeLine(j + "\t" + currKmer + "\t" + bgModel.getMarkovProb(i, j));
          }
          else if (bgModel instanceof FrequencyBackgroundModel) {          	
          	lblfw.writeLine(j + "\t" + currKmer + "\t" + ((FrequencyBackgroundModel)bgModel).getFrequency(i, j));
          }
          else if (bgModel instanceof CountsBackgroundModel) {
          	lblfw.writeLine(j + "\t" + currKmer + "\t" + ((CountsBackgroundModel)bgModel).getKmerCount(i, j));
          }
        }
      }
      lblfw.writeLine("");
    } 
    finally { 
      if (lblfw != null) {
        lblfw.closeFile();
      }
    }
  }
}
