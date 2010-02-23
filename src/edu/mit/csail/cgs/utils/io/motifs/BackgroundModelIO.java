/**
 * 
 */
package edu.mit.csail.cgs.utils.io.motifs;

import java.io.IOException;
import java.text.ParseException;
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

  public static final String BG_LINE_COUNTS_REG_EX = "^([\\d]+)\\s*([ACGTacgt]+)\\s*(\\d[\\.?]([\\d+])\\s*";
  public static final String BG_LINE_PROBS_REG_EX = "^([\\d]+)\\s*([ACGTacgt]+)\\s*([01]\\.([\\d+])\\s*";

  public static final Pattern BG_LINE_COUNTS_PATTERN = Pattern.compile(BG_LINE_COUNTS_REG_EX);
  public static final Pattern BG_LINE_PROBS_PATTERN = Pattern.compile(BG_LINE_PROBS_REG_EX);

  private static void initBackgroundModel(BackgroundModel model, String[] lines) throws ParseException {
    // initialize background model object
    Matcher bgLineMatcher = null;
    if (model instanceof CountsBackgroundModel) {
      bgLineMatcher = BG_LINE_COUNTS_PATTERN.matcher("");
    }
    else {
      bgLineMatcher = BG_LINE_PROBS_PATTERN.matcher(""); 
    }
    int lineIndex = 0;
    for (int i = 0; i < model.getMaxKmerLen(); i++) {
      for (int j = 0; j < Math.pow(4, i + 1); j++) {
        double total = 0.0;
        bgLineMatcher.reset(lines[lineIndex]);
        if (bgLineMatcher.matches()) {
          int intVal = Integer.valueOf(bgLineMatcher.group(1));
          String mer = bgLineMatcher.group(2).toUpperCase();
          if ((intVal != j) || !mer.equals(BackgroundModel.intToSeq(j, (i + 1)))) {
            throw new ParseException("Expected index " + j + " and kmer " + BackgroundModel.intToSeq(j, (i + 1))
                + ", but got " + lines[lineIndex], 0);
          }
          double prob = Double.valueOf(bgLineMatcher.group(3));
          total = total + prob;
          model.setModelVal(mer, prob);
          lineIndex++;
          model.setModelVal(mer, prob);
        }
        else {
          throw new ParseException("Incorrectly formatted line: " + lines[lineIndex], 0);
        }
      }
    }
  }


  /**
   * 
   * @param lines
   * @return
   * @throws ParseException
   */
  private static int checkModelOrder(String[] lines) throws ParseException {
    // determine order of model and check that the number of lines is valid
    int order = 0;
    int numLines = 4;
    while (numLines < lines.length) {
      order++;
      numLines = numLines + (int) Math.pow(4, order + 1);
    }
    if (numLines != lines.length) {
      throw new ParseException("Parsed " + lines.length
          + " non-comment lines, but number of lines must be 4 + 4^2 + 4^3 + ... + 4^(k+1)", 0);
    }
    return order;
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

      int order = BackgroundModelIO.checkModelOrder(lines);

      // construct the Counts Background Model object of that order
      CountsBackgroundModel bgModel = new CountsBackgroundModel(order);

      // initialize the model
      BackgroundModelIO.initBackgroundModel(bgModel, lines);

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

      int order = BackgroundModelIO.checkModelOrder(lines);

      // construct the Markov Background Model object of that order
      MarkovBackgroundModel bgModel = new MarkovBackgroundModel(order);

      // initialize the model
      BackgroundModelIO.initBackgroundModel(bgModel, lines);
      
      String[] normKmers = MarkovBackgroundModel.isModelNormalized(bgModel);
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

      int order = BackgroundModelIO.checkModelOrder(lines);

      // construct the Frequency Background Model object of that order
      FrequencyBackgroundModel bgModel = new FrequencyBackgroundModel(order);

      // initialize the model
      BackgroundModelIO.initBackgroundModel(bgModel, lines);

      int normKmers = FrequencyBackgroundModel.isModelNormalized(bgModel);
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
  

  //Print the background model to a file
  public static void printToFile(BackgroundModel bgModel, String filename) throws IOException {   
    LineByLineFileWriter lblfw = null;
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(filename);      
      
      for (int i = 1; i <= bgModel.getMaxKmerLen(); i++) {
        for (int j = 0; j < Math.pow(4, i); j++) {
          String currKmer = BackgroundModel.intToSeq(j, i);
          lblfw.writeLine(j + "\t" + currKmer + "\t" + bgModel.getModelVal(i, j));  
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
