/**
 * 
 */
package edu.mit.csail.cgs.utils.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * @author rca
 *
 */
public class LineByLineFileReader {

  public static final String[] DEFAULT_COMMENT_PREFIXES = new String[] { "#", "//" }; 
  
  private BufferedReader reader = null;
  
  public LineByLineFileReader() {
    //do nothing
  }
  
  /**
   * Open the specified file in a buffered reader
   * @param filename the name of the file to open
   * @throws IOException
   */
  public void openFile(String filename) throws IOException {
    FileReader fileReader = new FileReader(filename);
    reader = new BufferedReader(fileReader);
  }
  
  
  /**
   * Close the buffered reader (if it's open)
   * @throws IOException
   */
  public void closeFile() throws IOException {
    if (reader != null) {
      reader.close();
      reader = null;
    }
  }
  
  
  /**
   * Read the next line from the file
   * @return the next line from the file, or null if EOF is reached
   * @throws IOException
   */
  public String readLine() throws IOException {    
    return reader.readLine();
  }
  
  
  /**
   * Read the next line that isn't a comment, as indicated by starting with one
   * of the specified prefixes
   *
   * @param commentPrefixes prefixes indicating a comment line
   * @return the next non-comment line, or null if EOF is reached
   * @throws IOException
   */
  public String readLine(String[] commentPrefixes) throws IOException {
    String currLine = reader.readLine();

    while ((currLine != null) && LineByLineFileReader.isComment(currLine, commentPrefixes)) {
      currLine = reader.readLine();
    }
    
    return currLine;
  }
  
  
  /**
   * Read an entire file in one pass and return an array of lines
   * No lines are treated as comments
   * @param filename
   * @return
   * @throws IOException
   */
  public static String[] readFile(String filename) throws IOException {
    return LineByLineFileReader.readFile(filename, new String[] {});
  }
  
  
  /**
   * Read an entire file in one pass and return an array of the lines
   * @param filename the name of the file to read
   * @param commentPrefixes prefixes indicating a comment line
   * @return
   * @throws IOException
   */
  public static String[] readFile(String filename, String[] commentPrefixes) throws IOException {
    ArrayList<String> lines = new ArrayList<String>();
    LineByLineFileReader lblfr = null;
    try {
      lblfr = new LineByLineFileReader();
      lblfr.openFile(filename);
      String currLine = lblfr.readLine(commentPrefixes);
      while (currLine != null) {
        lines.add(currLine);
        currLine = lblfr.readLine();
      }
    }
    finally {
      if (lblfr != null) {
        lblfr.closeFile();
      }
    }
    return lines.toArray(new String[] {""});
  }  
  
  
  /**
   * Check whether the line is a comment, as indicated by starting with any of
   * the specified prefixes 
   * @param line
   * @param commentPrefixes
   * @return true if a comment, false if not a comment
   */
  public static boolean isComment(String line, String[] commentPrefixes) {
    for (String commentPrefix : commentPrefixes) {
      if (line.startsWith(commentPrefix)) {
        return true;
      }
    }
    return false;
  }
}
