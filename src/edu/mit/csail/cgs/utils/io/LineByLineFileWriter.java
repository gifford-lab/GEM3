/**
 * 
 */
package edu.mit.csail.cgs.utils.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author rca
 *
 */
public class LineByLineFileWriter {
  private BufferedWriter writer = null;
  
  public LineByLineFileWriter() {
    //do nothing
  }
  
  /**
   * Open the specified file in a buffered Writeer
   * @param filename the name of the file to open
   * @throws IOException
   */
  public void openFile(String filename) throws IOException {
    FileWriter fileWriter = new FileWriter(filename);
    writer = new BufferedWriter(fileWriter);
  }
  
  
  /**
   * Close the buffered Writeer (if it's open)
   * @throws IOException
   */
  public void closeFile() throws IOException {
    if (writer != null) {
      writer.close();
      writer = null;
    }
  }
  
  
  /**
   * Write the next line from the file
   * @return the next line from the file
   * @throws IOException
   */
  public void writeLine(String line) throws IOException {    
    writer.write(line);
    writer.newLine();
  }
  
  
  /**
   * Write an entire file at once
   * @param filename the name of the file to write
   * @param lines the lines of the file
   * @throws IOException
   */
  public static void writeFile(String filename, String[] lines) throws IOException {
    LineByLineFileWriter lblfw = null;
    try {
      lblfw = new LineByLineFileWriter();
      lblfw.openFile(filename);
      for (String line : lines) {
        lblfw.writeLine(line);
      }
    }
    finally {
      if (lblfw != null) {
        lblfw.closeFile();
      }
    }    
  }  
}
