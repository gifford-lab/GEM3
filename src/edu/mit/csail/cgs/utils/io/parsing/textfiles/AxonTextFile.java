package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;

public class AxonTextFile extends RowsAndColumns {

    public AxonTextFile (String filename, 
                           RowsColumnsHandler handler) throws IOException, FileNotFoundException {        
        super(filename,handler);
        System.err.println("NEW ATF " + filename);
    }
    public AxonTextFile(InputStreamReader isr,
                          RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }

    public void processFileOpen() throws IOException {
        String firstline = getReader().readLine();
        if (!firstline.matches("^ATF.*")) {
            throw new IOException ("File missing ATF marker at start");
        }
        String secondline = getReader().readLine();
        String hinfo[] = secondline.split("\\s+",-1);
        int hlines = Integer.parseInt(hinfo[0]);
        int colcount = Integer.parseInt(hinfo[1]);
        String line;
        System.err.println("Parsing " + hlines + " head lines");
        while ((hlines-- > 0) && ((line = getReader().readLine()) != null)) {
            System.err.println(": " + line);
            processHeaderLine(line);
        }
        super.processFileOpen();
    }

    

}
