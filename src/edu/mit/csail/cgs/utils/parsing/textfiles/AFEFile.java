package edu.mit.csail.cgs.utils.parsing.textfiles;

import java.io.*;

public class AFEFile extends RowsAndColumns {
    public AFEFile (String filename, 
                         RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(filename,handler);
    }
    public AFEFile(InputStreamReader isr,
                        RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }
    
    public void processFileOpen() throws IOException {
        for (int i = 0; i < 9; i++) {
            processHeaderLine(getReader().readLine());
        }
    }
}
