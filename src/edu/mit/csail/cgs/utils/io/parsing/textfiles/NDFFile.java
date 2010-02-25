package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;

public class NDFFile extends RowsAndColumns {
    public NDFFile (String filename, 
                         RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(filename,handler);
    }
    public NDFFile(InputStreamReader isr,
                        RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }        

}
