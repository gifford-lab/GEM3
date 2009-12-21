package edu.mit.csail.cgs.utils.parsing.textfiles;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;

public class TDTFile extends RowsAndColumns {
    public TDTFile (String filename, 
                    RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(filename,handler);
    }

    public TDTFile(InputStreamReader isr,
                   RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }

    public void processFileOpen() throws IOException {
        super.processFileOpen();
    }

}
