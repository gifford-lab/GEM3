package edu.mit.csail.cgs.utils.parsing.textfiles;

import java.io.*;

public class GALFile extends AxonTextFile {
    public GALFile (String filename, 
                    RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(filename,handler);
    }
    public GALFile(InputStreamReader isr,
                        RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }
    
    public String[] fixColumnLabels(String[] fields) {
        for (int i = 0; i < fields.length; i++) {
            if (fields[i].matches("^[Dd]escription$")) {
                fields[i] = "Sequence";
            }
        }
        return fields;
    }

}
