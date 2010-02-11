package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;

public class GPRFile extends AxonTextFile {
    public GPRFile (String filename, 
                         RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(filename,handler);
    }
    public GPRFile(InputStreamReader isr,
                        RowsColumnsHandler handler) throws IOException {
        super(isr,handler);
    }
    
    public String[] fixColumnLabels(String[] fields) {
        for (int i = 0; i < fields.length; i++) {
            fields[i] = fields[i].replaceAll("650","635");
            fields[i] = fields[i].replaceAll("550","532");
            fields[i] = fields[i].replaceAll("\\s\\(635\\/532\\)","");
            //            fields[i] = "Sequence";
        }
        return fields;
    }

}
