package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;

public class PairFile extends PairedRowsAndColumns {

    public PairFile(String fone, String ftwo, PairedRowsColumnsHandler handler) throws IOException, FileNotFoundException {
        super(fone,ftwo,handler);        
    }
    public PairFile(InputStreamReader isrone,
                    InputStreamReader isrtwo,
                    PairedRowsColumnsHandler handler) throws IOException {
        super(isrone, isrtwo,handler);
    }
    public void processFileOpen() throws IOException {
        processHeaderLine(getReaderOne().readLine(),
                          getReaderTwo().readLine());
    }
    public void processHeaderLine(String lineone, String linetwo) {
        getHandler().processHeaderLine(lineone,linetwo);
    }
    

}
