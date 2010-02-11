package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;

/* Base class for parsing a file with rows and columns */

public abstract class RowsAndColumns {

    private BufferedReader reader;
    private RowsColumnsHandler handler;

    public RowsAndColumns (String filename, 
                           RowsColumnsHandler handler) throws IOException, FileNotFoundException {
        reader = new BufferedReader(new FileReader(filename));
        this.handler = handler;
    }

    public RowsAndColumns(InputStreamReader isr,
                          RowsColumnsHandler handler) throws IOException {
        reader = new BufferedReader(isr);
        this.handler = handler;
    }

    public BufferedReader getReader() {return reader;}
    public RowsColumnsHandler getHandler() {return handler;}

    public void parse() throws IOException {
        processFileOpen();
        String line = reader.readLine();
        processLabelLine(line);
        line = reader.readLine();
        while (line != null) {
            processBodyLine(line);
            line = reader.readLine();
        }
        processFileClose();
        reader.close();
    }

    public void processFileOpen() throws IOException { }
    public void processHeaderLine(String headerline) {
        handler.processHeaderLine(headerline);
    }
    public void processLabelLine(String labelline) {
        String[] fields = labelline.split("\\t",-1);
        for (int i = 0; i < fields.length; i++) {
            fields[i] = fields[i].replaceAll("^\"","").replaceAll("\"$","");
        }
        fields = fixColumnLabels(fields);
        handler.processLabelLine(fields);
    }
    public String[] fixColumnLabels(String[] fields) {
        return fields;
    }
    public void processBodyLine(String bodyline) {
        String[] fields = bodyline.split("\\t",-1);
        for (int i = 0; i < fields.length; i++) {
            fields[i] = fields[i].replaceAll("^\"","").replaceAll("\"$","");
        }
        handler.processBodyLine(fields);
    }
    public void processFileClose() throws IOException {
        handler.processFileClose();
    }
}
