package edu.mit.csail.cgs.utils.parsing.textfiles;

import java.io.*;

public abstract class PairedRowsAndColumns {
    private BufferedReader readerone, readertwo;
    private PairedRowsColumnsHandler handler;

    public PairedRowsAndColumns (String filenameone, 
                                 String filenametwo,
                                 PairedRowsColumnsHandler handler) throws IOException, FileNotFoundException {
        readerone = new BufferedReader(new FileReader(filenameone));
        readertwo = new BufferedReader(new FileReader(filenametwo));
        this.handler = handler;
    }

    public PairedRowsAndColumns(InputStreamReader isrone,
                                InputStreamReader isrtwo,
                                PairedRowsColumnsHandler handler) throws IOException {
        readerone = new BufferedReader(isrone);
        readertwo = new BufferedReader(isrtwo);
        this.handler = handler;
    }

    public BufferedReader getReaderOne() {return readerone;}
    public BufferedReader getReaderTwo() {return readertwo;}
    public PairedRowsColumnsHandler getHandler() {return handler;}

    public void parse() throws IOException {
        processFileOpen();
        processLabelLine(readerone.readLine(),
                         readertwo.readLine());
        String lone, ltwo;
        lone = readerone.readLine();
        ltwo = readertwo.readLine();
        while (lone != null && ltwo != null) {
            processBodyLine(lone,ltwo);
            lone = readerone.readLine();
            ltwo = readertwo.readLine();            
        }
        processFileClose();
        readerone.close();
        readertwo.close();
    }

    public void processFileOpen() throws IOException {
        handler.processFileOpen();
    }
    
    public void processLabelLine(String lineone,
                                 String linetwo) {
        String[] fieldsone = lineone.split("\\t",-1);
        String[] fieldstwo = linetwo.split("\\t",-1);
        fieldsone = fixColumnLabels(fieldsone);
        fieldstwo = fixColumnLabels(fieldstwo);
        handler.processLabelLine(fieldsone, fieldstwo);
    }
    public String[] fixColumnLabels(String[] fields) {
        return fields;
    }
    public void processBodyLine(String lineone, String linetwo) {
        String[] fieldsone = lineone.split("\\t",-1);
        String[] fieldstwo = linetwo.split("\\t",-1);
        handler.processBodyLine(fieldsone,fieldstwo);
    }

    public void processFileClose() throws IOException {
        handler.processFileClose();
    }


}
