package edu.mit.csail.cgs.tools.microarray;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.utils.parsing.textfiles.*;

/* takes name of input file(s) on command line.
   Produces FASTA formatted output on STDOUT */

public class SeqFromNDF extends NDFHandler {
    public static void main(String args[]) throws Exception {
        for (String fname : args) {
            SeqFromNDF handler = new SeqFromNDF();
            NDFFile file = new NDFFile(fname,handler);
            file.parse();
        }

    }

    private int idcol, seqcol;
    public SeqFromNDF() {
        super();
        idcol = -1;
        seqcol = -1;
    }
    public void setIDCol(int i) {idcol = i;}
    public void setSequenceCol(int i) {seqcol = i;}
    public String getIDLabel() {return "PROBE_ID";}
    public String getSequenceLabel() {return "PROBE_SEQUENCE";}

    public void processLabelLine(String[] headers) {
        super.processLabelLine(headers);
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals(getIDLabel())) {setIDCol(i);}
            if (headers[i].equals(getSequenceLabel())) {setSequenceCol(i);}
        }
        if (idcol == -1) {throw new RuntimeException("Couldn't get ID in " + headers);}
        if (seqcol == -1) {throw new RuntimeException("Couldn't get Sequence in " + headers);}
    }

    public void processBodyLine(String[] fields) {
        System.out.println(">" + fields[idcol] + "\n" + fields[seqcol]);
    }


}
