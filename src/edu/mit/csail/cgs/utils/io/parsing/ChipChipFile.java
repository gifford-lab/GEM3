package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;
import java.io.*;

public abstract class ChipChipFile {
    private String dset, method, ip, wce, condition, chrom, type, filename, directory;

    public ChipChipFile(String filename, String type) {
        this.filename = filename;
        this.type = type;
        parseName();
    }
    public void parseName() {
        String[] pathpieces = filename.split(File.separator);
        String[] pieces = pathpieces[pathpieces.length-1].split("\\.");
        dset = pieces[0];
        method = pieces[1];
        chrom = pieces[2];
        int vs = method.indexOf(" vs ");
        int in = method.indexOf(" in ");
        ip = method.substring(0,vs);
        wce = method.substring(vs+4,in);
        condition = method.substring(in+4);
        if (filename.lastIndexOf(File.separator) == -1) {
            directory = ".";
        } else {
            directory = filename.substring(0,filename.lastIndexOf(File.separator));
        }
    }    
    public String getBaseFilename() {return filename;}
    public String getDset() {return dset;}
    public String getMethod() {return method;}
    public String getIP() {return ip;}
    public String getWCE() {return wce;}
    public String getCondition() {return condition;}
    public String getChrom() {return chrom;}
    public String getType() {return type;}
    public String getExpt() {return ip + " vs " + wce + " in " + condition;}
    public String getDirectory() {return directory;}
}
