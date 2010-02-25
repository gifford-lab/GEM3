package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.io.*;
import java.util.Map;
import java.util.HashMap;

public class PairHandler extends PairedRowsColumnsHandler {

    public Map<String,String> headersone, headerstwo;

    public PairHandler() {
        headersone = new HashMap<String,String>();
        headerstwo = new HashMap<String,String>();
    }


    public void processHeaderLine(String lineone, String linetwo) {
        lineone = lineone.replaceAll("^#\\s*","");
        String[] pieces = lineone.split("\\t");
        for (int i = 0; i < pieces.length; i++) {
            String[] kv = pieces[i].split("=",2);
            if (kv.length == 2) {
                headersone.put(kv[0],kv[1]);
            }
        }
        pieces = linetwo.split("\\t");
        for (int i = 0; i < pieces.length; i++) {
            String[] kv = pieces[i].split("=",2);
            if (kv.length == 2) {
                headerstwo.put(kv[0],kv[1]);
            }
        }
    }
    public String getHeaderOne(String key) {return headersone.get(key);}
    public String getHeaderTwo(String key) {return headerstwo.get(key);}


}
