package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.util.*;

public class ATFHandler extends RowsColumnsHandler {
    
    private Map<String,String> headers;
    public ATFHandler() {
        super();
        headers = new HashMap<String,String>();        
    }
    public Map<String,String> getHeaders() {return headers;}
    public String getHeader(String key) {return headers.get(key);}
    public void processHeaderLine(String line) {
        System.err.println("Parsing header line " + line);
        line.replaceAll("^\\\"","");
        line.replaceAll("\\\"$","");
        String f[] = line.split("=",-1);
        for (int i = 0; i < f.length; i++) {
            f[i].replaceAll("^\\s+","");
            f[i].replaceAll("\\s+$","");
        }
        if (f.length >= 2) {
            headers.put(f[0],f[1]);
        } else {
            System.err.println("Too few fields in " + line);
        }
    }
    
}
