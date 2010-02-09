package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.util.*;

public class AFEHandler extends RowsColumnsHandler {
    private Map<String,String> headers;
    private String paramscols[];    
    public AFEHandler () {
        super();
        headers = new HashMap<String,String>();        
        paramscols = null;
    }
    
    public void processHeaderLine(String line) {
        if (line.matches("^TYPE.*")) {return;}
        if (line.matches("^\\*.*")) {return;}
        if (line.matches("^FEPARAMS.*")) {
            paramscols = line.split("\\t",-1);
            return;
        } 
        if (line.matches("^DATA.*") && 
            (paramscols != null)) {
            String vals[] = line.split("\\t",-1);
            for (int i = 1; i < vals.length; i++) {
                headers.put(paramscols[i],
                            vals[i]);
            }
            paramscols = null;
        }
    }
    public String getHeader(String key) {return headers.get(key);}
}
