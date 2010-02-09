package edu.mit.csail.cgs.utils.io.parsing.textfiles;

import java.util.Map;
import java.util.Set;
import java.util.HashMap;

public class RowsColumnsHandler {

    private Map<String,Integer> columns;
    public RowsColumnsHandler() {
        columns = new HashMap<String,Integer>();
    }
    
    public void processHeaderLine(String line) {}
    public void processLabelLine(String[] fields) {
        for (int i = 0; i < fields.length; i++) {
            columns.put(fields[i],i);
        }
    }
    public void processBodyLine(String[] fields) {}
    public void processFileClose() {}
    public Map<String,Integer> getColumns() {return columns;}
    public Integer getColumn(String c) {return columns.get(c);}
    public String joinStringArray(String[] a) {
        StringBuffer out = new StringBuffer();
        for (int i = 0; i < a.length; i++) {
            if (i > 0) {
                out.append("\t");
            }
            out.append(a[i]);
        }
        return out.toString();
    }
}
