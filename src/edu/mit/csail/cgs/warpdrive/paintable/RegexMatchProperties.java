package edu.mit.csail.cgs.warpdrive.paintable;

import java.util.Map;

public class RegexMatchProperties extends PaintableProperties {
    public String[] expressions;
    public String[] labels;

    public RegexMatchProperties (int n) {
        expressions = new String[n];
        labels = new String[n];
    }
    public RegexMatchProperties(Map<String,String> map) {
        int n = map.size();
        expressions = new String[n];
        labels = new String[n];
        int i = 0;
        for (String label : map.keySet()) {
            labels[i] = label;
            expressions[i] = map.get(label);
            i++;
        }
    }
    public void addRegex(String regex) {
        addRegex(regex,regex);
    }
    public void addRegex(String label, String regex) {
        int i = firstFree();
        if (i < 0) {
            expand();
            i = firstFree();
        }
        expressions[i] = regex;
        labels[i] = label;
    }
    private int firstFree() {
        for (int i = 0; i < expressions.length; i++) {
            if (expressions[i] == null) {
                return i;
            }
        }
        return -1;
    }
    private void expand() {
        String[] newe = new String[expressions.length * 2];
        String[] newl = new String[labels.length * 2];
        for (int i = 0; i < expressions.length; i++) {
            newe[i] = expressions[i];
            newl[i] = labels[i];
        }
        expressions = newe;
        labels = newl;
    }
}