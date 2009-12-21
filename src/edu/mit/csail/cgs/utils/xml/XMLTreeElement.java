/*
 * Created on Mar 1, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.utils.xml;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Vector;

public class XMLTreeElement {
    
    private static final String istr = "  ";
    
    public XMLTreeElement parent;
    public String name;
    public Map<String,String> values;
    public String content;
    
    public Vector<XMLTreeElement> children;
    
    public XMLTreeElement(String n, XMLTreeElement p) { 
        parent = p;
        name = n;
        values = new LinkedHashMap<String,String>();
        content = null;
        children = new Vector<XMLTreeElement>();

        if(parent != null) { parent.children.add(this); }
    }
    
    public LinkedList<XMLTreeElement> collectElements(String n) { 
        LinkedList<XMLTreeElement> elmts = new LinkedList<XMLTreeElement>();
        gatherElements(n, elmts);
        return elmts;
    }
    
    private void gatherElements(String n, LinkedList<XMLTreeElement> elmts) { 
        if(name.equals(n)) { elmts.addLast(this); }
        for(XMLTreeElement child : children) { child.gatherElements(n, elmts); }
    }
    
    public void printTree(PrintStream ps) { 
        printTree(ps, 0);
    }
    
    private void printTree(PrintStream ps, int indent) { 
        StringBuilder idb = new StringBuilder();
        for(int i = 0; i < indent; i++) { idb.append(istr); }
        String indentStr = idb.toString();
        
        ps.println(indentStr + name + " {" + content + "} --> " + children.size() + " children.");
        for(String key : values.keySet()) { 
            ps.println(indentStr + istr + key + "=>\"" + values.get(key) + "\"");
        }
        for(int i = 0; i < children.size(); i++) { 
            children.get(i).printTree(ps, indent+1);
        }
    }
}
