/*
 * Created on Mar 1, 2007
 */
package edu.mit.csail.cgs.utils.xml;

import java.util.HashSet;
import java.util.Set;
import java.util.Stack;

import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;

public class XMLTreeHandler extends DefaultHandler {

    private Set<String> tags;
    private Stack<StringBuilder> charBuilders;
    private Stack<Boolean> immediate;
    private XMLTreeElement root, currentElement;
    
    public XMLTreeHandler()  {
        super();
        tags = new HashSet<String>();
        charBuilders = new Stack<StringBuilder>();
        root = new XMLTreeElement("#root#", null);
        currentElement = root;
        immediate = new Stack<Boolean>();
        immediate.push(true);
    }
    
    public XMLTreeElement getTree() { return root; }
    public void addTag(String tag) { tags.add(tag); }
    
    public void characters(char[] array, int start, int length) {
        if(!charBuilders.isEmpty()) {
            StringBuilder sb = charBuilders.peek();
            for(int i = 0; i < length; i++ ) { sb.append(array[start+i]); }
        }
    }

    public void startElement(String uri, String localName, String qName, Attributes atts) {
        if(tags.contains(localName)) { 
            XMLTreeElement newCurrent = new XMLTreeElement(localName, currentElement);
            currentElement = newCurrent;
            immediate.push(true);
        } else { 
            immediate.push(false); 
        }
        
        charBuilders.push(new StringBuilder());
    }

    public void endElement(String uri, String localName, String qName) {
        if(currentElement != null) { 
            if(currentElement.name.equals(localName)) { 
                if(!charBuilders.isEmpty()) { 
                    currentElement.content = charBuilders.peek().toString(); 
                }
                currentElement = currentElement.parent;
            } else {
                if(immediate.get(0)) { 
                    String value = !charBuilders.isEmpty() ? charBuilders.peek().toString() : "";
                    String key = localName;
                    currentElement.values.put(key, value);
                }
            }
        }
        charBuilders.pop();
        immediate.pop();
    }
}
