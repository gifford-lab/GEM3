package edu.mit.csail.cgs.tools.microarray;

import java.util.*;
import java.io.*;
import org.xml.sax.*;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.helpers.XMLReaderFactory;

public class XMLtoFASTA extends DefaultHandler {


    public void startElement(String uri, String localName, String qName, Attributes atts) {
        if (localName.equals("reporter")) {
            System.out.println(">" + atts.getValue("name") + "\n" + atts.getValue("active_sequence"));
        }
    }

    public static void main(String args[]) throws Exception {
        XMLtoFASTA x2f = new XMLtoFASTA();
        XMLReader oParser = XMLReaderFactory.createXMLReader(); 
        oParser.setContentHandler(x2f);

        for (String fname : args) {
            FileReader reader = new FileReader(fname);            
            InputSource source = new InputSource(reader);
            oParser.parse(source);
        }

    }



}
