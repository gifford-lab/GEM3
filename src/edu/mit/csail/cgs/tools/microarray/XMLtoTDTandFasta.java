package edu.mit.csail.cgs.tools.microarray;

import java.io.*;
import java.net.*;
import java.util.*;
import org.xml.sax.*;
import org.xml.sax.helpers.*;
import edu.mit.csail.cgs.utils.xml.XMLTreeHandler;
import edu.mit.csail.cgs.utils.xml.XMLTreeElement;

public class XMLtoTDTandFasta extends DefaultHandler {

    private PrintWriter fasta, tdt;
    ProbeInfo lastprobe;
    private List<ProbeInfo> probes;
    private Set<Double> x, y;

    public XMLtoTDTandFasta() {
        super();
    }

    public void handleFile(String fname) throws IOException, SAXException {
        String tdtname = fname.replaceAll(".xml$",".tdt");
        String fastaname = fname.replaceAll(".xml$",".fasta");
        XMLReader oParser = XMLReaderFactory.createXMLReader(); 
        oParser.setContentHandler(this);
        EntityResolver resolver = new EntityResolver() {
                public InputSource resolveEntity(String publicId, String systemId) {                    
                    System.out.println("Resolve: \"" + publicId + "\", \"" + systemId + "\"");                                           
                    //String dtdPath = "edu.mit.csail.cgs.tools.chipchip/";
                    String dtdPath = "edu/mit/csail/cgs/tools/chipchip/";
                    String[] recs = { "GEMLPattern.dtd"};
                    String[] paths = recs;
                    for (int i = 0; i < recs.length; i++) { 
                        if (systemId.endsWith(recs[i])) { 
                            URL url = ClassLoader.getSystemClassLoader().getResource(dtdPath + paths[i]);
                            System.err.println(systemId + " ends with " + recs[i] + " url is " + url);
                            System.err.println("path was " + dtdPath + paths[i]);
                            try {
                                return new InputSource(url.openStream());
                            } catch (IOException e) {
                                System.err.println("Couldn't open " + url);
                                e.printStackTrace();
                            } catch (NullPointerException e) {
                                System.err.println("Doesn't exist " + url);
                                e.printStackTrace(System.err);
                            }                     
                        }
                    }
                    System.out.println("No match.");
                    return null;
                }
            };
        oParser.setEntityResolver(resolver);
        fasta = new PrintWriter(fastaname);
        tdt = new PrintWriter(tdtname);
        probes = new ArrayList<ProbeInfo>();
        x = new HashSet<Double>();
        y = new HashSet<Double>();
        tdt.println("Row\tColumn\tID\tSequence\tName\tControlType");
        oParser.parse(new InputSource(new FileReader(fname)));        
        fasta.close();
        Map<Double,Integer> cols, rows;
        cols = getGrid(y);
        rows = getGrid(x);
        System.err.println("Got " + cols.size() + " cols and " + rows.size() + " rows");
        for (ProbeInfo p : probes) {
            tdt.println((rows.get(p.x) + 1)+ "\t" + (cols.get(p.y) + 1) + "\t" + p.id + "\t" + p.seq + "\t" + p.name + "\t" + p.type);
        }
        tdt.close();
    }

    public Map<Double,Integer> getGrid(Set<Double> coords) {
        List<Double> list = new ArrayList<Double>();
        list.addAll(coords);
        Collections.sort(list);
        HashMap<Double,Integer> out = new HashMap<Double,Integer>();
        for (int i = 0; i < coords.size(); i++) {
            out.put(list.get(i),i);
        }
        return out;
    }

    public void startElement(java.lang.String uri, java.lang.String localName, java.lang.String qName, Attributes attributes) {
        if (localName.equals("reporter")) {
            lastprobe = new ProbeInfo();
            lastprobe.id = attributes.getValue("name");
            lastprobe.name = attributes.getValue("systematic_name");
            lastprobe.seq = attributes.getValue("active_sequence");
            lastprobe.type = attributes.getValue("control_type");
            if (lastprobe.name == null) {
                lastprobe.name = "";
            }
            if (lastprobe.type == null) {
                lastprobe.type = "";
            }
            if (lastprobe.seq == null) {
                lastprobe.seq = "";
            }
            if (lastprobe.seq.length() > 0) {
                fasta.println(">" + lastprobe.id + "\n" + lastprobe.seq);
            }
            probes.add(lastprobe);
        }
        if (localName.equals("position")) {            
            Double xval = Double.parseDouble(attributes.getValue("x"));
            Double yval = Double.parseDouble(attributes.getValue("y"));
            x.add(xval);
            y.add(yval);
            lastprobe.x = xval;
            lastprobe.y = yval;
        }
    }

    public static void main(String args[]) throws Exception {
        XMLtoTDTandFasta converter = new XMLtoTDTandFasta();
        for (int i = 0; i < args.length; i++) {
            converter.handleFile(args[i]);
        }
    }    
}

class ProbeInfo {
    public String name, seq, type, id;
    public Double x, y;
}