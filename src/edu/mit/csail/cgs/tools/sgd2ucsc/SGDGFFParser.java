/*
 * Created on Jan 26, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.sgd2ucsc;

import java.util.*;
import java.net.URLDecoder;
import java.sql.*;
import java.io.*;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.*;
import org.biojava.utils.ParserException;
import org.biojava.bio.seq.StrandedFeature;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Pair;

public class SGDGFFParser {

    public static void main(String[] args) { 
        File inputFile = args.length > 0 ? 
                new File(args[0]) : 
                new File("C:\\Documents and Settings\\tdanford\\Desktop\\sacCer1.gff");
                
        try {
            SGDGFFParser parser = new SGDGFFParser();
            BinCalculator bincalc = new BinCalculator();
            parser.parseInputFile(inputFile);
            
            for(String k : parser.geneFeatures.keySet()) { 
                try {
                    GeneFeatures gf = parser.geneFeatures.get(k);
                    StringBuffer line = new StringBuffer();
                    ArrayList<Integer> starts = new ArrayList<Integer>();
                    ArrayList<Integer> ends = new ArrayList<Integer>();
                    for (GFFRecord cds : gf.cds) {
                        starts.add(cds.getStart());
                        ends.add(cds.getEnd());
                    }
                    if (starts.size() == 0) {
                        starts.add(gf.gene.getStart());
                        ends.add(gf.gene.getEnd());
                    }

                    Collections.sort(starts);
                    Collections.sort(ends);
                    String strand = "";
                    if (gf.gene.getStrand() == StrandedFeature.NEGATIVE) {
                        strand = "-";
                    } else if (gf.gene.getStrand() == StrandedFeature.POSITIVE) {
                        strand = "+";
                    }
                    /* dumb that we have to put this back on, but that's the UCSC standard format */
                    String chrom = "chr" + Genome.fixYeastChrom(gf.gene.getSeqName());
               
                    line.append(gf.id + "\t");
                    line.append(chrom + "\t"); // do we need to get rid of roman numerals here ?
                    line.append(strand + "\t");
                    line.append(gf.gene.getStart() + "\t");
                    line.append(gf.gene.getEnd() + "\t");
                    line.append(starts.get(0) + "\t");
                    line.append(ends.get(ends.size()-1) + "\t");
                    line.append(starts.size() + "\t");
                    for (int i = 0; i < starts.size(); i++) {
                        if (i == 0) {
                            line.append(starts.get(i));
                        } else {
                            line.append("," + starts.get(i));
                        }
                    }
                    line.append("\t");
                    for (int i = 0; i < ends.size(); i++) {
                        if (i == 0) {
                            line.append(ends.get(i));
                        } else {
                            line.append("," + ends.get(i));
                        }
                    }
                    line.append("\t");
                
                    // no protein ID
                    System.out.println(line);
                } catch (Exception e) {
                    System.err.println(e.toString());
                    e.printStackTrace();
                }
            }
            for (GFFRecord record : parser.otherRecords) {
                if (record.getFeature().equals("intron")) {
                    continue;
                }

                try {
                    StringBuffer line = new StringBuffer();
                    line.append(bincalc.getBinFromRange(record.getStart(),
                                                        record.getEnd()) + "\t");
                    String strand = "";
                    if (record.getStrand() == StrandedFeature.NEGATIVE) {
                        strand = "-";
                    } else if (record.getStrand() == StrandedFeature.POSITIVE) {
                        strand = "+";
                    }
                    /* dumb that we have to put this back on, but that's the UCSC standard format */
                    String chrom = "chr" + Genome.fixYeastChrom(record.getSeqName());
                    Map<String,List<String>> attrs = decodeAttrMap(record);
                    String name = "";
                    if (attrs != null && attrs.containsKey("Name") && attrs.get("Name").size() > 0) {
                        name = attrs.get("Name").get(0);
                    } else if (attrs != null && attrs.containsKey("ID") && attrs.get("ID").size() > 0) {
                        name = attrs.get("ID").get(0);
                    }



                    line.append(chrom + "\t");
                    line.append(record.getStart() + "\t" + record.getEnd() + "\t");
                    line.append(name + "\t");
                    line.append(record.getScore() + "\t");
                    line.append(strand + "\t");
                    line.append(record.getFeature());
                    System.out.println(line);
                } catch (Exception e) {                    
                    System.err.println(e.toString());
                    e.printStackTrace();
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public static void printKeyValues(String k, List<String> v, PrintStream ps) { 
        ps.print(k + " => ");
        boolean first = true;
        for(String val : v) { 
            if(!first) { ps.print(", "); }
            first = false;
            ps.print(val);
        }
        ps.println();
    }
    
    public Map<String,GeneFeatures> geneFeatures;
    public LinkedList<GFFRecord> otherRecords;
    
    public SGDGFFParser() { 
        geneFeatures = new HashMap<String,GeneFeatures>();
        otherRecords = new LinkedList<GFFRecord>();
    }
    
    public void parseInputFile(File inputFile) throws IOException {
        geneFeatures.clear();
        otherRecords.clear();

        try { 
            GFFEntrySet gffEntries = GFFTools.readGFF(inputFile);

            Iterator itr = gffEntries.lineIterator();
            int count = 0;
            int intronFeatures = 0;
            LinkedList<GFFRecord> cdsRecs = new LinkedList<GFFRecord>();

            while(itr.hasNext()) { 
                Object val = itr.next();
                if(val instanceof GFFRecord) { 
                    GFFRecord rec = (GFFRecord)val;
                    count += 1;

                    if(rec.getFeature().endsWith("gene")) { 
                        GeneFeatures gf = new GeneFeatures(rec);
                        geneFeatures.put(gf.id, gf);
                    } else if(rec.getFeature().equals("CDS")) {
                        cdsRecs.addLast(rec);
                    } else { 
                        otherRecords.add(rec);
                    }
                }
            }

            for(GFFRecord rec : cdsRecs) { 
                Map<String,List<String>> attrs = decodeAttrMap(rec);
                if(geneFeatures.containsKey(attrs.get("Parent").get(0))) { 
                    geneFeatures.get(attrs.get("Parent").get(0)).addCDS(rec, attrs);                                
                } else { 
                    System.err.println("Unknown CDS Parent: " + attrs.get("Parent").get(0));
                }                
            }

            for(String k : geneFeatures.keySet()) { 
                GeneFeatures gf = geneFeatures.get(k);
                if(gf.cds != null && gf.cds.size() > 1) { 
                    intronFeatures++;
                }
            }

            System.err.println("# GFF Records: " + count);
            System.err.println("# Gene Feature Sets: " + geneFeatures.size());
            System.err.println("\t# Intron-Features: " + intronFeatures);

        } catch (ParserException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }        
    }
    
    public static Pair<String,List<String>> decodeKeyValues(String str) { 
        int index = str.indexOf("=");
        if(index == -1) { return null; }
        String k = str.substring(0, index);
        LinkedList<String> vlist = new LinkedList<String>();
        try {
            String vstr = URLDecoder.decode(str.substring(index+1, str.length()), "UTF-8");
            String[] array = vstr.split(",");
            for(int i = 0; i < array.length; i++) { vlist.addLast(array[i]); }
            return new Pair<String,List<String>>(k, vlist);
            
        } catch (UnsupportedEncodingException e) {
            System.err.println("BAD STRING " + str);
            e.printStackTrace();
            return null;
        } catch (Exception e) {
            System.err.println("BAD STRING " + str);
            return null;
        }

    }
    
    public static Map<String,List<String>> decodeAttrMap(GFFRecord rec) { 
        Map<String,List<String>> map = new HashMap<String,List<String>>();
        for(Object ko : rec.getGroupAttributes().keySet()) { 
            String str = (String)ko;
            Pair<String,List<String>> vals = decodeKeyValues(str);
            if(vals != null) { 
                map.put(vals.getFirst(), vals.getLast());
            }
        }
        return map;
    }
    
    public static class GeneFeatures { 
        
        public String id;
        public GFFRecord gene;
        public Map<String,List<String>> geneAttrs;
        public LinkedList<GFFRecord> cds;
        
        public GeneFeatures(GFFRecord g) {
            if(!g.getFeature().endsWith("gene")) { throw new IllegalArgumentException(); }
            gene = g;
            geneAttrs = decodeAttrMap(g);
            cds = new LinkedList<GFFRecord>();
            if(!geneAttrs.containsKey("ID")) { throw new IllegalArgumentException(); }
            id = geneAttrs.get("ID").get(0);
        }
        
        public void addCDS(GFFRecord rec, Map<String,List<String>> attrs) { 
            if(!rec.getFeature().equals("CDS")) { throw new IllegalArgumentException(); }
            if(attrs == null) { attrs = decodeAttrMap(rec); } 
            if(!attrs.containsKey("Parent")) { throw new IllegalArgumentException(); }
            if(!attrs.get("Parent").contains(id)) { throw new IllegalArgumentException(); }
            
            cds.addLast(rec);
        }        
    }
    
}
