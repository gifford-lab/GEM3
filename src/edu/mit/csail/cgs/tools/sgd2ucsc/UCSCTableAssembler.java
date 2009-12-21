/*
 * Created on Jan 29, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.sgd2ucsc;

import java.util.*;
import java.io.*;
import java.sql.*;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.StrandedFeature;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;

public class UCSCTableAssembler {
    
    public static void main(String[] args) { 
        File gffFile = args.length > 0 ? new File(args[0]) : 
            new File("C:\\Documents and Settings\\tdanford\\Desktop\\sacCer1.gff");
        File dbxrefFile = args.length > 0 ? new File(args[1]) : 
            new File("C:\\Documents and Settings\\tdanford\\Desktop\\dbxref.tab");
        
        try {
            Genome genome = Organism.findGenome("SGDv1");
            UCSCTableAssembler uta = new UCSCTableAssembler(gffFile, dbxrefFile);

            //uta.printSGDLines(System.out);
            uta.insertIntoDB(genome);
            
            SgdToNameTable namer = new SgdToNameTable(gffFile);
            namer.populateTable(genome);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private static BinCalculator binCalc;
    private static Map<String,String> gfftype2ucsc;
    private static int[] romvals;

    static {
        romvals = null;
        binCalc = new BinCalculator();
        gfftype2ucsc = new HashMap<String,String>();
        gfftype2ucsc.put("telomere", "Telomeric Region");
        gfftype2ucsc.put("tRNA", "tRNA");
        gfftype2ucsc.put("snoRNA", "snoRNA");
        gfftype2ucsc.put("snRNA", "snRNA");
        gfftype2ucsc.put("rRNA", "rRNA");
        gfftype2ucsc.put("ncRNA", "RNA");
        gfftype2ucsc.put("centromere", "CEN");
        gfftype2ucsc.put("transposable_element", "Transposon");
    }

    private SGDGFFParser gffParser;
    private UniProtGeneMap uniprotMap;

    private LinkedList<SGDOther> sgdOthers;
    private LinkedList<SGDGene> sgdGenes;
    
    public UCSCTableAssembler(File gff, File dbxref) throws IOException { 
        uniprotMap = new UniProtGeneMap(dbxref);
        gffParser = new SGDGFFParser();
        gffParser.parseInputFile(gff);
        sgdOthers = new LinkedList<SGDOther>();
        sgdGenes = new LinkedList<SGDGene>();
        populateSGDLists();
    }
    
    public void printSGDLines(PrintStream ps) { 
        for(SGDGene gene : sgdGenes) { ps.println(gene.toString()); }
        for(SGDOther other : sgdOthers) { ps.println(other.toString()); }
    }
    
    public void insertIntoDB(Genome g) { 
        try {
            java.sql.Connection cxn = g.getUcscConnection();
            insertIntoDB(cxn);
            DatabaseFactory.freeConnection(cxn);
            
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public void insertIntoDB(Connection cxn) throws SQLException {
        cxn.setAutoCommit(false);
        
        Statement s = cxn.createStatement();
        s.executeUpdate("delete from sgdGene");
        s.executeUpdate("delete from sgdOther");
        cxn.commit();
        s.close();
        
        PreparedStatement genePS = prepareSGDGeneInsert(cxn);
        PreparedStatement otherPS = prepareSGDOtherInsert(cxn);
        
        for(SGDGene gene : sgdGenes) { gene.insertIntoDB(genePS); }
        System.out.println("Inserted " + sgdGenes.size() + " sgdGene Entries.");
        cxn.commit();
        
        for(SGDOther other : sgdOthers) { other.insertIntoDB(otherPS); }
        System.out.println("Inserted " + sgdOthers.size() + " sgdOther Entries.");
        cxn.commit();
        
        genePS.close();
        otherPS.close();
        
        cxn.setAutoCommit(true);
    }
    
    public void populateSGDLists() { 
        sgdGenes.clear();
        sgdOthers.clear();
        
        LinkedList<SGDGene> exonGenes = new LinkedList<SGDGene>();
        
        for(String geneKey : gffParser.geneFeatures.keySet()) { 
            SGDGFFParser.GeneFeatures gf = gffParser.geneFeatures.get(geneKey);
            
            if(gf.gene.getFeature().equals("transposable_element_gene")) { 
                String id = gf.geneAttrs.get("ID").get(0);
                String uniprot = uniprotMap.containsGeneName(id) ? uniprotMap.getUniprot(id) : "n/a";
                SGDGene g = new SGDGene(gf, uniprot);
                
                if(g.exonCount > 0) { exonGenes.addLast(g); }
                sgdGenes.addLast(g);
                //System.out.println(g);
                
            } else if(gf.gene.getFeature().equals("pseudogene")) {                 
                SGDOther ot = new SGDOther(gf);
                
                sgdOthers.addLast(ot);
                //System.out.println(ot);
                
            } else { 
                String classification = gf.geneAttrs.get("orf_classification").get(0);
                
                if(classification.equals("Verified")) {
                    String id = gf.geneAttrs.get("ID").get(0);
                    String uniprot = uniprotMap.containsGeneName(id) ? uniprotMap.getUniprot(id) : "n/a";
                    SGDGene g = new SGDGene(gf, uniprot);
                    
                    if(g.exonCount > 0) { exonGenes.addLast(g); }
                    sgdGenes.addLast(g);
                    //System.out.println(g);
                    
                } else { 
                    SGDOther ot = new SGDOther(gf);
                    
                    sgdOthers.addLast(ot);
                    //System.out.println(ot);
                }
            }
        }
        
        HashSet<String> otherTypes = new HashSet<String>();
        HashSet<String> unrecognized = new HashSet<String>();
        
        for(GFFRecord gff : gffParser.otherRecords) {
            if(gfftype2ucsc.containsKey(gff.getFeature())) { 
                SGDOther other = new SGDOther(gff);
                otherTypes.add(gff.getFeature());
                
                sgdOthers.addLast(other);
                //System.out.println(other);
            } else { 
                //System.err.println("Unrecognized Other type: \"" + gff.getFeature() + "\"");
                unrecognized.add(gff.getFeature());
            }
        }
        
        System.out.println("# Other records: " + gffParser.otherRecords.size());
        System.out.println("\tRecognized Other Types: " + otherTypes);
        System.out.println("\tUnrecognized Other Types: " + unrecognized);
        
        System.out.println("# Exon-Genes: " + exonGenes.size());
    }
    
    public static String fixChrom(String chrom) {
        if (romvals == null) {
            romvals = new int[Character.getNumericValue('Z')];
            romvals[Character.getNumericValue('X')] = 10;
            romvals[Character.getNumericValue('V')] = 5;
            romvals[Character.getNumericValue('I')] = 1;
        }
        String chr = chrom;
        if (chr.matches("^[cC][hH][rR].*")) {
            chr = chr.substring(3);
        }
        if(chr.matches("^Mito$")) { return "M"; }
        if(chr.matches("^2-micron$")) { return "2micron"; }
        
         if(chr.matches("^[XVI]+$")) {
             int val = 0, pos = 1, curval, lastval, buffer; 
             char cur, last;
             boolean random = false;
             if (chr.matches("_random$")) {
                 random = true;
                 chr.replaceFirst("_random$","");
             }            
             last = chr.charAt(0);
             lastval = romvals[Character.getNumericValue(last)];
             buffer = lastval;
             while (pos < chr.length()) {
                 cur = chr.charAt(pos);
                 curval = romvals[Character.getNumericValue(cur)];
                 if (curval > lastval) {
                     val += curval - lastval;
                     buffer = 0;
                 } else if (cur != last) {
                     val += buffer;
                     buffer = curval;
                 } else {
                     buffer += curval;
                 }
                 last = cur;
                 lastval = curval;
                 pos++;
             }
             val += buffer;
             if (random) {
                 return Integer.toString(val) + "_random";
             } else {
                 return Integer.toString(val);
             }
         } else 
        if (chr.matches("^[1234567890MUXY]+(_random)?[LRh]?$")) {
            return chr;
        } else {
            throw new NumberFormatException("Can't fix chrom name " + chrom + "," + chr);
        }
    }

    public static String translateChromName(String seqName) { 
        return "chr" + fixChrom(seqName);
    }
    
    public static String translateStrand(boolean strand) { 
        return strand ? "+" : "-";
    }
    
    public static PreparedStatement prepareSGDOtherInsert(Connection cxn) throws SQLException {
        return cxn.prepareStatement("insert into sgdOther (bin, chrom, chromStart, chromEnd, " +
                "name, score, strand, type) values (?, ?, ?, ?, ?, 0, ?, ?)");
    }
    
    public static class SGDOther {
        
        public int bin;
        public String chrom;
        public int chromStart, chromEnd;
        public String name;
        public boolean strand;
        public String type;
        
        public void insertIntoDB(PreparedStatement ps) throws SQLException { 
            ps.setInt(1, bin);
            ps.setString(2, chrom);
            ps.setInt(3, chromStart);
            ps.setInt(4, chromEnd);
            ps.setString(5, name);
            ps.setString(6, translateStrand(strand));
            ps.setString(7, type);
            
            ps.executeUpdate();
        }
        
        public SGDOther(SGDGFFParser.GeneFeatures gf) { 
            if(gf.gene.getFeature().equals("pseudogene")) { 
                type = "CDS:pseudogene"; 
            } else { 
                type = "Dubious:CDS";
            }
            bin = binCalc.getBinFromRange(gf.gene.getStart(), gf.gene.getEnd());
            chrom = translateChromName(gf.gene.getSeqName());
            chromStart = gf.gene.getStart();
            chromEnd = gf.gene.getEnd();
            name = gf.geneAttrs.get("ID").get(0);
            strand = gf.gene.getStrand() == StrandedFeature.POSITIVE;
        }
        
        public SGDOther(GFFRecord gff) {
            Map<String,List<String>> attrs = SGDGFFParser.decodeAttrMap(gff);
            bin = binCalc.getBinFromRange(gff.getStart(), gff.getEnd());
            chrom = translateChromName(gff.getSeqName());
            chromStart = gff.getStart();
            chromEnd = gff.getEnd();
            name = attrs.get("Name").get(0);
            strand = gff.getStrand() == StrandedFeature.POSITIVE;
            
            if(!gfftype2ucsc.containsKey(gff.getFeature())) { throw new IllegalArgumentException(gff.getFeature()); }
            type = gfftype2ucsc.get(gff.getFeature());
        }
        
        public String toString() { 
            StringBuilder sb = new StringBuilder();
            sb.append("SGDOther");
            sb.append(" " + chrom + ":" + chromStart + "-" + chromEnd);
            sb.append(" " + name);
            sb.append(" " + translateStrand(strand));
            sb.append(" \"" + type + "\"");
            return sb.toString();
        }
    }
    
    public static PreparedStatement prepareSGDGeneInsert(Connection cxn) throws SQLException { 
        return cxn.prepareStatement("insert into sgdGene (name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, " +
                "exonCount, exonStarts, exonEnds, proteinID) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    }
    
    public static class SGDGene { 

        public String name, chrom;
        public boolean strand;
        public int txStart, txEnd, cdsStart, cdsEnd;
        public int exonCount;
        public Vector<Integer> exonStarts, exonEnds;
        public String protID;
        
        public void insertIntoDB(PreparedStatement ps) throws SQLException { 
            ps.setString(1, name);
            ps.setString(2, chrom);
            ps.setString(3, translateStrand(strand));
            ps.setInt(4, txStart);
            ps.setInt(5, txEnd);
            ps.setInt(6, cdsStart);
            ps.setInt(7, cdsEnd);
            ps.setInt(8, exonCount);
            
            String ess = getExonStartsString();
            Reader esreader = new StringReader(ess);
            ps.setCharacterStream(9, esreader, ess.length());
            
            String ees = getExonEndsString();
            Reader eereader = new StringReader(ees);
            ps.setCharacterStream(10, eereader, ees.length());
            
            ps.setString(11, protID);
            
            ps.executeUpdate();
        }
        
        public SGDGene(SGDGFFParser.GeneFeatures gf, String uniprot) {
            name = gf.geneAttrs.get("ID").get(0);
            chrom = translateChromName(gf.gene.getSeqName());
            txStart = gf.gene.getStart();
            txEnd = gf.gene.getEnd();
            protID = uniprot;
            exonStarts = new Vector<Integer>();
            exonEnds = new Vector<Integer>();
            strand = gf.gene.getStrand() == StrandedFeature.POSITIVE;
            
            exonCount = 0; 
            
            if(gf.cds.size() > 0) { 
                cdsStart = txEnd;
                cdsEnd = txStart;
                for(GFFRecord cds : gf.cds) { 
                    int es = cds.getStart(), ee = cds.getEnd();
                    exonStarts.add(es); exonEnds.add(ee);
                    exonCount++;
                    cdsStart = Math.min(cdsStart, es);
                    cdsEnd = Math.max(cdsEnd, ee);
                }         
            } else { 
                cdsStart = txStart; cdsEnd = txEnd;
            }
        }
        
        public String getExonStartsString() { 
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < exonCount; i++) { sb.append(exonStarts.get(i) + ","); }
            return sb.toString();
        }
        
        public String getExonEndsString() { 
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < exonCount; i++) { sb.append(exonEnds.get(i) + ","); }
            return sb.toString();
        }
        
        public String toString() { 
            StringBuilder sb = new StringBuilder();
            sb.append("SGDGene");
            sb.append(" " + chrom + ":" + txStart + "-" + txEnd + " " + cdsStart + ":" + cdsEnd);
            sb.append(" " + name + "/" + protID);
            sb.append(" " + translateStrand(strand));
            sb.append(" " + exonCount);
            for(int i = 0; i < exonCount; i++) { sb.append(" " + exonStarts.get(i) + "-" + exonEnds.get(i)); }
            return sb.toString();
        }
    }
}
