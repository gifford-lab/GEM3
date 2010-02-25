/*
 * Created on Aug 14, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.homologene;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 * 
 * Built to parse the homologene.data file that comes with the 
 * HomoloGene data distribution.  
 */
public class HomoloGeneAssignments {

    public static void main(String[] args) {
        try {
            HomoloGeneAssignments assigns = new HomoloGeneAssignments(new File(args[0]));
            String line = null;
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

            System.out.print(">"); System.out.flush();
            while((line = br.readLine()) != null) {
                line = line.trim();
                if(line.length() > 0) { 
                    String[] array = line.split("\\s+");
                    String hg = null, tax = null, gene = null;
                    if(!array[0].equals("-")) { hg = array[0]; }
                    if(!array[1].equals("-")) { tax = array[1]; }
                    if(!array[2].equals("-")) { gene = array[2]; }
                    
                    Set<Entry> entries = assigns.selectEntries(hg, tax, gene);
                    for(Entry e : entries) { 
                        System.out.println(e.toString());
                    }
                    System.out.println("# Entries: " + entries.size());
                    System.out.println();
                }
                
                System.out.print(">"); System.out.flush();
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }        
    }
    
    private Vector<Entry> entries;

    public HomoloGeneAssignments(File f) throws IOException { 
        parse(f);
    }
    
    public HomoloGeneAssignments() throws IOException { 
        parse(new File("homologene.data"));
    }
    
    public Set<Entry> selectEntries(String hgID, String taxID, String geneID) { 
        HashSet<Entry> selected = new HashSet<Entry>();
        
        for(Entry e : entries) { 
            boolean isSelected = true;
            if(isSelected && hgID != null && !e.homologeneID.equals(hgID)) { isSelected = false; }
            if(isSelected && taxID != null && !e.taxID.equals(taxID)) { isSelected = false; }
            if(isSelected && geneID != null && !e.geneID.equals(geneID)) { isSelected = false; }
            
            if(isSelected) { selected.add(e); }
        }
        
        return selected;
    }
    
    private void parse(File f) throws IOException { 
        entries = new Vector<Entry>();
        
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = null;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                Entry e = new Entry(line);
                entries.add(e);
            }
        }
        
        br.close();
    }
    
    public static class Entry { 
        public String homologeneID;
        public String taxID;
        public String geneID, geneSymb;
        public String protein, proteinAcc;
        
        // note to self: 'geneID' here is the Entrez gene identifier. The 'geneSymb'
        // is the "common" gene-name.
        
        public Entry(String line) { 
            String[] array = line.split("\t");
            if(array.length != 6) { 
                throw new IllegalArgumentException("Line \"" + line + "\" doesn't have six toekns."); 
            }
            
            homologeneID = array[0];
            taxID = array[1];
            geneID = array[2]; geneSymb = array[3];
            protein = array[4]; proteinAcc = array[5];
        }
        
        public int hashCode() { 
            int code = 17;
            code += homologeneID.hashCode(); code *= 37;
            code += taxID.hashCode(); code *= 37;
            code += geneID.hashCode(); code *= 37;
            code += geneSymb.hashCode(); code *= 37;
            code += protein.hashCode(); code *= 37;
            code += proteinAcc.hashCode(); code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Entry)) { return false; }
            Entry e = (Entry)o;
            if(!homologeneID.equals(e.homologeneID)) { return false; }
            if(!taxID.equals(e.taxID)) { return false; }
            if(!geneID.equals(e.geneID)) { return false; }
            if(!geneSymb.equals(e.geneSymb)) { return false; }
            if(!protein.equals(e.protein)) { return false; }
            if(!proteinAcc.equals(e.proteinAcc)) { return false; }
            return true;
        }
        
        public String toString() { 
            return geneID + " \"" + geneSymb + "\" [" + taxID + "] --> " + homologeneID; 
        }
    }
}
