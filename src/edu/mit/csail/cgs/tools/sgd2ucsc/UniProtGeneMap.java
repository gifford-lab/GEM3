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

public class UniProtGeneMap {
    
    public static void main(String[] args) { 
        File input = args.length > 0 ? new File(args[0]) : 
            new File("C:\\Documents and Settings\\tdanford\\Desktop\\dbxref.tab");
        UniProtGeneMap map = new UniProtGeneMap(input);
    }

    public static String idType = "UniProt/Swiss-Prot ID";
    
    private Map<String,String> gene2uniprot;

    public UniProtGeneMap(File f) { 
        DBXRefParser xrefParser = new DBXRefParser(f);
        gene2uniprot = new HashMap<String,String>();
        int duplicates = 0;
        for(DBXRefParser.XRefLine line : xrefParser.lines) {
            if(line.type.equals(idType)) { 
                if(gene2uniprot.containsKey(line.scFeature)) { duplicates += 1; }
                gene2uniprot.put(line.scFeature, line.id);
            }
        }
        
        System.out.println("Gene->Uniprot: " + gene2uniprot.size() + " entries.");
        System.out.println("\t#Duplicates: " + duplicates);
    }
    
    public boolean containsGeneName(String gn) { return gene2uniprot.containsKey(gn); }
    public String getUniprot(String geneName) { return gene2uniprot.get(geneName); }
    
    public Set<String> getGeneNames(String uniprot) { 
        HashSet<String> gns = new HashSet<String>();
        for(String k : gene2uniprot.keySet()) { 
            if(gene2uniprot.get(k).equals(uniprot)) { 
                gns.add(k);
            }
        }
        return gns;
    }
    
}
