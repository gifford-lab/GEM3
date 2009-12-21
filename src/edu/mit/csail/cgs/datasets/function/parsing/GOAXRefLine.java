/*
 * Created on Feb 12, 2007
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;

import edu.mit.csail.cgs.utils.Pair;

/**
1.  Database from which master entry of this IPI entry has been taken. 
    One of either: 
    SP (UniProt/Swiss-Prot), TR (UniProt/TrEMBL), ENSEMBL (Ensembl),
    REFSEQN (RefSeq NP data set), REFSEQX (Refseq XP data set), TAIR 
    (TAIR Protein data set) or HINV (H-Invitational Database).
2.  UniProt accession number or Ensembl ID or RefSeq ID or 
    TAIR Protein ID or HINV (H-Invitational Database).
3.  International Protein Index identifier (see section 5).
4.  Supplementary UniProt/Swiss-Prot entries associated with 
    this IPI entry.
5.  Supplementary UniProt/TrEMBL entries associated with this 
    IPI entry.
6.  Supplementary Ensembl entries associated with this IPI entry.
7.  Supplementary RefSeq NP entries associated with this IPI entry.
8.  Supplementary RefSeq XP entries associated with this IPI entry.
9.  Supplementary TAIR Protein entries associated with this IPI 
    entry.
10. Supplementary H-Inv Protein entries associated with this IPI 
    entry.
11. Protein identifiers (cross reference to EMBL/Genbank/DDBJ 
    nucleotide databases).
12. List of HGNC number, Genew official gene symbol couples 
    (separated by a semi-colon ';') associated with this IPI entry.
13. List of NCBI Entrez Gene gene number, Entrez Gene Default Gene 
    Symbol couples (separated by a semi-colon ';') associated with 
    this IPI entry.
14. UNIPARC identifier associated with the sequence of this IPI 
    entry.
15. UniGene identifiers associated with this IPI entry. 
 */


public class GOAXRefLine {

    private String database;
    private String dbID;
    private String ipiID;
    private String suppSwissProt, suppTrembl;
    private String suppEnsembl, suppRefseqNP, suppRefseqXP, suppTair, suppHInv;
    private String protIDs;
    private LinkedList<Pair<String,String>> speciesGeneSymbols;
    private LinkedList<Pair<String,String>> entrezGenes;
    private String uniparc;
    private String unigene;
    
    public GOAXRefLine(String line) { 
        String[] array = line.split("\t");
        
        database = array[0]; 
        dbID = array[1];
        ipiID = array[2];
        suppSwissProt = array[3].length() > 0 ? array[3] : null;
        suppTrembl = array[4].length() > 0 ? array[4] : null;
        suppEnsembl = array[5].length() > 0 ? array[5] : null;
        suppRefseqNP = array[6].length() > 0 ? array[6] : null;
        suppRefseqXP = array[7].length() > 0 ? array[7] : null;
        
        //suppTair = array[8];
        suppTair = null;
        
        suppHInv = array[8].length() > 0 ? array[8] : null;
        protIDs = array[9];
        
        String[] symarray = array[10].split(";");
        speciesGeneSymbols = new LinkedList<Pair<String,String>>();
        for(int i = 0; i < symarray.length; i++) { 
        	if(symarray[i].length() > 0) {
        		String[] parray = symarray[i].split(",");
        		speciesGeneSymbols.addLast(new Pair<String,String>(parray[0], parray[1]));
        	}
        }
        
        String[] entarray = array[11].split(";");
        entrezGenes = new LinkedList<Pair<String,String>>();
        for(int i = 0; i < entarray.length; i++) {
        	if(entarray[i].length() > 0) {
        		String[] parray = entarray[i].split(",");
        		entrezGenes.addLast(new Pair<String,String>(parray[0], parray[1]));
        	}
        }
        
        uniparc = array.length >= 13 && array[12].length() > 0 ? array[12] : null;
        unigene = array.length >= 14 && array[13].length() > 0 ? array[13] : null;
        
        /*
        System.out.println("db: \"" + database + "\"");
        System.out.println("dbID: \"" + dbID + "\"");
        System.out.println("ipiID: \"" + ipiID + "\"");
        System.out.println("suppSwissProt: \"" + suppSwissProt + "\"");
        System.out.println("suppTrembl: \"" + suppTrembl + "\"");
        System.out.println("suppEnsembl: \"" + suppEnsembl + "\"");
        System.out.println("suppRefseqNP: \"" + suppRefseqNP + "\"");
        System.out.println("suppRefseqXP: \"" + suppRefseqXP + "\"");
        System.out.println("suppTair: \"" + suppTair + "\"");
        System.out.println("suppHInv: \"" + suppHInv + "\"");
        System.out.println("protIDs: \"" + protIDs + "\"");
        
        System.out.print("gene Symbols:");
        for(String gs : speciesGeneSymbols) { System.out.print(" " + gs); }
        System.out.println();
        
        System.out.print("entrezGenes:");
        for(String e : entrezGenes) { System.out.print(" " + e); }
        System.out.println();
        
        System.out.println("uniparc: \"" + uniparc + "\"");
        System.out.println("unigene: \"" + unigene + "\"");
        System.out.println();
        */
    }
    
    public Vector<String> parseName(String n) {
        Vector<String> names = new Vector<String>();
        String[] array = n.split(";");
        for(int i = 0; i < array.length; i++) { 
            if(array[i].length() > 0) { 
                if(array[i].contains(":")) { 
                    int colon = array[i].indexOf(":");
                    String na = array[i].substring(colon+1, array[i].length());
                    if(na.length()> 0) { names.add(na); }
                } else { 
                    names.add(array[i]);
                }
            }
        }
        return names;
    }
    
    public Collection<String> getAllNames() { 
        TreeSet<String> names = new TreeSet<String>();
        names.add(dbID);
        names.add(ipiID);
        
        for(Pair<String,String> pair : speciesGeneSymbols) {
            names.add(pair.getLast()); 
        }
        
        for(Pair<String,String> pair : entrezGenes) { 
            names.add(pair.getFirst()); 
            names.add(pair.getLast());
        }
        
        if(unigene != null) { names.addAll(parseName(unigene)); }
        if(protIDs != null) { names.addAll(parseName(protIDs)); } 

        return names;
    }
    
    public String getDatabase() { return database; }
    public String getDatabaseID() { return dbID; }
    public Collection<Pair<String,String>> getEntrezIDs() { return entrezGenes; }
    public Collection<Pair<String,String>> getGeneSymbols() { return speciesGeneSymbols; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof GOAXRefLine)) { return false; }
        GOAXRefLine line = (GOAXRefLine)o;
        if(!database.equals(line.database)) { return false; }
        if(!dbID.equals(line.dbID)) { return false; }
        if(!ipiID.equals(line.ipiID)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += database.hashCode(); code *= 37;
        code += dbID.hashCode(); code *= 37;
        code += ipiID.hashCode(); code *= 37;
        return code;
    }
}
