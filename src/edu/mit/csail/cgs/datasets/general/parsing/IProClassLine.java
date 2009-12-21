/*
 * Created on Nov 9, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.general.parsing;

/*
 * iproclass.tb

This table includes the following IDs (or ACs) delimited by tab:

1.  UniProt_ac
2.  UniProt_id
3.  EntrezGene
4.  RefSeq
5.  GIID
6.  PDB
7.  PFAM
8.  GO
9.  PIRSF
10. IPI
11. UniRef_100
12. UniRef_90
13. UniRef_50
14. UniParc
15. PIR-PSD
16. Taxon ID
17. OMIM
18. UniGene
19. Ensemble ID
20. PMID

 */

public class IProClassLine {

    public String[] array;
    public String uniprot_ac, uniprot_id, entrez, refseq, giid, pdb, pfam;
    public String go, pirsf, ipi, uniref_100, uniref_90, uniref_50, uniparc;
    public String pir_psd, taxon, omim, unigene, ensemble, pmid;
    
    public IProClassLine(String line) { 
        array = line.split("\\t");
        uniprot_ac = array[0]; 
        uniprot_id = array[1]; 
        entrez = array[2]; 
        refseq = array[3]; 
        giid = array[4];
        pdb = array[5]; 
        pfam = array[6];
        go = array[7]; 
        pirsf = array[8]; 
        ipi = array[9];
        uniref_100 = array[10]; 
        uniref_90 = array[11]; 
        uniref_50 = array[12]; 
        uniparc = array[13];
        pir_psd = array[14]; 
        taxon = array[15];
        omim = array[16]; 
        unigene = array[17]; 
        ensemble = array[18]; 
        pmid = array[19];
    }
    
    public int hashCode() { 
        int code = 17;
        for(int i = 0; i < array.length; i++) { 
            code += array[i].hashCode(); code *= 37;
        }
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof IProClassLine)) { return false; }
        IProClassLine line = (IProClassLine)o;
        for(int i =0; i < array.length; i++) { 
            if(!array[i].equals(line.array[i])) { return false; }
        }
        return true;
    }

}
