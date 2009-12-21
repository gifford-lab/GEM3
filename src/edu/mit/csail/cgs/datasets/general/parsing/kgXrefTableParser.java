/*
 * Created on May 19, 2006
 */
package edu.mit.csail.cgs.datasets.general.parsing;

import java.io.*;
import java.util.*;
import java.sql.SQLException;
import java.text.*;

import edu.mit.csail.cgs.datasets.function.DatabaseFunctionLoader;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class kgXrefTableParser {
    

//  kgXref (
//      kgID varchar(40) NOT NULL default '',
//      mRNA varchar(40) default NULL,
//      spID varchar(40) default NULL,
//      spDisplayID varchar(40) default NULL,
//      geneSymbol varchar(40) default NULL,
//      refseq varchar(40) default NULL,
//      protAcc varchar(40) default NULL,
//      description varchar(255) default NULL,
//  }
    
    public static class Entry { 
        private String kgID, mRNA, spID, sDisplayID, geneSymbol, refseq, protAcc, description;
        
        public Entry(String line) { 
            String[] array = line.split("\t");
            if(array.length != 8) { 
                throw new IllegalArgumentException("\"" + line + "\" " + array.length); 
            }
            
            kgID = array[0];
            mRNA = array[1];
            spID = array[2];
            sDisplayID = array[3];
            geneSymbol = array[4];
            refseq = array[5];
            protAcc = array[6];
            description = array[7];
        }
        
        public String get_kgID() { return kgID; }
        public String get_mRNA() { return mRNA; }
        public String get_spID() { return spID; }
        public String get_sDisplayID() { return sDisplayID; }
        public String get_geneSymbol() { return geneSymbol; }
        public String get_refseq() { return refseq; }
        public String get_protAcc() { return protAcc; }
        public String get_description() { return description; }
        
        public int hashCode() { 
            int code = 17;
            code += kgID.hashCode(); code *= 37;
            code += mRNA.hashCode(); code *= 37;
            code += spID.hashCode(); code *= 37;
            code += sDisplayID.hashCode(); code *= 37;
            code += geneSymbol.hashCode(); code *= 37;
            code += refseq.hashCode(); code *= 37;
            code += protAcc.hashCode(); code *= 37;
            code += description.hashCode(); code *= 37;
            return code; 
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Entry)) { return false; }
            Entry e = (Entry)o;
            if(!kgID.equals(e.kgID)) { return false; }
            if(!mRNA.equals(e.mRNA)) { return false; }
            if(!spID.equals(e.spID)) { return false; }
            if(!sDisplayID.equals(e.sDisplayID)) { return false; }
            if(!geneSymbol.equals(e.geneSymbol)) { return false; }
            if(!refseq.equals(e.refseq)) { return false; }
            if(!protAcc.equals(e.protAcc)) { return false; }
            if(!description.equals(e.description)) { return false; }
            return true;
        }
    }
    
    private Vector<Entry> entries;
 
    public kgXrefTableParser(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        entries = new Vector<Entry>();
        int count = 0;
        String line = null;
        while((line = br.readLine()) != null) { 
            if(line.length() > 0) { 
                Entry e = new Entry(line);
                entries.add(e);
                count += 1;
            }
        }
        br.close();
        System.out.println("kgXref Table, # Entries: " + count);
    }
    
    public int size() { return entries.size(); }
    public Entry getEntry(int i) { return entries.get(i); }

    public Map<String,? extends Collection<String>> build_spID2refseq_map() { 
        HashMap<String,LinkedList<String>> map = new HashMap<String,LinkedList<String>>();
        int duplicates = 0;
        for(Entry e : entries) { 
            String spID = e.get_spID();
            String refseq = e.get_refseq();
            if(refseq.length() > 0) { 
                if(!map.containsKey(spID)) { map.put(spID, new LinkedList<String>()); }
                map.get(spID).add(refseq);
            }
        }
        return map;
    }

    public Map<String,? extends Collection<String>> build_kgID2refseq_map() { 
        HashMap<String,LinkedList<String>> map = new HashMap<String,LinkedList<String>>();
        int duplicates = 0;
        for(Entry e : entries) { 
            String kgID = e.get_kgID();
            String refseq = e.get_refseq();
            if(refseq.length() > 0) { 
                if(!map.containsKey(kgID)) { map.put(kgID, new LinkedList<String>()); }
                map.get(kgID).add(refseq);
            }
        }
        return map;
    }

    public Map<String,? extends Collection<String>> build_kgID2spID_map() { 
        HashMap<String,Set<String>> map = new HashMap<String,Set<String>>();
        int duplicates = 0;
        for(Entry e : entries) { 
            String kgID = e.get_kgID();
            //String refseq = e.get_refseq();
            String spID = e.get_spID();
            if(spID.length() > 0) { 
                if(!map.containsKey(kgID)) { map.put(kgID, new HashSet<String>()); }
                map.get(kgID).add(spID);
            }
        }
        return map;
    }

}
