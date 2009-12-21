/*
 * Created on May 12, 2006
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class OBOParser {
    
    public static void main(String[] args) { 
        try {
            OBOParser parser = new OBOParser(new File(args[0]));
            for(OboObject obj : parser.terms) { 
                Map<String,LinkedList<String>> values = obj.getKeyValueMap();
                if(values.containsKey("id")) { 
                    String id = values.get("id").getFirst();
                    System.out.println(id);
                }
            }
            System.out.println(String.format("%d terms", parser.terms.size()));
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }
       
    private OboObject header;
    private Vector<OboObject> terms;

    public OBOParser(File f) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(f)); 
        OboObject h = new OboObject(br);
        if(!h.isInitialized()) { throw new IllegalArgumentException(f.getName()); }
        header = h;
        terms = new Vector<OboObject>();
        do { 
            h = new OboObject(br);
            if(h.isInitialized() && h.getTermTag().equals("Term")) {
                terms.add(h);
            }
        } while(h.isInitialized());
        br.close();
    }

    public int size() { return terms.size(); }
    public OboObject getOboObject(int i) { return terms.get(i); }
    public Collection<String> getHeaderValues(String key) { return header.getValues(key); }
    
    public ProvisionalGOTerm[] createGOTerms() {
        ProvisionalGOTerm[] array = new ProvisionalGOTerm[terms.size()];
        
        Map<String,Integer> indexMap = new HashMap<String,Integer>();
        
        int i = 0;
        for(OboObject oo : terms) { 
            array[i] = new ProvisionalGOTerm(oo);
            indexMap.put(array[i].getID(), i);
            i += 1;
        }
        
        for(i = 0; i < array.length; i++) { 
            OboObject base = array[i].getObject();
            Collection<String> isaValues = base.getValues("is_a");
            
            //System.out.println(array[i].getID() + " ==> ");
            //for(String isa : isaValues) { System.out.println("\t" + isa); }
            
            for(String isa : isaValues) { 
                String[] isa_array = isa.split("!");
                String isaID = isa_array[0].trim();
                if(!indexMap.containsKey(isaID)) { 
                    throw new IllegalArgumentException(isaID); 
                }
                int pgIndex = indexMap.get(isaID);
                array[i].addParent(array[pgIndex]);
                //System.out.println("\t* " + array[pgIndex].getID());
            }
        }
        
        Arrays.sort(array);
        /*
        for(int j = 0; j < array.length; j++) { 
            System.out.println(j + ": (" + array[j].insertPriority + ") : " + array[j].toString());
        }
        */
        
        return array;
    }
    
}
