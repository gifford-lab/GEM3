/*
 * Created on Nov 9, 2006
 */
package edu.mit.csail.cgs.datasets.general.parsing;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class IProClassParser implements Iterator<IProClassLine> {
    
    public static Map<String,Set<String>> createMap(File f, int i1, int i2) { 
        HashMap<String,Set<String>> map = new HashMap<String,Set<String>>();
        try {
            IProClassParser parser = new IProClassParser(f);
            while(parser.hasNext()) { 
                IProClassLine line = parser.next();
                String v1 = line.array[i1], v2 = line.array[i2];
                if(!map.containsKey(v1)) { map.put(v1, new HashSet<String>()); }
                map.get(v1).add(v2);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return map;
    }
    
    private BufferedReader br;
    private IProClassLine nextLine;
    
    public IProClassParser(File f) throws IOException { 
        br = new BufferedReader(new FileReader(f));
        String line = br.readLine();
        if(line != null) { nextLine = new IProClassLine(line); }
    }

    public boolean hasNext() {
        return nextLine != null;
    }

    public IProClassLine next() {
        IProClassLine ret = nextLine;
        String line;
        
        try {
            line = br.readLine();
        } catch (IOException e1) {
            line = null;
            e1.printStackTrace();
        }
        
        if(line != null) { 
            nextLine = new IProClassLine(line);
        } else {
            nextLine = null;
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }    
}
