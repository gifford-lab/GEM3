/*
 * Created on Dec 19, 2005
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.Collection;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * @author tdanford
 */
public class MIPSCategORF { 
    private String orf;
    private SortedSet<String> categs;
    private int code;
    
    public MIPSCategORF(String inputLine) {
        orf = null;
        categs = new TreeSet<String>();
        code = 17;
        addLine(inputLine);
    }        
    
    private void recalcHashCode() { 
        code = 17;
        code += orf.hashCode(); code *= 37;
        for(String v : categs) { 
            code += v.hashCode(); code *= 37;
        }
    }
    
    public void addLine(String line) { 
        String[] entries = line.split("\\|");
        if(orf == null) { 
            orf = entries[0].toUpperCase();
        } else { 
            if(!orf.equals(entries[0].toUpperCase())) { 
                throw new IllegalArgumentException(orf + " --> [" + line + "]");
            }
        }

        categs.add(entries[1]);
        recalcHashCode();
    }
    
    public String getORF() { return orf; }
    public Collection<String> getCategories() { return categs; }

    public boolean hasCategory(String name) { 
        return categs.contains(name);   
    }

    public boolean hasAnyCategory(Set<String> names) { 
        for(String n : names) { 
            if(categs.contains(n)) { return true; }
        }
        return false;
    }
    
    public int hashCode() { return code; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof MIPSCategORF)) { return false; }
        MIPSCategORF other = (MIPSCategORF)o;
        if(!orf.equals(other.orf)) { return false; }
        if(!categs.equals(other.categs)) { return false; }
        return true;
    }
}