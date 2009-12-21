/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * @author tdanford
 */
public interface RecurrenceEquation {
    
    public double getValue(Arguments a); 
    
    public static class Arguments { 
        
        private int[] values;
        
        public Arguments(Collection<Integer> a) { 
            values = new int[a.size()];
            int i = 0; 
            for(int v : a) { 
                values[i++] = v;
            }
        }
        
        public Arguments(int[] a) { 
            values = (int[])a.clone();
        }
        
        public int size() { return values.length; }
        public int getArgument(int i) { return values[i]; }
        
        public int hashCode() { 
            int code = 17;
            for(int i = 0; i < values.length; i++) { 
                code += values[i]; 
            }
            code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Arguments)) { return false; }
            Arguments a = (Arguments)o;
            if(values.length != a.values.length) { return false; }
            for(int i = 0; i < values.length; i++) { 
                if(values[i] != a.values[i]) { return false; }
            }
            return true;
        }
        
        public String toString() { 
            StringBuilder sb = new StringBuilder();
            sb.append("{ ");
            for(int i = 0; i < values.length; i++) { 
                sb.append(values[i] + " ");
            }
            sb.append("}");
            return sb.toString();
        }
    }
}
