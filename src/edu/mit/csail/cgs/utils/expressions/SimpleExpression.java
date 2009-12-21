/*
 * Created on Apr 13, 2005
 */
package edu.mit.csail.cgs.utils.expressions;

import java.util.Set;
import java.util.HashSet;

/**
 * @author tdanford
 */
public class SimpleExpression implements Expression {

    private String value;
    
    public String getValue() { return value; }
    
    public SimpleExpression(String v) { value = v; }
    public String toString() { return value; }
    public boolean isCompound() { return false; }
    
    public Expression substitute(String token, String newValue) { 
        if(value.equals(token)) { 
            return new SimpleExpression(newValue);
        } else { 
            return this;
        }
    }
    
    public Set<String> findFreeTerms() { 
        HashSet<String> s = new HashSet<String>();
        s.add(value);
        return s;
    }
}
