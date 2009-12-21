/*
 * Created on Apr 13, 2005
 */
package edu.mit.csail.cgs.utils.expressions;

import java.util.Set;

/**
 * @author tdanford
 */
public interface Expression {
    public boolean isCompound();
    public Expression substitute(String token, String substitution);
    public Set<String> findFreeTerms();
}
