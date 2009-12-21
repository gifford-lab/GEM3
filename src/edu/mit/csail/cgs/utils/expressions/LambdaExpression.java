/*
 * Created on May 15, 2005
 */
package edu.mit.csail.cgs.utils.expressions;

import java.util.*;

/**
 * @author tdanford
 */
public class LambdaExpression extends CompoundExpression {
    
    private CompoundExpression params, body;
    private Vector<String> paramNames;
    private Set<String> paramSet;
    
    public LambdaExpression(CompoundExpression params, CompoundExpression body) { 
        super(new SimpleExpression("lambda"));
        this.params = params;
        this.body = body;
        exprs.add(params);
        exprs.add(body);
        paramSet = new HashSet<String>();
        paramNames = new Vector<String>();
        for(Expression e : params.exprs) { 
            if(!(e instanceof SimpleExpression)) { 
                throw new IllegalArgumentException();
            }
            String token = ((SimpleExpression)e).getValue();
            if(paramSet.contains(token)) { throw new IllegalArgumentException(); }
            paramSet.add(token);
            paramNames.add(token);
        }
    }
    
    public Expression substitute(String token, String value) {
        CompoundExpression newBody = body;
        if(!(paramSet.contains(token))) { 
            newBody = (CompoundExpression)newBody.substitute(token, value);
        }
        return new LambdaExpression(params, newBody);
    }
    
    public Expression call(Vector<Expression> args) { 
        if(args.size() != paramNames.size()) { throw new IllegalArgumentException(); }
        Expression e = body;
        for(int i = 0; i < args.size(); i++) {
            e = e.substitute(paramNames.get(i), args.get(i).toString());
        }
        return e;
    }
    
    public Set<String> findFreeTerms() { 
        Set<String> s = super.findFreeTerms();
        HashSet<String> ns = new HashSet<String>();
        for(String term : s) { 
            if(!paramSet.contains(term)) { ns.add(term); }
        }
        return ns;
    }
}
