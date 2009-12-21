/*
 * Created on Apr 13, 2005
 */
package edu.mit.csail.cgs.utils.expressions;

import java.util.*;

/**
 * @author tdanford
 */
public class CompoundExpression implements Expression {
    protected Vector<Expression> exprs;
    
    public CompoundExpression(Expression head) { 
        exprs = new Vector<Expression>();
        exprs.add(head);
    }
    
    public CompoundExpression(Expression head, Collection<Expression> args) { 
        exprs = new Vector<Expression>();
        exprs.add(head);
        exprs.addAll(args);
    }
    
    public CompoundExpression(Collection<Expression> args) { 
        exprs = new Vector<Expression>(args);
    }
    
    public Expression substitute(String token, String val) { 
        Vector<Expression> expr = new Vector<Expression>();
        for(Expression e : exprs) { 
            expr.add(e.substitute(token, val));
        }
        return new CompoundExpression(expr);
    }
    
    public Set<String> findFreeTerms() { 
        HashSet<String> s = new HashSet<String>();
        for(Expression e : exprs) { 
            s.addAll(e.findFreeTerms());
        }
        return s;
    }
    
    public int size() { return exprs.size(); }
    public Expression getHead() { return exprs.get(0); }
    public int getNumArgs() { return exprs.size()-1; }

	public Expression getTail() {
		if(exprs.size() <= 1) { return null; }
		return new CompoundExpression(getArgExprList());
	}
	

	public List<Expression> getExprs() {
		return new LinkedList<Expression>(exprs);
	}
	
    public List<Expression> getArgExprList() { 
        LinkedList<Expression> lst = new LinkedList<Expression>();
        for(int i = 1; i < exprs.size(); i++) { 
            lst.addLast(exprs.get(i));
        }
        return lst;
    }
    public Expression getExpr(int i) { return exprs.get(i); }
    public Expression getArg(int i) { return exprs.get(i+1); }
    public boolean isCompound() { return true; }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        for(int i = 0; i < exprs.size(); i++) { 
            Expression e = exprs.get(i);
            sb.append(e.toString());
            if(i < exprs.size()-1) { sb.append(" "); }
        }
        sb.append(")");
        return sb.toString();
    }
}

