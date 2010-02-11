/*
 * Author: tdanford
 * Date: Jun 12, 2008
 */
package edu.mit.csail.cgs.utils.io.parsing.sexprs;

import java.util.ArrayList;
import java.util.Collection;

public class CompoundExpr implements SExpr {
	
	private ArrayList<SExpr> args;
	
	public CompoundExpr(Collection<SExpr> as) { 
		args = new ArrayList<SExpr>(as);
	}
	
	public CompoundExpr(SExpr... exprs) {
		args = new ArrayList<SExpr>();
		for(int i = 0; i < exprs.length; i++) { 
			args.add(exprs[i]);
		}
	}

	public SExpr car() { return args.get(0); }
	public CompoundExpr cdr() { return new CompoundExpr(args.subList(1, args.size())); }
	
	public boolean isCompound() {
		return true;
	}
	
	public void addSubExpr(SExpr e) { 
		args.add(e);
	}

	public int length() {
		return args.size();
	}

	public SExpr subExpr(int i) {
		return args.get(i);
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("(");
		for(SExpr e : args) { 
			if(sb.length() > 1) { sb.append(" "); }
			sb.append(e.toString());
		}
		sb.append(")");
		return sb.toString();
	}
	
	public int hashCode() { 
		int code = 17;
		for(SExpr e : args) { 
			code += e.hashCode(); code *= 37;
		}
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof CompoundExpr)) { 
			return false;
		}
		
		CompoundExpr c = (CompoundExpr)o;
		if(args.size() != c.args.size()) { return false; }
		for(int i = 0; i < args.size(); i++) { 
			if(!args.get(i).equals(c.args.get(i))) { 
				return false;
			}
		}
		return true;
	}
}