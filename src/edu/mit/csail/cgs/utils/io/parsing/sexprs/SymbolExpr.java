/*
 * Author: tdanford
 * Date: Jun 12, 2008
 */
package edu.mit.csail.cgs.utils.io.parsing.sexprs;

public class SymbolExpr implements SExpr {
	
	private String value;

	public SymbolExpr(String v) { 
		value = v;
	}
	
	public int hashCode() { return value.hashCode(); }
	
	public boolean equals(Object o) { 
		if(!(o instanceof SymbolExpr)) { return false; }
		SymbolExpr tes = (SymbolExpr)o;
		return tes.value.equals(value);
	}

	public boolean isCompound() {
		return false;
	}

	public int length() {
		return 1;
	}
	
	public String getValue() { return value; }
	public String toString() { return value; }

	public SExpr subExpr(int i) {
		throw new IllegalArgumentException(String.valueOf(i));
	} 
	
}