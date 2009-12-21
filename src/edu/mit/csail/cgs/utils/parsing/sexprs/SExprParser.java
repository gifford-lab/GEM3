/*
 * Author: tdanford
 * Date: Jun 12, 2008
 */
package edu.mit.csail.cgs.utils.parsing.sexprs;

import java.util.*;

public class SExprParser {
	
	public static void main(String[] args) { 
		/** 
		 * Testing...
		 */
		
		String str = "xy xyz (and e (or a    b) \n(not c) d) f";
		SExprParser parser = new SExprParser();
		LinkedList<SExpr> exprs = parser.parseSExprs(str);
		
		int i = 0;
		for(SExpr e : exprs) { 
			System.out.println(i + ": " + e);
			i++;
		}
	}
	
	private SExprParser childParser;
	
	public SExprParser() { 
		childParser = null;
	}
	
	public LinkedList<SExpr> parseSExprs(String str) { 
		StringLexer lexer = new StringLexer();
		return parseSExprs(lexer.lexString(str));
	}
	
	public LinkedList<SExpr> parseSExprs(LinkedList<Token> tokens) { 
		LinkedList<SExpr> exprs = new LinkedList<SExpr>();
		while(!tokens.isEmpty()) { 
			exprs.addLast(parseHeadSExpr(tokens));
		}
		return exprs;
	}
	
	public SExpr parseHeadSExpr(LinkedList<Token> tokens) {
		if(tokens.isEmpty()) { throw new IllegalArgumentException(); }
		Token head = tokens.removeFirst();
		
		if(head.isSExprStart()) {
			LinkedList<SExpr> subExprs = new LinkedList<SExpr>();
			if(childParser == null) { childParser = new SExprParser(); }
			
			while(!tokens.isEmpty() && !tokens.getFirst().isSExprEnd()) { 
				subExprs.addLast(childParser.parseHeadSExpr(tokens));
			}
			
			if(tokens.isEmpty()) { 
				throw new IllegalArgumentException("Missing s-expr end token.");
			} else { 
				tokens.removeFirst();
			}
			
			return new CompoundExpr(subExprs);
			
		} else if (head.isSExprEnd()) { 
			throw new IllegalArgumentException("Unexpected s-expr end token.");
		} else if (head.isWhitespace()) { 
			throw new IllegalArgumentException();
		} else { 
			return new SymbolExpr(head.value);
		}
	}
}

class StringLexer { 
	
	public LinkedList<Token> lexString(String str) { 
		int pos = 0;
		int start = -1;
		LinkedList<Token> tokens = new LinkedList<Token>();
		while(pos < str.length()) {
			char val = str.charAt(pos);
			if(Character.isWhitespace(val)) { 
				while(pos < str.length() && Character.isWhitespace(str.charAt(pos))) { 
					pos++;
				}
				//tokens.add(new Token(" "));
			} else if (isSExprStart(val)) {
				tokens.add(new Token(String.valueOf(str.charAt(pos++))));
				
			} else if (isSExprEnd(val)) { 
				tokens.add(new Token(String.valueOf(str.charAt(pos++))));
				
			} else if (isQuoteStart(val)) { 
				start = pos;
				while(pos < str.length() && 
						!isQuoteEnd(str.charAt(pos))) { 
					pos++;
				}
				
				String tokString = str.substring(start, pos);
				Token t = new Token(tokString);
				t.quoted = true;
				tokens.add(t);
				
			} else {
				start = pos;
				while(pos < str.length() && !Character.isWhitespace(str.charAt(pos)) && 
						!isSExprStart(str.charAt(pos)) && !isSExprEnd(str.charAt(pos))) { 
					pos++;
				}
				
				String tokString = str.substring(start, pos);
				Token t = new Token(tokString);
				tokens.add(t);
			}
		}
		return tokens;
	}
	
	public boolean isEscape(char c) { 
		return c == '\\';
	}
	
	public boolean isQuoteStart(char c) { 
		return c == '"';
	}
	
	public boolean isQuoteEnd(char c) { 
		return c == '"';
	}
	
	public boolean isSExprStart(char c) { 
		return c == '(';
	}
	
	public boolean isSExprEnd(char c) { 
		return c == ')';
	}
}

class Token { 
	public boolean quoted;
	public String value;
	
	public Token(String v) { value = v; quoted = false; }
	
	public boolean isWhitespace() { 
		for(int i = 0; i < value.length(); i++) { 
			if(!Character.isWhitespace(value.charAt(i))) { 
				return false;
			}
		}
		return true;
	}

	public boolean isSymbol() { 
		return !isWhitespace() && !isSExprStart() && !isSExprEnd(); 
	}
	
	public boolean isSExprStart() { return value.equals("("); }
	public boolean isSExprEnd() { return value.equals(")"); }
}

