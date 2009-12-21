/*
 * Created on Apr 13, 2005
 */
package edu.mit.csail.cgs.utils.expressions;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class ExpressionParser {
    
    public static void main(String[] args) { 
        for(int i = 0; i < args.length; i++) { 
            System.out.println("Arg " + i + ":");
            ExpressionParser ep = new ExpressionParser(args[i]);
            ep.print(System.out);
            System.out.print("\n");
        }
    }

    public static Expression parseLeadingExpression(String str) { 
        ParserState state = new ParserState(str);
        return parseExpression(state);
    }
    
    private static Expression parseExpression(ParserState state) {
        Expression e = null;
        state.eatWhitespace();
		boolean quoted = false;
		if(state.isQuoteChar()) {
			state.advance();
			quoted = true;
		}

        if(state.isCompoundStart()) { 
            e = parseCompoundExpression(state);
        } else { 
            String token = state.getToken();
            e = new SimpleExpression(token);
        }
        state.eatWhitespace();

		if(quoted) {
			LinkedList<Expression> c = new LinkedList<Expression>();
			c.addLast(new SimpleExpression("quote"));
			c.addLast(e);
			e = new CompoundExpression(c);
		}

        return e;
    }
    
    private static CompoundExpression parseCompoundExpression(ParserState state) { 
        if(!state.isCompoundStart()) { 
            throw new IllegalArgumentException();
        }
        state.advance();
        LinkedList<Expression> args = new LinkedList<Expression>();
        Expression header = parseExpression(state);
        while(!state.isCompoundEnd()) { 
            args.addLast(parseExpression(state));
        }
        state.advance();
        CompoundExpression ce = null;
        if(header.toString().equals("lambda")) { 
            CompoundExpression params = (CompoundExpression)args.get(0);
            CompoundExpression body = (CompoundExpression)args.get(1);
            ce = new LambdaExpression(params, body);
        } else { 
            ce = new CompoundExpression(header, args);
        }
        return ce;
    }

    private Vector<Expression> exprList;
    
    public ExpressionParser(String parsable) { 
        ParserState state = new ParserState(parsable);
        exprList = new Vector<Expression>();
        while(!state.isFinished()) { 
            exprList.add(parseExpression(state));
        }
    }
    
    public List<Expression> getExprList() { return new LinkedList<Expression>(exprList); }
    public int size() { return exprList.size(); }
    public Expression getExpr(int i) { return exprList.get(i); }
    
    public void print(PrintStream ps) { 
        int i = 0;
        for(Expression expr : exprList) { 
            ps.print(i + ": ");
            i++;
            ps.println(expr.toString());
        }
    }
}

class ParserState { 
    private char[] array;
    private int pos;
    
    public ParserState(String str) { 
		//System.out.println("STATE: [" + str + "]");
        array = str.toCharArray();
        pos = 0;
    }
    
    public boolean isFinished() { return pos >= array.length; }
    public boolean isCompoundStart() { return array[pos] == '('; }
    public boolean isCompoundEnd() { return array[pos] == ')'; }
	public boolean isQuoteChar() { return array[pos] == '\''; }
    public boolean isCompoundChar() { return isCompoundStart() || isCompoundEnd(); }
    public boolean isWhitespace() { return Character.isWhitespace(array[pos]); }
    
    public void eatWhitespace() { 
        while(!isFinished() && isWhitespace()) { 
            advance();
        }
    }
    
    public char getCurrentChar() { return array[pos]; }
    public void advance() { pos++; }
    public void reset() { pos = 0; }
    
    public String getToken() { 
        StringBuilder sb = new StringBuilder();
        while(!isFinished() && !isWhitespace() && !isCompoundChar()) { 
            sb.append(getCurrentChar());
            advance();
        }
        if(!isFinished() && isCompoundStart()) { 
            throw new IllegalArgumentException();
        }
		//System.out.println("TOKEN: " + sb.toString());
        return sb.toString();
    }
}
