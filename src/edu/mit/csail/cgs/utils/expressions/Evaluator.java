package edu.mit.csail.cgs.utils.expressions;

import java.util.*;
import java.io.*;

public class Evaluator {

	public static void main(String[] args) {
		try {
			eval();
		} catch(IOException ie) {
			ie.printStackTrace(System.err);
		}
	}

	public static void eval() throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.print(">"); System.out.flush();
		String line;
		Evaluator baseEval = new Evaluator();
		while((line = br.readLine()) != null) {
			ExpressionParser ep = new ExpressionParser(line);
			List<Expression> lst = ep.getExprList();
			for(Expression expr : lst) {
				Object obj = baseEval.eval(expr);
				System.out.println(stringify(obj));
			}
			System.out.print(">");
			System.out.flush();
		}
	}

	public static String stringify(Object obj) {
		if(obj == null) { return "null"; }
		if(obj instanceof Boolean) {
			boolean v = ((Boolean)obj).booleanValue();
			if(v) { return "#t"; }
			return "#f";
		}
		return obj.toString();
	}

	private static Map<String,Object> opBindings;

	static {
		opBindings = new HashMap<String,Object>();
		opBindings.put("+", new Addition());
		opBindings.put("-", new Subtraction());
		opBindings.put("*", new Multiplication());
		opBindings.put("/", new Division());
		opBindings.put("=", new Equality());
		opBindings.put("!=", new Inequality());
		opBindings.put("<", new LessThan());
		opBindings.put(">", new GreaterThan());
		opBindings.put("car", new Car());
		opBindings.put("cdr", new Cdr());
		opBindings.put("cons", new Cons());
	}

	private Evaluator parent;
	private Map<String,Object> bindings;

	public Evaluator() {
		parent = null;
		bindings = new HashMap<String,Object>();
	}

	public Evaluator(Evaluator parent, Vector<String> params, Vector<Object> args) {
		this.parent = parent;
		bindings = new HashMap<String,Object>();
		for(int i = 0; i < params.size(); i++) {
			bindings.put(params.get(i), args.get(i));
		}
	}

	public Object lookup(String name) {
		if(bindings.containsKey(name)) { return bindings.get(name); }
		if(parent != null) { 
			return parent.lookup(name); 
		} else {
			if(opBindings.containsKey(name)) {
				return opBindings.get(name);
			}
		}
		return null;
	}

	public Closure createClosure(CompoundExpression args, Expression body) {
		return new Closure(this, args, body);	
	}

	public Object eval(Expression e) {
		if(e instanceof SimpleExpression) {
			SimpleExpression se = (SimpleExpression)e;
			String val = se.getValue();
			if(val.equals("null")) { return null; }
			if(val.equals("#t")) { return Boolean.TRUE; }
			if(val.equals("#f")) { return Boolean.FALSE; }
			if(val.matches("-?[0-9]+")) { return new Integer(val); }
			if(val.matches("-?[0-9]*\\.[0-9]+")) { return new Double(val); }
			return lookup(val);
		} else {
			CompoundExpression ce = (CompoundExpression)e;
			Expression headExpr = ce.getHead();	
			List<Expression> argExprs = ce.getArgExprList();
			String headString = headExpr.toString();

			if(headString.equals("eval")) {
				Expression ee = argExprs.get(0);
				Object res = eval(ee);

				// the following line amounts to an identity function 
				// if the ee is an Expression already, but it converts "lists" 
				// into "Expression"s otherwise, so that something like: 
				// (eval '(+ 1 2)) 
				// works the same as 
				// (eval (cons '+ (cons 1 (cons 2 null))))
				Expression arge = ExpressionParser.parseLeadingExpression(Evaluator.stringify(res));

				return eval(arge);
			}

			if(headString.equals("set!")) {
				SimpleExpression nameExpr = (SimpleExpression)argExprs.get(0);
				Object val = eval(argExprs.get(1));
				bindings.put(nameExpr.getValue(), val);	
				return null;
			}

			if(headString.equals("if")) {
				Boolean test = (Boolean)eval(argExprs.get(0));
				if(test.booleanValue()) {
					return eval(argExprs.get(1));
				} else {
					return eval(argExprs.get(2));
				}
			}
			
			if(headString.equals("quote")) {
				Expression quotedExpr = argExprs.get(0);
				return quotedExpr;
			}

			if(headString.equals("lambda")) {
				CompoundExpression larg = (CompoundExpression)argExprs.get(0);
				Expression barg = argExprs.get(1);
				return createClosure(larg, barg);	
			}

			Object headObj = eval(headExpr);

			Vector<Object> args = new Vector<Object>();
			for(Expression arg : argExprs) {
				args.add(eval(arg));
			}

			if(headObj instanceof Operator) {
				Operator op = (Operator)headObj;
				return op.operate(args);
			}

			Closure c = (Closure)headObj;
			Evaluator cEval = c.getEvaluator();
			Evaluator lEval = new Evaluator(this, c.getParams(), args);	
			return lEval.eval(c.getBody());
		}
	}
}

interface SpecialForm {
	public Object getForm(Evaluator eval, Expression expr);
}

class QuoteForm implements SpecialForm { 
	public QuoteForm() {}
	public Object getForm(Evaluator eval, Expression expr) {
		CompoundExpression ce = (CompoundExpression)expr;
		Expression arg = ce.getArg(0);
		if(arg instanceof SimpleExpression) { return arg.toString(); }
		return new ConsCell((CompoundExpression)arg);
	}
}

class LambdaForm implements SpecialForm {
	public LambdaForm() {}
	public Object getForm(Evaluator eval, Expression expr) {
		CompoundExpression ce = (CompoundExpression)expr;
		SimpleExpression headExpr = (SimpleExpression)ce.getHead();
		CompoundExpression args = (CompoundExpression)ce.getArg(0);
		Expression body = ce.getArg(1);
		return eval.createClosure(args, body);
	}
}

class IfForm implements SpecialForm {
	public IfForm() {}
	public Object getForm(Evaluator eval, Expression expr) {
		CompoundExpression ce = (CompoundExpression)expr;
		SimpleExpression ifTag = (SimpleExpression)ce.getHead();
		Expression test = ce.getArg(0);
		Expression trueExpr = ce.getArg(1);
		Expression falseExpr = ce.getArg(2);
		Object testObj = eval.eval(test);
		Boolean v = (Boolean)testObj;
		if(v.booleanValue()) { return eval.eval(trueExpr);
		} else {
			return eval.eval(falseExpr);	
		}
	}
}

class EvalForm implements SpecialForm { 
    public EvalForm() {}
    public Object getForm(Evaluator eval, Expression expr) { 
        CompoundExpression ce = (CompoundExpression)expr;
        Expression arg = ce.getArg(0);
        Object argValue = eval.eval(arg);
        String argString = argValue.toString();
        Expression parsedExpr = ExpressionParser.parseLeadingExpression(argString);
        return eval.eval(parsedExpr);
    }
}

interface Operator {
	public Object operate(Vector<Object> args);	
}

class Multiplication implements Operator {
	public Multiplication() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Integer(i1 * i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Double(i1 * i2);
		}

		throw new IllegalArgumentException();
	}
}
class Division implements Operator {
	public Division() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Integer(i1 / i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Double(i1 / i2);
		}

		throw new IllegalArgumentException();
	}
}
class Subtraction implements Operator {
	public Subtraction() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Integer(i1 - i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Double(i1 - i2);
		}

		throw new IllegalArgumentException();
	}
}

class GreaterThan implements Operator {
	public GreaterThan() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Boolean(i1 > i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Boolean(i1 > i2);
		}

		throw new IllegalArgumentException();
	}
}

class LessThan implements Operator {
	public LessThan() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Boolean(i1 < i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Boolean(i1 < i2);
		}

		throw new IllegalArgumentException();
	}
}

class Equality implements Operator {
	public Equality() {}
	public Object operate(Vector<Object> args) {
		return new Boolean(Evaluator.stringify(args.get(0)).equals(Evaluator.stringify(args.get(1))));
	}
}

class Inequality implements Operator {
	public Inequality() {}
	public Object operate(Vector<Object> args) {
		return new Boolean(!Evaluator.stringify(args.get(0)).equals(Evaluator.stringify(args.get(1))));
	}
}

class Addition implements Operator {
	public Addition() {}
	public Object operate(Vector<Object> args) {
		Object a1 = args.get(0), a2 = args.get(1);

		if(a1 instanceof Integer && a2 instanceof Integer) {
			int i1 = ((Integer)a1).intValue(), i2 = ((Integer)a2).intValue();
			return new Integer(i1 + i2);
		}

		if(a1 instanceof Double && a2 instanceof Double) {
			double i1 = ((Double)a1).doubleValue(), i2 = ((Double)a2).doubleValue();
			return new Double(i1 + i2);
		}

		throw new IllegalArgumentException();
	}
}

class Cons implements Operator {
	public Cons() {}
	public Object operate(Vector<Object> args) {
		return new ConsCell(args.get(0), args.get(1));
	}
}

class Cdr implements Operator {
	public Cdr() {}
	public Object operate(Vector<Object> args) {
		ConsCell cc = (ConsCell)args.get(0);
		return cc.getCdr();
	}
}

class Car implements Operator {
	public Car() {}
	public Object operate(Vector<Object> args) {
		ConsCell cc = (ConsCell)args.get(0);
		return cc.getCar();
	}
}

class ConsCell {
	private Object car, cdr;

	public ConsCell(Object head, Object tail) {
		car = head;
		cdr = tail;
	}

	public ConsCell(CompoundExpression expr) {
		Expression h = expr.getHead();
		if(h instanceof SimpleExpression) {
			car = ((SimpleExpression)h).getValue();
		} else {
			car = new ConsCell((CompoundExpression)h);
		}
		if(expr.getNumArgs() == 0) {
			cdr = null;
		} else {
			cdr = new ConsCell(new CompoundExpression(expr.getArgExprList()));
		}
	}
	 
	public LinkedList asList() {
		if(cdr == null) {
			LinkedList l = new LinkedList();
			l.addLast(car);
			return l;
		}

		if(!(cdr instanceof ConsCell)) { throw new IllegalArgumentException(); }
		ConsCell cc = (ConsCell)cdr;
		LinkedList ll = cc.asList();
		ll.addFirst(car);
		return ll;
	}

	public boolean isList() {
		if(cdr == null) { return true; }
		if(cdr instanceof ConsCell) { return ((ConsCell)cdr).isList(); }
		return false;
	}

	public Object getCar() { return car; }
	public Object getCdr() { return cdr; }

	public String toString() {
		StringBuilder sb = new StringBuilder();
		if(isList()) {
			LinkedList ll = asList();
			int i = 0;
			sb.append("(");
			for(Object o : ll) {
				sb.append(Evaluator.stringify(o));
				if(i < ll.size()-1) { sb.append(" "); }
				i++;
			}
			sb.append(")");
			return sb.toString();
		}

		sb.append("(");
		sb.append(Evaluator.stringify(car));
		sb.append(" . ");
		sb.append(Evaluator.stringify(cdr));
		sb.append(")");
		return sb.toString();
	}
}

class Closure {
	private Vector<String> params;
	private Evaluator env;
	private Expression body;
	
	public Closure(Evaluator env, CompoundExpression p, Expression b) {
		this.env = env;
		params = new Vector<String>();
		Set<String> paramSet = new HashSet<String>();
		for(Expression e : p.getExprs()) {
			SimpleExpression se = (SimpleExpression)e;
			if(!paramSet.contains(se.getValue())) {
				params.add(se.getValue());
			}
		}
		body = b;
	}

	public Vector<String> getParams() { return params; }
	public Expression getBody() { return body; }
	public Evaluator getEvaluator() { return env; }
}


