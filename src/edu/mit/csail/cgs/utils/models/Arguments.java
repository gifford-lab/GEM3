/*
 * Author: tdanford
 * Date: Nov 6, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.lang.reflect.*;
import java.util.*;

/**
 * Model-based command-line argument parsing!  
 * 
 * The class takes a Model subclass in the constructor.  When the parse() method is called,
 * it attempts to parse the given arguments and put them in an instance of the Model, which 
 * it then returns. 
 * 
 * If it sees a String field "key" in the Model, it looks to fill it from an instance of 
 * "--key value" on the command line.  
 * 
 * If it sees a Boolean field "flag" in the Model, it looks to fill it with 'true' if 
 * the value-less flag "--flag" appears on the command line, or 'false' otherwise.
 * 
 * Any remaining fields that are either key-less, or argument listed after the separator 
 * "--", it tries to stuff into a field called 'args' in the Model, if one exists and has the 
 * type String[].  (The name of this field can be modified by passing an additional argument
 * to the constructor.)  
 * 
 * @author tdanford
 *
 * @param <ArgModel>
 */
public class Arguments<ArgModel extends ArgumentModel> {
	
	public static void main(String[] args) { 
		Arguments<Test> a = new Arguments<Test>(Test.class);
		Test t = a.parse(args);
		if(!t.checkArgs()) { 
			System.err.println(t.getArgErrors()); 
		} else { 
			System.out.println(t.toString());
		}
	}
	
	public static class Test extends ArgumentModel {
		
		public static String[] required = { "bar" };
		
		public String[] args;
		public String foo, bar;
		public Integer x, y;
		public Boolean flag;
		
		public Integer y() { return 2*x; }
	}

	private String remainingArgsName;
	private Class<ArgModel> modelClass;
	private ModelFieldAnalysis<ArgModel> analysis;
	
	public Arguments(Class<ArgModel> c) { 
		analysis = new ModelFieldAnalysis<ArgModel>(c);
		modelClass = c;
		remainingArgsName = "args";
	}
	
	public Arguments(Class<ArgModel> c, String rem) { 
		analysis = new ModelFieldAnalysis<ArgModel>(c);
		modelClass = c;
		remainingArgsName = rem;		
	}

	public ArgModel parse(String[] args) {
		
		try {
			ArgModel model = modelClass.newInstance();

			/*
			 * We need to make sure that the Boolean fields get a "false" value by 
			 * default.  Everything else will (presumably) get a null, unless there is an 
			 * argument-class-specific default constructor that sets up alternate default 
			 * values.  
			 */
			for(Field f : analysis.findTypedFields(Boolean.class)) { 
				f.set(model, Boolean.FALSE);
			}

			ArrayList<String> remlist = new ArrayList<String>();

			argloop: for(int i = 0; i < args.length; i++) { 

				if(args[i].startsWith("--")) {
					
					if(args[i].equals("--")) {
						// We've seen a divider, so we need to parse the rest of the list
						// and put *all* the values into the 'remainder' list.
						
						for(int j = i + 1; j < args.length; j++) { 
							remlist.add(args[j]);
						}
						
						break argloop;  // quits the for-loop altogether.
					} else {
						String key = args[i].substring(2, args[i].length());
						boolean isFlag = (i == args.length-1) || args[i+1].startsWith("--"); 
						
						if(isFlag) { 
							// This parses the "--flag" without a corresponding "value" entry.
							
							Field f = analysis.findField(key);
							if(Model.isSubclass(f.getType(), Boolean.class)) { 
								f.set(model, (Boolean)true);
							}
							
						} else { 
							// We're parsing a "--key value" pair in this block.
							
							String value = args[++i];  
							// the ++i in that line, makes sure that we don't double-parse the 
							// args[i+1] element (the 'value'). 
							
							Field f = analysis.findField(key);
							if(f != null) { 
								Class type = f.getType();
								try { 
									if(Model.isSubclass(type, String.class)) { 
										f.set(model, value);

									} else if (Model.isSubclass(type, Integer.class)) {
										f.set(model, Integer.parseInt(value));

									} else if (Model.isSubclass(type, Double.class)) { 
										f.set(model, Double.parseDouble(value));
									}
								} catch(NumberFormatException e) { 
									throw new IllegalArgumentException(key);
								}
							}
						}
					}
					
				} else { 
					// This is an entry (a) before a "--" divider, but *not* preceded by a 
					// --key element.  So we put it in the "remainder" list, to be handled 
					// at the end.
					
					remlist.add(args[i]);
				}
			}

			/*
			 * Any non-keyed, or trailing, arguments were put in the 'remlist' variable.
			 * At this point, we package them all up in an array, and put them in the 
			 * contents of the remainingArgsName field (assuming that field is a String[] type). 
			 */
			Field remf = analysis.findField(remainingArgsName);
			if(remf != null && Model.isSubclass(remf.getType(), String[].class)) {
				String[] remarray = remlist.toArray(new String[remlist.size()]);
				remf.set(model, remarray);
			}			

			/*
			 * This block looks for any non-static methods, that
			 * (a) have the same name as a field
			 * (b) have a return-type that is a subclass of the corresponding field type
			 * (c) take no arguments. 
			 * 
			 * For each such method that it finds, it calls the method and sets the value 
			 * of the field to be the returned value from the method invocation.
			 *  
			 * This lets us define argument-specific transformations -- for instance, 
			 * loading a genome object from a genome string, etc.  
			 */
			Class modelClass = analysis.getModelClass();
			Method[] methods = modelClass.getMethods();
			for(int i = 0; i < methods.length; i++) { 
				if((methods[i].getModifiers() & Modifier.STATIC) == 0) { 
					if(methods[i].getParameterTypes().length == 0) { 
						String methodName = methods[i].getName();
						Field field = analysis.findField(methodName);
						if(field != null && Model.isSubclass(methods[i].getReturnType(), field.getType())) { 
							Object value = methods[i].invoke(model, null);
							if(value != null) { 
								field.set(model, value);
							}
						}
					}
				}
			}

			return model;
			
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		}

		return null;
	}
}

