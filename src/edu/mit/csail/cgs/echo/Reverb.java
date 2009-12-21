/*
 * Created on Feb 16, 2007
 */
package edu.mit.csail.cgs.echo;

import java.util.*;

import edu.mit.csail.cgs.echo.components.ConstantGenerator;
import edu.mit.csail.cgs.echo.gui.EchoGUIReverb;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Generator;

/**
 * @author Timothy Danford
 * 
 * Reverb is one of the core classes of Echo -- It is a template, 
 * which handles the constants that work as parameters for the EchoProcessors, 
 * and the other Reverb objects that are the inputs to the processor.
 * 
 * At runtime, it is compiled into an EchoProcessor object, and handed
 * off to the scheduler.
 */
public class Reverb implements EchoInputInterface {

	private Map<String,EchoConstant> params;
	private Map<String,Reverb> inputs;

	private Map<String,Integer> inputIndices;
	private Map<String,Integer> paramIndices;

	private SelfDescribingVerb verb;
	private EchoProcessor proc;

	private EchoGUIReverb guiPeer;

	public Reverb(SelfDescribingVerb p) { 
		params = new HashMap<String,EchoConstant>();
		inputs = new HashMap<String,Reverb>();
		inputIndices = new HashMap<String,Integer>();
		paramIndices = new HashMap<String,Integer>();

		guiPeer = null;
		verb = p;
		proc = null;

		String[] inputNames = verb.getInputNames();
		for(int i = 0; inputNames != null && i < inputNames.length; i++) { 
			inputIndices.put(inputNames[i], i);
		}

		String[] paramNames = verb.getParameterNames();
		for(int i = 0; paramNames != null && i < paramNames.length; i++) { 
			paramIndices.put(paramNames[i], i);
		}
	}

	public SelfDescribingVerb getVerb() { return verb; }
	public EchoGUIReverb getPeer() { return guiPeer; }
	public void setPeer(EchoGUIReverb p) { guiPeer = p; }

	public String[] findMatchingInputNames(Reverb r) {
		if(verb.getOutputClass() == null) { return null; }

		String[] inputNames = r.getVerb().getInputNames();
		if(inputNames == null) { return null; }

		EchoType[] inputClasses = r.getVerb().getInputClasses();
		Vector<String> admNames = new Vector<String>();

		for(int i = 0; i < inputNames.length; i++) {
            if(r.isLegalInput(inputNames[i], this)) { 
				admNames.add(inputNames[i]);
			}
		}

		String[] array = admNames.toArray(new String[admNames.size()]);
		return array;
	}

	public EchoConstant getParam(String n) { 
		return params.containsKey(n) ? params.get(n) : null;
	}

	public Reverb getInput(String n) { 
		return inputs.containsKey(n) ? inputs.get(n) : null;
	}

	public void setInput(String n, Reverb r) throws EchoException {
		if(!isLegalInput(n, r)) { 
			System.err.println(verb.toString() + ": " + n + " --> " + inputIndices.keySet());
			String outputName = r.getVerb().getOutputClass().getName();
			String inputName = 
				inputIndices.containsKey(n) ? verb.getInputClasses()[inputIndices.get(n)].getName() : "unknown";
				throw new EchoException("Input " + outputName + " doesn't match \"" + n + "\" (" + inputName + ")");
		}

		inputs.put(n, r);

		if(verb instanceof DependentSelfDescribingVerb) { 
			DependentSelfDescribingVerb v = (DependentSelfDescribingVerb)verb;
			v.setInput(n, r.getVerb().getOutputClass()); 
		}
	}

	public void clearInput(String n) throws EchoException { 
		if(!inputIndices.containsKey(n)) { throw new EchoException("No such input \"" + n + "\""); }
		if(inputs.containsKey(n)) { inputs.remove(n); }
		if(verb instanceof DependentSelfDescribingVerb) { 
			DependentSelfDescribingVerb v = (DependentSelfDescribingVerb)verb;
			v.clearInput(n);
		}
	}

	public void clearInput(Reverb r) { 
		Set<String> names = new HashSet<String>();
		for(String inputName : inputs.keySet()) { 
			if(inputs.get(inputName).equals(r)) { 
				names.add(inputName);
			}
		}

		for(String n : names) { 
			inputs.remove(n);
		}
	}

	public boolean isLegalInput(String n, Reverb r) { 
		if(inputIndices.containsKey(n)) { 
			int idx = inputIndices.get(n);
			EchoType[] inputClasses = getVerb().getInputClasses();
			EchoType inputClass = inputClasses[idx];
			EchoType outputClass = r.getEvalType();
			
			if(!outputClass.isSubType(inputClass) && 
					!isExpandingSubType(outputClass, inputClass)) { 
				return false;
			}
		} else { 
			return false;
		}
		return true;
	}
	
	private static boolean isExpandingSubType(EchoType output, EchoType input) { 
		if(input instanceof SequenceType) { 
			SequenceType st = (SequenceType)input;
			return output.isSubType(st.getInnerType());
		}
		return false;
	}

	public static boolean isSubclass(Class subc, Class superc) {
		return superc.isAssignableFrom(subc);
	}

	public EchoType getEvalType() { return getVerb().getOutputClass(); }

	public String[] getParamNames() { return getVerb().getParameterNames(); }

	public boolean isLegalParam(String name, EchoConstant sdconst) {
		if(!paramIndices.containsKey(name)) { return false; }
		int index = paramIndices.get(name);
		EchoType constClass = sdconst.getConstantClass();
		EchoType verbParamClass = verb.getParameterClasses()[index];
		if(!constClass.isSubType(verbParamClass)) { return false; }
		return true;
	}


	public void setParam(String name, EchoConstant sdconst) throws EchoException {
		if(!isLegalParam(name, sdconst)) { 
			throw new EchoException("Parameter \"" + name + "\" doesn't match class " + 
					sdconst.getConstantClass().getName());    
		}
		params.put(name, sdconst);

		if(verb instanceof DependentSelfDescribingVerb) { 
			DependentSelfDescribingVerb v = (DependentSelfDescribingVerb)verb;
			v.setParameter(name, sdconst.getConstantClass());
		}
	}

	public void clearParam(String name) throws EchoException { 
		if(!paramIndices.containsKey(name)) { throw new EchoException("No such parameter \"" + name + "\""); }
		if(params.containsKey(name)) { params.remove(name); }

		if(verb instanceof DependentSelfDescribingVerb) { 
			DependentSelfDescribingVerb v = (DependentSelfDescribingVerb)verb;
			v.clearParameter(name);
		}
	}

	public void clearParam(EchoConstant c)  { 
		Set<String> names = new HashSet<String>();
		for(String pn : params.keySet()) { 
			if(params.get(pn).equals(c)) { 
				names.add(pn);
			}
		}

		for(String n : names) { params.remove(n); }
	}

	public Set<String> findMissingInputs() { 
		HashSet<String> names = new HashSet<String>();
		String[] inputNames = verb.getInputNames();
		for(int i = 0; inputNames != null && i < inputNames.length; i++) { 
			if(!inputs.containsKey(inputNames[i])) { names.add(inputNames[i]); }
		}
		return names;
	}

	public Set<String> findMissingParams() { 
		HashSet<String> names = new HashSet<String>();
		String[] paramNames = verb.getParameterNames();
		for(int i = 0; paramNames != null && i < paramNames.length; i++) { 
			if(!params.containsKey(paramNames[i])) { names.add(paramNames[i]); }
		}
		return names;    	
	}

	public Vector<String> findMatchingInputs(Reverb r, boolean listAll) { 
		TreeSet<String> names = new TreeSet<String>();

		EchoType outputClass = r.getEvalType();
		for(String name : inputIndices.keySet()) { 
			int index = inputIndices.get(name);
			if(listAll || !inputs.containsKey(name)) {
				EchoType verbParamClass = verb.getParameterClasses()[index];
				if(outputClass.isSubType(verbParamClass)) { 
					names.add(name);
				}
			}
		}

		return new Vector<String>(names);
	}

	public void reset() { 
		proc = null;
	}

	public void evaluate() throws EchoException {
		if(proc == null) { 
			Set<String> missingParams = findMissingParams();
			Set<String> missingInputs = findMissingInputs();

			if(!missingParams.isEmpty()) {
				System.err.println("\n" + verb.getClass().getName() + " => Parameters expected: ");
				String[] na = verb.getParameterNames();
				for(int i = 0; na != null && i < na.length; i++) { System.err.println("\t" + i + ": " + na[i]); }
				throw new EchoException("Missing parameters: " + missingParams); 
			}
			if(!missingInputs.isEmpty()) { throw new EchoException("Missing Inputs: " + missingInputs); }

			Vector<InputQueue> inpqs = new Vector<InputQueue>();
			String[] inputNames = verb.getInputNames();
			EchoType[] inputClasses = verb.getInputClasses();
			
			for(int i = 0; inputNames != null && i < inputNames.length; i++) { 
				Reverb r = inputs.get(inputNames[i]);
				SelfDescribingVerb v = r.getVerb();
				
				r.evaluate();
				
				boolean expand = true;
				if((v instanceof Generator || v instanceof Expander || v instanceof SelfDescribingOutput) && 
						!shouldExpand(v.getOutputClass(), inputClasses[i])) { 
					expand = false;
				}

				InputQueue q = r.getQueue(expand);
				inpqs.add(q);
			}

			proc = new EchoProcessor(verb, inpqs);

			Map<String,Object> values = new HashMap<String,Object>();
			for(String paramName : params.keySet()) { 
				values.put(paramName, params.get(paramName).evaluate());
			}

			proc.init(values);
		}
	}
	
	public static boolean shouldExpand(EchoType output, EchoType input) { 
		if(output.equals(input)) { return true; }
		if(input instanceof SequenceType) { 
			SequenceType st = (SequenceType)input;
			/*
			if(st.getInnerType().equals(output)) { 
				return false;
			}
			*/
			
			if(output.isSubType(st.getInnerType())) { 
				return false;
			}
		}
		return true;
	}

	public InputQueue getQueue(boolean expand) { 
		return proc.createOutput(expand);
	}

	public EchoProcessor getProcessor() { return proc; }
}
