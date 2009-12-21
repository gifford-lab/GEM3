package edu.mit.csail.cgs.echo;

import java.util.*;

import edu.mit.csail.cgs.ewok.types.SelfDescribingOutput;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.Closeable;

/**
 * @author Timothy Danford
 *
 * EchoProcessor is the glue which combines verbs together in an 
 * actual, functioning arrangement.  It is a replacement for the 
 * role that the java.util.Iterator implementations (ExpanderIterator,
 * MapperIterator, etc) play in the original ewok framework.
 *
 * It is an asynchronous data manipulator, whereas the Iterators
 * are synchronous data manipulators.
 */
public class EchoProcessor {
	
	private Vector<InputQueue> inputs;
	private Vector<EchoQueue> outputs;
	private SelfDescribingVerb verb;

    private boolean hasNonExpandingQueues, isFinishedGenerator;
	
	public EchoProcessor(SelfDescribingVerb v, Vector<InputQueue> qs) { 
		verb = v;
		inputs = new Vector<InputQueue>(qs);
		outputs = new Vector<EchoQueue>();
		hasNonExpandingQueues = false;
        isFinishedGenerator = false;
	}
	
	public EchoQueue createOutput() { 
		return createOutput(true);
	}
	
	public EchoQueue createOutput(boolean exp) {
		EchoQueue q = new EchoQueue(exp);
		outputs.add(q);
		hasNonExpandingQueues = hasNonExpandingQueues || !exp;		
		return q;
	}
	
	public void init(Map<String,Object> params) {
		
		for(EchoQueue q : outputs) { q.reset(); }
		verb.init(params);
		
        /*
		if(verb instanceof Generator) {
			inputs.clear();
			Generator gen = (Generator)verb;
			Iterator itr = gen.execute();
			IteratorInputQueue q = new IteratorInputQueue(itr);
			inputs.add(q);
		}
        */
	}
	
	public SelfDescribingVerb getVerb() { return verb; }
	
	public void finish() { 
        System.out.println("\t\tFinishing : " + verb.getClass().getName());
        
        if(verb instanceof Sink) { 
			Sink s = (Sink)verb;
			s.finish();
		}

		if(verb instanceof MultiSink) { 
			MultiSink s = (MultiSink)verb;
			s.finish();
		}

		if(verb instanceof SelfDescribingOutput) { 
			SelfDescribingOutput outputSink = (SelfDescribingOutput)verb;
			Collection vals = outputSink.getValues();
			
			for(EchoQueue q : outputs) {
				if(q.isExpanding()) { 
					for(Object val : vals) { 
						q.addValue(val);
					}
				} else { 
					q.addValue(vals.iterator());
				}
			}
		}
        
        int c = 0;
		for(EchoQueue q : outputs) { q.finish(); c++; }
        
        System.out.println("\t\tFinished " + c + " queues.");
	}
	
	public boolean isFinished() {
        if((verb instanceof Generator)) { 
            return isFinishedGenerator;

        } else if(verb instanceof MultiSink) { 
			// MultiSinks have a different semantics than any other verb -- 
			// they are only finished when *all* their input queues are finished.
			for(InputQueue inp : inputs) { 
				if(!inp.isFinished()) { 
					return false; 
				}
			}
			return true;
		
        } else { 
			// Normal verbs are finished when *any* of their input queues
			// are finished.
			for(InputQueue inp : inputs) { 
				if(inp.isFinished()) { 
					return true; 
				}
			}
			return false;
		}
	}
	
	public boolean isReady() {
		// isFinished() should be checked *before* calling this method.
        
		if(verb instanceof Generator) { 
            return true;
            
        } else if(verb instanceof MultiSink) { 
			// Once again, MultiSinks have different semantics here.  A multisink 
			// is 'ready' if *any* of its input queues are non-empty.
			for(InputQueue inp : inputs) { 
				if(!inp.isEmpty()) { 
					return true; 
				}
			}
			return false;
		} else { 
			// Normal verbs are only ready if *all* their input queues are non-empty.
			for(InputQueue inp : inputs) { 
				if(inp.isEmpty()) { 
					return false; 
				}
			}
			return true;
		}
	}
	
	private Collection convertIterator(Iterator itr) { 
		LinkedList lst = new LinkedList();
		while(itr.hasNext()) { lst.addLast(itr.next()); }
		return lst;
	}
	
    private void outputIterator(Iterator expanded) { 
        if(hasNonExpandingQueues) { 
            Collection expValues = convertIterator(expanded);
            
            for(EchoQueue output : outputs) { 
                if(output.isExpanding()) { 
                    for(Object v : expValues) { 
                        output.addValue(v);
                    }
                } else { 
                    output.addValue(expValues.iterator());
                }
            }
        } else { 
            while(expanded.hasNext()) { 
                Object outval = expanded.next();
                for(EchoQueue output : outputs) { 
                    output.addValue(outval);
                }
            }
        }
    }
    
	public void process() { 
        if(verb instanceof Generator) {
            Generator gen = (Generator)verb;
            Iterator expanded = gen.execute();
            outputIterator(expanded);
            
            isFinishedGenerator = true;            
		}
        
        if(verb instanceof BiCombiner) { 
            Object inval1 = inputs.get(0).getFirstValue();
            Object inval2 = inputs.get(1).getFirstValue();
            Object outval = ((BiCombiner)verb).execute(inval1, inval2);
            for(EchoQueue output : outputs) { 
                output.addValue(outval);
            }                       
        }
        
        if(verb instanceof Sink) { 
            Object inval = inputs.get(0).getFirstValue();
            Sink s = (Sink)verb;
            s.consume(inval);            
        }

		if(verb instanceof MultiSink) { 
			String[] inputNames = verb.getInputNames();
			MultiSink ms = (MultiSink)verb;
			for(int i = 0; i < inputNames.length; i++) { 
				if(!inputs.get(i).isEmpty()) { 
					ms.consume(inputNames[i], inputs.get(i).getFirstValue()); 
				}
			}
		}

		if(verb instanceof Expander) { 
			Object val = inputs.get(0).getFirstValue();
			Iterator expanded = ((Expander)verb).execute(val);
			outputIterator(expanded);
		}
		
		if(verb instanceof Mapper) { 
            Object inval = inputs.get(0).getFirstValue();
            Object outval = ((Mapper)verb).execute(inval);
            for(EchoQueue output : outputs) { 
                output.addValue(outval);
            }           
		}
		
		if(verb instanceof Filter && !(verb instanceof Mapper)) {
			Object inval = inputs.get(0).getFirstValue();
			Object outval = ((Filter)verb).execute(inval);
			if(outval != null) { 
				for(EchoQueue output : outputs) { 
					output.addValue(outval);
				}			
			}
		}
		
		if(verb instanceof Combiner) { 
			Object inval1 = inputs.get(0).getFirstValue();
			Object inval2 = inputs.get(1).getFirstValue();
			Object outval = ((Combiner)verb).execute(inval1, inval2);
			for(EchoQueue output : outputs) { 
				output.addValue(outval);
			}						
		}
	}
}
