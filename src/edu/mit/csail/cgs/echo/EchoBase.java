/*
 * Created on Feb 16, 2007
 */
package edu.mit.csail.cgs.echo;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.echo.components.ChromRegionWrapper;
import edu.mit.csail.cgs.echo.components.TabbedFileRegionSink;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.types.ValueWrapper;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * @author Timothy Danford
 */
public class EchoBase {
    
    public static void main(String[] args) { 
        try {
            
            if(!Reverb.isSubclass(Gene.class, Region.class)) { 
                System.err.println("Gene doesn't appear to be a sub-class of Region");
            }
            
            Genome genome = Organism.findGenome("sacCer1");
            EchoConstant c1 = new EchoConstant(new ValueWrapper(genome));
            EchoConstant c2 = new EchoConstant(new ValueWrapper(
                    "C:\\Documents and Settings\\tdanford\\Desktop\\test_echo.txt"));
            
            Reverb r1 = new Reverb(new ChromRegionWrapper());
            r1.setParam("Genome", c1);
            
            Reverb r2 = new Reverb(new RefGeneGenerator(genome));
            r2.setParam("Genome", c1);
            r2.setInput("regions", r1);
            
            Reverb fs = new Reverb(new TabbedFileRegionSink());
            fs.setParam("Filename", c2);
            fs.setInput("regions", r2);
            
            EchoBase base = new EchoBase();
            base.addConstant(c1);
            base.addConstant(c2);
            base.addReverb(r1);
            base.addReverb(r2);
            base.addReverb(fs);
            
            Vector<Reverb> targets = new Vector<Reverb>();
            targets.add(fs);
            
            base.runEchoProgram(targets);
            
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (EchoException e) {
            e.printStackTrace();
        }
    }
	
	private Vector<Reverb> reverbs;
	private Vector<EchoConstant> constants;

	public EchoBase() { 
		reverbs = new Vector<Reverb>();
		constants = new Vector<EchoConstant>();
    }
    
    public void clear() { 
        reverbs.clear();
        constants.clear();
    }
    
    public void remove(Reverb r) { 
    	reverbs.remove(r);
    }
    
    public void remove(EchoConstant c) { 
    	constants.remove(c);
    }
    
    public void addEchoObject(Object obj) { 
        if(obj instanceof Reverb) { addReverb((Reverb)obj); } 
        if(obj instanceof EchoConstant) { addConstant((EchoConstant)obj); }
    }

	public void addReverb(Reverb r) { 
		reverbs.add(r);
	}

	public void addConstant(EchoConstant c) { 
		constants.add(c);
	}
    
    public void runEchoProgram() throws EchoException { 
        runEchoProgram(reverbs);
    }

	public void runEchoProgram(Vector<Reverb> targets) throws EchoException {
        for(int i = 0; i < constants.size(); i++) { 
            constants.get(i).init();
        }
        
        for(Reverb r : reverbs) { 
        	r.reset();
        }
        
        for(Reverb r : targets) { 
            try { 
                r.evaluate();
                System.out.println("Evaluated: " + r.getVerb().getClass().getName());
            } catch(EchoException ee) { 
                System.err.println("Couldn't evaluate reverb: " + r.getVerb().getClass().getName());
                System.err.println("Message: " + ee.getMessage());
            }
        }
        
		Vector<EchoProcessor> procs = new Vector<EchoProcessor>(); 

        for(Reverb r : reverbs) { 
			EchoProcessor proc = r.getProcessor();
			if(proc != null) { procs.add(proc); }
		}
		
		EchoScheduler sched = new EchoScheduler(procs);
		Thread t = new Thread(sched);
		t.start();
	}
}
