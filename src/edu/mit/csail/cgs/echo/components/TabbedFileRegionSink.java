
package edu.mit.csail.cgs.echo.components;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;

public class TabbedFileRegionSink implements SelfDescribingVerb, Sink<Region> { 

	private static EchoType[] paramClassArray = { new ClassType(File.class) };
	private static String[] paramNameArray = { "File" };

    private static final EchoType[] inputClasses = { new ClassType(Region.class) };
    private static final String[] inputNames = { "Regions" };
    
	private Map<String,Object> params;
    private PrintStream ps;
	
	public TabbedFileRegionSink() { 
		params = new HashMap<String,Object>();
        ps = null;
	}
	
	public void init(Map pms) { 
		params = new HashMap<String,Object>(pms); 
		init();
	}
    
    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }
	
	public EchoType[] getParameterClasses() { return paramClassArray; }
	public String[] getParameterNames() { return paramNameArray; }
    
    public EchoType getOutputClass() { return null; }

	public void init() { 
		File output = (File)params.get(paramNameArray[0]);
		try { 
			ps = new PrintStream(new FileOutputStream(output));
		} catch(IOException ie) { 
			ps = null;
		}
	}

	public void finish() { 
		if(ps != null) { ps.close(); }
	}

	public void consume(Region r) { 
		PrintStream sps = ps != null ? ps : System.out;
		sps.print(r.getChrom() + "\t" + r.getStart() + "\t" + r.getEnd());
		if(r instanceof NamedRegion) { 
			sps.print("\t" + ((NamedRegion)r).getName());
		}
		
		if(r instanceof Gene) { 
			sps.print("\t" + ((Gene)r).getID());
		}
		sps.println();
	}
	
	public void consume(Iterator<Region> itr) { 
		init();
		while(itr.hasNext()) { 
			Region r = (Region)itr.next();
			consume(r);
		}
		finish();
	}
}

