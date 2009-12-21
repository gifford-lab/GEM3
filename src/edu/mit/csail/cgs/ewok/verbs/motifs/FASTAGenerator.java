package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.Iterator;
import java.util.Map;
import java.io.File;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.ewok.verbs.Generator;
import edu.mit.csail.cgs.ewok.types.*;

public class FASTAGenerator implements Generator<Pair<String,String>> { 

	public static final EchoType[] inputTypes = null; 
	public static final String[] inputNames = null;
	public static final EchoType outputType = new PairType(new ClassType(String.class), new ClassType(String.class));
	public static final EchoType[] paramTypes = { new ClassType(File.class) };
	public static final String[] paramNames = { "File" };

	private File file;
	private FASTALoader loader;

	public FASTAGenerator() { file = null; loader = new FASTALoader(); }
	public FASTAGenerator(File f) { file = f; loader = new FASTALoader();  }

	public void init(Map<String,Object> params) { 
		file = (File)params.get(paramNames[0]);
	}

	public EchoType[] getInputClasses() { return inputTypes; }
	public String[] getInputNames() { return inputNames; }
	public EchoType getOutputClass() { return outputType; }
	public EchoType[] getParameterClasses() { return paramTypes; }
	public String[] getParameterNames() { return paramNames; }

	public Iterator<Pair<String,String>> execute() { 
		return loader.execute(file);
	}
}
