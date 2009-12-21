package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;

public class Parser { 

	private File file;

	public Parser(File f) { 
		file = f;
	}

	public DirectedGraph parseDirectedGraph() throws IOException { 
		DirectedGraph dg = new DirectedGraph();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = null;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			if(line.length() > 0) { 
				String[] array = line.split("\\s+");
				if(array.length == 3 && array[1].equals("-->")) { 
					String v1 = array[0], v2 = array[2];
					dg.addVertex(v1);
					dg.addVertex(v2);
					dg.addEdge(v1, v2);
				} else { 
					System.err.println("Unknown line format: \"" + line + "\"");
				}
			}
		}

		br.close();
		return dg;
	}
}
