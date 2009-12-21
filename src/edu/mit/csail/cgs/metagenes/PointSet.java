package edu.mit.csail.cgs.metagenes;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

/**
 * @author tdanford
 */
public class PointSet<PointClass extends Point> {
	
	private Vector<PointClass> points;

	public PointSet() { 
		points = new Vector<PointClass>();
	}
	
	public PointSet(Genome g, File f) throws IOException { 
		Mapper<String,Point> pointMapper = new DefaultPointParser(g);
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null) { 
			PointClass pc = (PointClass)pointMapper.execute(line);
			points.add(pc);
		}
		br.close();
	}
	
	public PointSet(Mapper<String,PointClass> pointMapper, File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null) { 
			PointClass pc = pointMapper.execute(line);
			points.add(pc);
		}
		br.close();
	}
	
	public void addPoint(PointClass p) { 
		points.add(p);
	}
	
	public int size() { return points.size(); }
	public PointClass point(int i) { return points.get(i); }

	public static Pattern pointPattern = Pattern.compile("(?:chr)?([^:]+):(\\d+)");
	
	public static class DefaultPointParser implements Mapper<String,Point> { 
		
		private Genome genome;
		public DefaultPointParser(Genome g) { genome = g; }
		
		public Point execute(String line) { 
			String[] a = line.split("\\s+");
			Matcher m = pointPattern.matcher(a[0]);
			if(!m.matches()) { throw new IllegalArgumentException(line); }
			String chrom = m.group(1);
			int loc = Integer.parseInt(m.group(2));
			return new Point(genome, chrom, loc);
		}
	}
	
}
