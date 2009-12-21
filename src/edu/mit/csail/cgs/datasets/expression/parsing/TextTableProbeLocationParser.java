package edu.mit.csail.cgs.datasets.expression.parsing;

import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.expression.ExpressionInserter;
import edu.mit.csail.cgs.datasets.expression.ExpressionLoader;
import edu.mit.csail.cgs.datasets.expression.ProbePlatform;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;

public class TextTableProbeLocationParser {
	
	public static void main(String[] args) { 
		try {
			Genome g = Organism.findGenome(args[0]);
			File inputfile = new File(args[1]);
			String platformName = args[2];
			
			TextTableProbeLocationParser parser = 
				new TextTableProbeLocationParser(g, inputfile);
			
			System.out.println("# Probes: " + parser.getNumProbes());
			System.out.println("# Locations: " + parser.getNumRegions());
			
			ExpressionLoader loader = new ExpressionLoader();
			ProbePlatform pp = loader.loadPlatform(platformName);
			
			System.out.println("Platform ID: " + pp.getDBID());
			
			ExpressionInserter inserter = new ExpressionInserter();
			parser.insertIntoDB(pp, inserter);
			
			System.out.println("Finished inserting probe locations.");
			
			inserter.insertProbePlatformGenomeMapping(pp, g);
			System.out.println("Inserted probe_platform_to_genome entry.");
			
			inserter.close();
			loader.close();
			
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	private static Pattern locpatt = Pattern.compile("chr([\\w\\d]+):(\\d+)-(\\d+)");
	
	private Map<String,Set<Region>> locations;
	private Genome genome;
	private int regionCount;

	public TextTableProbeLocationParser(Genome g, File f) throws IOException { 
		locations = new HashMap<String,Set<Region>>();
		genome = g;
		regionCount = 0;
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		int lineCount = 0;
		while((line = br.readLine()) != null) { 
			try { 
				String[] array = line.split("\\s+");
				String pname = array[0];
				String locstring = array[1];
				Matcher m = locpatt.matcher(locstring);
				if(m.matches()) { 
					String chrom = m.group(1);
					int start = Integer.parseInt(m.group(2));
					int end = Integer.parseInt(m.group(3));
					int chromID = genome.getChromID(chrom);
					Region region = null;
					
					if(start <= end) { 	
						region = new StrandedRegion(genome, chrom, start, end, '+');
					} else { 
						region = new StrandedRegion(genome, chrom, end, start, '-');
					}
					if(!locations.containsKey(pname)) { 
						locations.put(pname, new HashSet<Region>());
					}
					locations.get(pname).add(region);
					regionCount += 1;
				}
			} catch(Exception e) { 
				System.out.println("An error \"" + e.getMessage() + "\" was encountered for " + 
						" line \"" + line + "\" (" + lineCount + ")");
			}
			lineCount += 1;
		}
		br.close();
	}
	
	public int getNumRegions() { return regionCount; }
	public int getNumProbes() { return locations.size(); }
	
	public void insertIntoDB(ProbePlatform pp, ExpressionInserter inserter) throws SQLException { 
		inserter.insertProbeLocations(locations, pp);
	}
}
