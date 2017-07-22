package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;

public class ChIPseqSignal {

	public static void main(String[] args) {
		final int MAXREAD = 1000000;
		final boolean dev = false;
		Genome genome = CommonUtils.parseGenome(args);
		String outPrefix = Args.parseString(args, "out", "out");
		String infoFile = Args.parseString(args, "info", null);
		if (infoFile == null){
			System.out.println("Need --info: information about the data.");
			System.exit(-1);
		}
		
		ArrayList<String> text = CommonUtils.readTextFile(infoFile);
		ArrayList<String> expts = new ArrayList<String>();
		ArrayList<String> data_locations = new ArrayList<String>();
		ArrayList<Boolean> isReadDB = new ArrayList<Boolean>();
		ArrayList<Integer> radius = new ArrayList<Integer>();
		// expt name | expand_radius | readdb name
		for (String txt: text){
			if ( ! (txt.equals("") || txt.startsWith("#")) ){
				String[] f = txt.trim().split("\t");
				if (f.length == 4){		// new format with data format
					expts.add(f[0]);
					radius.add(Integer.parseInt(f[1]));
					isReadDB.add(f[2].equalsIgnoreCase("RDB"));
					data_locations.add(f[3]);
				}
				if (f.length == 3){		// old format
					expts.add(f[0]);
					radius.add(Integer.parseInt(f[1]));
					isReadDB.add(true);
					data_locations.add(f[2]);
				}
			}
		}
		
		String regionFile = Args.parseString(args, "coords", null);
		if (regionFile == null){
			System.out.println("Need --coords: the coordinates for getting ChIP-seq occupancy signals.");
			System.exit(-1);
		}
		ArrayList<Point> all_regions = new ArrayList<Point>();
		text = CommonUtils.readTextFile(regionFile);
		for (String t: text){
			String[] f = t.split("\t");
			if (f[0].contains("#") || f[0].contains("Position"))
				continue;
			all_regions.add(Point.fromString(genome, f[0]));
		}
		all_regions.trimToSize();

		int[][]signals = new int[all_regions.size()][expts.size()];
		for (int i=0;i<expts.size();i++){
			if (isReadDB.get(i)){		// readdb
				String readdb_name = data_locations.get(i);
				int rr = radius.get(i);
				
				List<ChipSeqLocator> rdbexpts = new ArrayList<ChipSeqLocator>();
				String[] pieces = readdb_name.trim().split(";");
	            if (pieces.length == 2) {
	            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1]));
	            } else if (pieces.length == 3) {
	            	rdbexpts.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
	            } else {
	                throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + readdb_name);
	            }
	            DeepSeqExpt ip = new DeepSeqExpt(genome, rdbexpts, "readdb", -1);
	            ReadCache ipCache = new ReadCache(genome, expts.get(i), null, null);
	            
				// cache sorted start positions and counts of all positions
				long tic = System.currentTimeMillis();
				System.err.print("Loading "+ipCache.getName()+" data from ReadDB ... \t");
				List<String> chroms = genome.getChromList();
				if (dev){
					chroms = new ArrayList<String>();
					chroms.add("19");
				}
				// load  data into cache by chroms or smaller chunks.
				for (String chrom: chroms ){
					int length = genome.getChromLength(chrom);
					Region wholeChrom = new Region(genome, chrom, 0, length-1);
					int count = ip.countHits(wholeChrom);
					ArrayList<Region> chunks = new ArrayList<Region>();
					// if there are too many reads in a chrom, read smaller chunks
					if (count>MAXREAD){
						int chunkNum = count/MAXREAD*2+1;
						int chunkLength = length/chunkNum;
						int start = 0;
						while (start<=length){
							int end = Math.min(length, start+chunkLength-1);
							Region r = new Region(genome, chrom, start, end);
							start = end+1;
							chunks.add(r);
						}
					}else
						chunks.add(wholeChrom);
	
					for (Region chunk: chunks){
						Pair<ArrayList<Integer>,ArrayList<Float>> hits = ip.loadStrandedBaseCounts(chunk, '+');
						ipCache.addHits(chrom, '+', hits.car(), hits.cdr());
						hits = ip.loadStrandedBaseCounts(chunk, '-');
						ipCache.addHits(chrom, '-', hits.car(), hits.cdr());
					}
				} // for each chrom
	
				ipCache.populateArrays(true);
				ip.closeLoaders();
				ip=null;
				System.gc();
				ipCache.displayStats();
				System.out.println(CommonUtils.timeElapsed(tic));
	            
				// now get the data from the cache
	            for (int j=0;j<all_regions.size();j++){
	            	Region region = all_regions.get(j).expand(rr);
	            	List<StrandedBase> bases = ipCache.getStrandedBases(region, '+');
	            	bases.addAll(ipCache.getStrandedBases(region, '-'));
	            	signals[j][i] = (int)StrandedBase.countBaseHits(bases);
	            }
			}
			else{		// BAM
				List<File> files = new ArrayList<File>();
				files.add(new File(data_locations.get(i)));
				DeepSeqExpt ip = new DeepSeqExpt(genome, files, true, "SAM", -1);
				int rr = radius.get(i);
				// now get the data from the BAM file
	            for (int j=0;j<all_regions.size();j++){
	            	Region region = all_regions.get(j).expand(rr);
	            	signals[j][i] = ip.countHits(region);
	            }
	            ip.closeLoaders();
			}
		}
		StringBuilder sb = new StringBuilder("#Site\t");
		for (int i=0;i<expts.size();i++){
			sb.append(expts.get(i)).append("\t");
		}
		CommonUtils.replaceEnd(sb, '\n');
		for (int j=0;j<all_regions.size();j++){
			sb.append(all_regions.get(j).toString()).append("\t");
			for (int i=0;i<expts.size();i++){
				sb.append(signals[j][i]).append("\t");
			}
			CommonUtils.replaceEnd(sb, '\n');
		}
		CommonUtils.writeFile("0_Read_signals."+outPrefix+".txt", sb.toString());
	}
}
