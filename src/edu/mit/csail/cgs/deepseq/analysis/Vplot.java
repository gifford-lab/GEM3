package edu.mit.csail.cgs.deepseq.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.analysis.TFBS_SpaitialAnalysis.Site;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;

public class Vplot {

	public static void main(String[] args) {
		Genome genome = CommonUtils.parseGenome(args);
		int radius = Args.parseInteger(args, "x_width", 1000)/2;
		int maxReadLength = Args.parseInteger(args, "y_length", 500);
		int xPixel = Args.parseInteger(args, "x_pixel", 500);
		int yPixel = Args.parseInteger(args, "y_pixel", 50);
		int binL = maxReadLength/yPixel;
		int binW = radius*2/xPixel;
		int min = Args.parseInteger(args, "min", 30);
		String outName = Args.parseString(args, "out", "out");
		ArrayList<String> lines = CommonUtils.readTextFile(Args.parseString(args, "atac", null));
		ArrayList<Region> reads = new ArrayList<Region>();
		for (String l: lines){
			String f[] = l.split("\t");
			reads.add(new Region(genome, f[0].replace("chr", ""), Integer.parseInt(f[1]), Integer.parseInt(f[5])));
		}
		reads.trimToSize();
		Collections.sort(reads);
		ArrayList<Point> starts = new ArrayList<Point>();
		for (Region r: reads)
			starts.add(r.startPoint());
		starts.trimToSize();
		
//		String gemFileName = Args.parseString(args, "gem", null);
//		List<GPSPeak> gpsPeaks =null;
//	    try{
//	    	gpsPeaks= GPSParser.parseGPSOutput(gemFileName, genome);
//	    }
//	    catch(IOException e){
//	    	System.err.println("Error reading/parsing GPS file: "+gemFileName);
//	    }
	    
		ArrayList<StrandedPoint> sps = new ArrayList<StrandedPoint>();
		lines = CommonUtils.readTextFile(Args.parseString(args, "piq", null));
		for (String l: lines){
			if (l.contains("score"))
				continue;
			String[] q = l.split(",");
			sps.add(new StrandedPoint(genome, q[1].replaceAll("\"", "").replaceAll("chr", ""), Integer.parseInt(q[2]), '+'));
		}
		lines = CommonUtils.readTextFile(Args.parseString(args, "piqrc", null));
		for (String l: lines){
			if (l.contains("score"))
				continue;
			String[] q = l.split(",");
			sps.add(new StrandedPoint(genome, q[1].replaceAll("\"", "").replaceAll("chr", ""), Integer.parseInt(q[2]), '-'));
		}
	    
	    int maxWidth = 0;
	    int numSites = 0;
	    StringBuilder mSB = new StringBuilder();	// matrix
	    StringBuilder iSB = new StringBuilder();	// info
	    mSB.append("#").append(yPixel).append("\t").append(xPixel).append("\n");
	    iSB.append("#y_pixel:").append(yPixel).append("\tx_pixel:").append(xPixel).append("\n");
		for (StrandedPoint p: sps){
			Region r = p.expand(radius);
			ArrayList<Region> inRange = new ArrayList<Region>();
			int id = Collections.binarySearch(starts, p);
			if (id < 0) // if key not found
				id = -(id + 1);
			if (id == starts.size())
				continue;
			for (int i=id; i>=0; i--){
				if (r.contains(starts.get(i)))
					inRange.add(reads.get(i));
				else
					break;
			}
			for (int i=id+1; i<starts.size(); i++){
				if (r.contains(starts.get(i)))
					inRange.add(reads.get(i));
				else
					break;
			}
			if (inRange.size()<min)
				continue;			
			
			// Passed the minimun read count cutoff
			numSites++;
			int[][] matrix = new int[xPixel][yPixel];
//			System.out.println("\n============");
			System.out.println(p);
			int numPoint = inRange.size();
			for (Region ri: inRange){
				int offset = (p.getStrand()=='-'?-1:1) * (ri.startPoint().offset(p));
				if (offset==radius){
					numPoint--;
					continue;
				}
				int length = ri.getWidth();
				if (maxWidth<length)
					maxWidth = length;
				if (length>=maxReadLength)
					length = maxReadLength-1;
				matrix[(offset+radius)/binW][length/binL]++;
//				if (maxWidth<ri.getWidth())
//					maxWidth = ri.getWidth();
			}
			
			//output
			iSB.append(numSites).append("\t")
				.append(p.toString()).append("\t")
				.append(r.toString()).append("\t").append(numPoint).append("\n");
			for (int i=0; i<matrix.length; i++){
				for (int j=0; j<matrix[i].length;j++){
					mSB.append(matrix[i][j]).append("\t");
				}
			}
			CommonUtils.replaceEnd(mSB, '\n');
		}
		CommonUtils.writeFile(outName+".data.txt", mSB.toString());
		CommonUtils.writeFile(outName+".info.txt", iSB.toString());
		
		System.out.println("\nTotal number of sites: "+numSites);
		System.out.println("Longest read length: "+maxWidth);
	}

}
