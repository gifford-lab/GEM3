package edu.mit.csail.cgs.metagenes;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.ChipSeqExpander;
import edu.mit.csail.cgs.metagenes.swing.MetaFrame;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class MetaMaker {
	private static boolean batchRun = false;
	private static boolean cluster = false;
	public static TreeMap<String, Color> map = new TreeMap<String, Color>();
	public static TreeMap<String, Color> getColorMap(){
		if (map.isEmpty()){
			map.put("red",Color.red);
			map.put("darkred",Color.red.darker());
			map.put("green",Color.green);
			map.put("darkgreen",Color.green.darker());
			map.put("black",Color.black);
			map.put("blue",Color.blue);
			map.put("darkblue",Color.blue.darker());
			map.put("cyan",Color.cyan);
			map.put("darkcyan",Color.cyan.darker());
			map.put("gray",Color.gray);
			map.put("magenta",Color.magenta);
			map.put("darkmagenta",Color.magenta.darker());
			map.put("orange",Color.orange);
			map.put("darkorange",Color.orange.darker());
			map.put("pink",Color.pink);
			map.put("darkpink",Color.pink.darker());
			map.put("darkgray",Color.darkGray);
		}
		return map;
	}
	
	public static void main(String[] args) {
		try {
			if(args.length < 2){ printError();}
			
			Genome gen = CommonUtils.parseGenome(args);
			
			// color of the plots
			getColorMap();
			
			if(Args.parseFlags(args).contains("showcolor")){
				System.out.println("Available colors:");
				for (String s:map.keySet())
					System.out.print(s+" ");
				System.out.println();
				System.exit(-1);
			}
			Color c = Color.blue;
			String newCol = Args.parseString(args, "color", "notFound");
			if (map.containsKey(newCol)){
				c = map.get(newCol);
			}
			else{
				System.out.println("The specified color is not defined!");
				System.out.println("Available colors:");
				for (String s:map.keySet())
					System.out.print(s+" ");
				System.out.println();
				System.exit(-1);
			}
			
			
			double peakMax = Args.parseDouble(args, "peakMax", 1.0);
			int winLen = Args.parseInteger(args,"win", 10000);
			int bins = Args.parseInteger(args,"bins", 100);
			int readExt = Args.parseInteger(args,"readext", 0);			// read extention, -1 for only five prime
			int readShift = Args.parseInteger(args,"readshift", 0);			// read extention, -1 for only five prime			
			double lineMin = Args.parseDouble(args,"linemin", 0);		// Min Color Value
			double lineMax = Args.parseDouble(args,"linemax", 100);		// Max Color Value
			int lineThick = Args.parseInteger(args,"linethick", 1);
			double pbMax = Args.parseDouble(args,"pbMax", 100);			// omit (NOT truncate the count) a position if the per base read count is higher 
			char strand = Args.parseString(args, "strand", "/").charAt(0);		// read strand
			boolean drawColorBar = !Args.parseFlags(args).contains("nocolorbar");
			String profilerType = Args.parseString(args, "profiler", "simplechipseq");	
			List<String> expts = (List<String>) Args.parseStrings(args,"expt");
			List<String> files = (List<String>) Args.parseStrings(args,"file");
			List<String> backs = (List<String>) Args.parseStrings(args,"back");
			List<String> peakFiles = (List<String>)Args.parseStrings(args, "peaks");
			String format = Args.parseString(args, "format", "SAM");
			String outName = Args.parseString(args, "out", "meta");
			if(Args.parseFlags(args).contains("batch")){batchRun=true;}
			if(Args.parseFlags(args).contains("cluster")){cluster=true;}
			
			if(gen==null || (expts.size()==0 && files.size()==0)){printError();}
	
			BinningParameters params = new BinningParameters(winLen, bins);
			System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
		
			PointProfiler profiler=null;
			boolean normalizeProfile=false;
			if(profilerType.equals("simplechipseq") || profilerType.equals("fiveprime")){
				//normalizeProfile=true;
				List<ChipSeqLocator> exptlocs = Args.parseChipSeq(args,"expt");
//				DeepSeqExpt ip = new DeepSeqExpt(gen, exptlocs, "readdb", -1);
				
				ArrayList<ChipSeqExpander> exptexps = new ArrayList<ChipSeqExpander>();
				for(ChipSeqLocator loc : exptlocs){
					System.out.println(loc.getExptName()+"\t"+loc.getReplicateString()+"\t"+loc.getAlignName());
					exptexps.add(new ChipSeqExpander(loc));
				}
				if (!files.isEmpty()){
					List<File> fs = new ArrayList<File>();
					for (String s: files)
						fs.add(new File(s));
					exptexps.add(new ChipSeqExpander(gen, fs, format));
				}
				System.out.println("Loading data...");
				if(profilerType.equals("fiveprime"))
					readExt = -1;
				profiler = new SimpleChipSeqProfiler(params, exptexps, readExt, readShift, pbMax,strand);
				
			}else if(profilerType.equals("simplechiapet")) {
				List<ChipSeqLocator> exptlocs = Args.parseChipSeq(args,"expt");
				ArrayList<ChipSeqExpander> exptexps = new ArrayList<ChipSeqExpander>();
				for(ChipSeqLocator loc : exptlocs){
					System.out.println(loc.getExptName()+"\t"+loc.getAlignName());
					exptexps.add(new ChipSeqExpander(loc,true));
				}
				System.out.println("Loading data...");
				profiler = new SimpleChipSeqProfiler(params, exptexps, readExt, readShift, pbMax,strand);
			}else if(profilerType.equals("chipseq5prime")){
				List<ChipSeqLocator> exptlocs = Args.parseChipSeq(args,"expt");
				ArrayList<ChipSeqExpander> exptexps = new ArrayList<ChipSeqExpander>();
				for(ChipSeqLocator loc : exptlocs){
					exptexps.add(new ChipSeqExpander(loc));
				}
				System.out.println("Loading data...");
				profiler = new ChipSeq5PrimeProfiler(params, exptexps, strand);
			}else if(profilerType.equals("chipseq")){
				normalizeProfile=true;
				ArrayList<ChipSeqLocator> exptlocs = (ArrayList<ChipSeqLocator>) Args.parseChipSeq(args,"expt");
				ArrayList<ChipSeqLocator> backlocs = backs.size()==0 ? null : (ArrayList<ChipSeqLocator>) Args.parseChipSeq(args,"back");
				System.out.println("Loading data...");
				profiler = new ChipSeqProfiler(params, gen, exptlocs,backlocs, 32, readExt);
			}else if(profilerType.equals("chipseqz")){
				normalizeProfile=true;
				ArrayList<ChipSeqLocator> exptlocs = (ArrayList<ChipSeqLocator>) Args.parseChipSeq(args,"expt");
				ArrayList<ChipSeqLocator> backlocs = backs.size()==0 ? null : (ArrayList<ChipSeqLocator>) Args.parseChipSeq(args,"back");
				System.out.println("Loading data...");
				profiler = new ChipSeqProfiler(params, gen, exptlocs,backlocs, 32, readExt, true);
			}else if(profilerType.equals("chipchip") || profilerType.equals("chipchipip") || profilerType.equals("chipchipwce")){
				normalizeProfile=true;
				ArrayList<ChipChipLocator> exptlocs = (ArrayList<ChipChipLocator>) Args.parseChipChip(gen, args, "expt");
				if(exptlocs.size()>0){
					System.out.println("Loading data...");
					if(profilerType.equals("chipchipip"))
						profiler = new ChipChipProfiler(params, gen, exptlocs.get(0), true, false);
					else if(profilerType.equals("chipchipwce"))
						profiler = new ChipChipProfiler(params, gen, exptlocs.get(0), false, true);
					else
						profiler = new ChipChipProfiler(params, gen, exptlocs.get(0));
				}
			}
			
			if(batchRun){
				System.out.println("Batch running...");
				
				if(peakFiles.size()>=1){
					MetaNonFrame nonframe = new MetaNonFrame(gen, params, profiler, normalizeProfile, peakMax);
					nonframe.setColor(c);
					nonframe.setDrawColorBar(drawColorBar);
					MetaProfileHandler handler = nonframe.getHandler();
					if(peakFiles.size()==1){
						System.out.println("Single set mode...");
						String peakFile = peakFiles.get(0);
						Vector<Point> points = nonframe.getUtils().loadPoints(new File(peakFile));
						handler.addPoints(points);
					}else{
						System.out.println("No --peaks option is found, use All TSS mode...");
						Iterator<Point> points = nonframe.getUtils().loadTSSs();
						handler.addPoints(points);
					}
					while(handler.addingPoints()){}
					if(cluster)
						nonframe.clusterLinePanel();
					//Set the panel sizes here...
					nonframe.setLineMin(lineMin);
					nonframe.setLineMax(lineMax);
					nonframe.setLineThick(lineThick);
					nonframe.saveImages(outName);
					nonframe.savePointsToFile(outName);					
				}else if(peakFiles.size()>1){
					System.out.println("Multiple set mode...");
					MetaNonFrameMultiSet multinonframe = new MetaNonFrameMultiSet(peakFiles, gen, params, profiler, true);
					for(int x=0; x<peakFiles.size(); x++){
						String pf = peakFiles.get(x);
						Vector<Point> points = multinonframe.getUtils().loadPoints(new File(pf));
						List<MetaProfileHandler> handlers = multinonframe.getHandlers();
						handlers.get(x).addPoints(points);
						while(handlers.get(x).addingPoints()){}
					}
					multinonframe.saveImage(outName);
					multinonframe.savePointsToFile(outName);
				}
				System.out.println("Finished");
				if(profiler!=null)
					profiler.cleanup();
			}else{
				System.out.println("Initializing Meta-point frame...");
				MetaFrame frame = new MetaFrame(gen, params, profiler, normalizeProfile);
				frame.setColor(c);
				frame.setLineMax(lineMax);
				frame.setLineMin(lineMin);
				frame.setLineThick(lineThick);
				frame.startup();
				if(peakFiles.size() > 0){
					MetaProfileHandler handler = frame.getHandler();
					for(String pf : peakFiles){
						Vector<Point> points = frame.getUtils().loadPoints(new File(pf));
						handler.addPoints(points);
					}
				}
				frame.setLineMax(lineMax);
				frame.setLineMin(lineMin);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void printError(){
		System.err.println("Usage: MetaMaker (--species <organism;genome> | --g <genome chrom.sizes file> \n" +
				"--peaks <peaks file name> anchor point coordinates, can be stranded\n" +
				"--out <output root name> \n" +
				"--expt <experiment names> : can be BAM file or readdb names, add more --expt for replicates \n" +
				"--win <profile width> --bins <num bins> \n" +
				"--profiler <simplechipseq/fiveprime/chipseq/chipseqz/chipchip> \n" +
				"--back <control experiment names (only applies to chipseq)> \n" +
				"--readext <read extension> \n" +
				"--linemin <min>  --linemax <max> \n" +
				"--pbmax <per base max>\n" +
				"--color <red/green/blue/...> \n" +
				"--showcolor [flag to display all the color codes] \n" +
				"--strand <+-/> to plot only reads from one strand\n" +
				"--cluster [flag to cluster in batch mode] \n" +
				"--batch [a flag to run without displaying the window]\n" +
				"--nocolorbar [flag to turn off colorbar in batch mode]\n");
		System.err.println("\nNote: if the peaks coordinates are stranded, the line for a minus-strand anchor point will be flipped \n");
		System.exit(1);
	}
}
