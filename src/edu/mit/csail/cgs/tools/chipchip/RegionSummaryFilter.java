package edu.mit.csail.cgs.tools.chipchip;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;

import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.ewok.verbs.CoveredFilter;
import edu.mit.csail.cgs.ewok.verbs.Filter;
import edu.mit.csail.cgs.ewok.verbs.FilterIterator;
import edu.mit.csail.cgs.ewok.verbs.GeneToPromoter;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.UniqueishGeneFilter;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.stats.ListUtil;

public class RegionSummaryFilter<X extends Region> implements Filter<X,Double> {

	public static final String[][][] mark_stages = { 
		{{ "H3K27me3", "ES Stage", "es" },
			{ "H3K27me3", "ES+2d Stage, before RA", "es2" },
			{ "H3K27me3", "ES+2d Stage, 8 hours post RA", "es2" },
			{ "H3K27me3", "2+1 day", "es2" },
			{ "H3K27me3", "Olig2 Stage", "olig2" },
			{ "H3K27me3", "Hb9 Stage", "hb9" }},

			{{ "H3K4me3", "ES Stage", "es" },
				{ "H3K4me3", "ES+2d Stage, before RA", "es2" },
				{ "H3K4me3", "ES+2d Stage, 8 hours post RA", "es2" },
				{ "H3K4me3", "2+1 day", "es2" },
				{ "H3K4me3", "Olig2 Stage", "olig2" },
				{ "H3K4me3", "Hb9 Stage", "hb9" }}};

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws FileNotFoundException 
	 */
	/*
	public static void main(String[] args) throws NotFoundException, FileNotFoundException {
		Genome g = Args.parseGenome(args).cdr();
		String n = Args.parseString(args, "exptname", "");
		String v = Args.parseString(args, "version", "");
		int upstream = Args.parseInteger(args, "upstream", 500);
		int downstream = Args.parseInteger(args, "downstream", 500);
		int startBuf = Args.parseInteger(args, "startbuff", -1);
		int endBuf = Args.parseInteger(args, "endbuff", -1);
		String geneFile = Args.parseString(args, "genefile", "");
		String regionFile = Args.parseString(args, "regionfile", "");
		String outFile = Args.parseString(args, "outfile", "");
		Iterator<Gene> geneIter = new GeneFileIterator(new File(geneFile));
		geneIter = new FilterIterator<Gene,Gene>(new UniqueishGeneFilter<Gene>(startBuf,endBuf),geneIter);
		Iterator<NamedStrandedRegion> promIter = new MapperIterator<Gene,NamedStrandedRegion>(new GeneToPromoter(upstream,downstream),geneIter);
		promIter = new FilterIterator<NamedStrandedRegion,NamedStrandedRegion>(new CoveredFilter<NamedStrandedRegion>(new File(regionFile)),promIter);

		int count = 0;
		for (int i=0; i<mark_stages.length; i++) {
			for (int j=0; j<mark_stages[i].length; j++) {
				count++;
			}
		}
		List<RegionSummaryFilter<NamedStrandedRegion>> filterList = new ArrayList<RegionSummaryFilter<NamedStrandedRegion>>();
		List<String> columns = new ArrayList<String>();
		count = 0;
		for (int i=0; i<mark_stages.length; i++) {
			for (int j=0; j<mark_stages[i].length; j++) {
				columns.add(mark_stages[i][j][0] + ":" + mark_stages[i][j][1]);
				filterList.add(new RegionSummaryFilter<NamedStrandedRegion>(g, markStageToName(mark_stages[i][j][0], mark_stages[i][j][1]), v));

				count++;
			}
		}

		PrintStream out = new PrintStream(outFile);
		out.print("Gene Names");
		for (String s : columns) {
			out.print("\t" + s);
		}
		out.println();
		while (promIter.hasNext()) {
			NamedStrandedRegion tmpProm = promIter.next();
			if (tmpProm!=null) {
				out.print(tmpProm.getName());
				for (RegionSummaryFilter<NamedStrandedRegion> f : filterList) {
					out.print("\t" + f.execute(tmpProm));
				}
				out.println();
			}
		}
		out.flush();
		out.close();
	}
	*/

	Genome genome;
	String name, version;
	ChipChipLocator loc;
	List<ChipChipData> dataList;

	public RegionSummaryFilter(Genome g, String n, String v, List<String> replicates) {
		this.genome = g;
		this.name = n;
		this.version = v;
		dataList = new ArrayList<ChipChipData>();
		for (String rep : replicates) {
			loc = new ChipChipLocator(genome, name, version, rep);
			dataList.add(loc.createObject());
		}
	}

	public Double execute(X a) {
		List<Double> ratioList = new ArrayList<Double>();
		for (ChipChipData data : dataList) {
			try {
				data.window(a.getChrom(), a.getStart(), a.getEnd());
			} catch (NotFoundException e) {
				e.printStackTrace();
			}

			for (int i=0; i<data.getCount(); i++) {
				ratioList.add(data.getRatio(i, 0));
			}
		}
		return ListUtil.average(ratioList);
	}

	public static double average(double[] arr) {
		double sum = 0.0d;
		int count = 0;
		for (int i=0; i<arr.length; i++) {
			if (!Double.isNaN(arr[i])) {
				sum += arr[i];
				count++;
			}
		}
		return sum / ((double)count);
	}

	public static String markStageToName(String mark, String stage) {
		if (mark.equals("H3K79me2") && (stage.equals("ES+2d Stage") || stage.equals("ES+2d Stage, 8 hours post RA") || stage.equals("2+1 day"))) {
			return "Mm " + mark + ":HBG3:" + stage + 
			" vs WCE:HBG3:" + stage;
		} else {
			return "Mm " + mark + ":HBG3:" + stage + 
			" vs H3:HBG3:" + stage;
		}
	}

}
