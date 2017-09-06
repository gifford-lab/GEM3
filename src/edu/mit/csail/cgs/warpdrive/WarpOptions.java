package edu.mit.csail.cgs.warpdrive;

import java.util.*;
import java.util.prefs.*;
import java.io.*;
import java.sql.SQLException;

import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.chipchip.AnalysisNameVersion;
import edu.mit.csail.cgs.datasets.chipchip.ExptNameVersion;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqExpt;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLoader;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAnalysis;
import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class WarpOptions {

	/**
	 * Constants for accessing settings
	 */
	private static final String WINDOW_WIDTH = "WINDOW_WIDTH";
	private static final String WINDOW_HEIGHT = "WINDOW_HEIGHT";
	private static final String WINDOW_IS_CENTERED = "WINDOW_IS_CENTERED";
	private static final String WINDOW_TOP_LEFT_X = "WINDOW_X";
	private static final String WINDOW_TOP_LEFT_Y = "WINDOW_Y";
	
	/**
	 * Constants for default values and limits on settings
	 */
	public static final int DEFAULT_WINDOW_WIDTH = 500;
	public static final int DEFAULT_WINDOW_HEIGHT = 650;
	public static final int MIN_WINDOW_WIDTH = 400;
	public static final int MAX_WINDOW_WIDTH = 2400;
	public static final int MIN_WINDOW_HEIGHT = 400;
	public static final int MAX_WINDOW_HEIGHT = 1600;
	public static final boolean DEFAULT_WINDOW_IS_CENTERED = true;
	public static final int DEFAULT_TOP_LEFT_X = 50;
	public static final int DEFAULT_TOP_LEFT_Y = 50;
	public static final int MIN_TOP_LEFT_X = 0;
	public static final int MAX_TOP_LEFT_X = 1000;
	public static final int MIN_TOP_LEFT_Y = 0;
	public static final int MAX_TOP_LEFT_Y = 700;
	
	// General display settings
	private int preferredWindowWidth;
	private int preferredWindowHeight;
	private boolean isWindowCentered;
	private int preferredWindowTopLeftX;
	private int preferredWindowTopLeftY;
	
	
	
    // General connection info
    public String species, genome;
    public String genomeString;			// not DB, genome file with chrom info

    // where to start the display.
    // Either use (chrom,start,stop), gene, position (which will be parsed
    // into either chrom/start/stop or gene), or a regionListFile
    public String chrom, gene, position, regionListFile;
    public int start, stop;

    // tracks to paint and their options
    public boolean hash, relative, seqletters, gccontent, pyrpurcontent, cpg, regexmatcher;
    public ArrayList<BindingScan> bindingScans;
    public ArrayList<String> genes, ncrnas, otherannots;
    public ArrayList<ExptNameVersion> agilentdata;
    public ArrayList<AnalysisNameVersion> bayesresults, agilentll, msp;
    public ArrayList<WeightMatrix> motifs;
    public ArrayList<WeightMatrixScan> motifscans;
    public ArrayList<ExptNameVersion> peakCallers;
    public ArrayList<Experiment> exprExperiments;
    public HashMap<String,String> chiapetExpts;
    public ArrayList<ChipSeqLocator> chipseqExpts;
    public ArrayList<ChipSeqLocator> pairedChipseqExpts;
    public ArrayList<ChipSeqLocator> chiapetArcs;
    public ArrayList<ChipSeqAnalysis> chipseqAnalyses;
    // filename to label mappings.  These are loaded from a file
    // and the data held statically
    public HashMap<String,String> regionTracks, regexes;
    
    // options for saving image to a file
    public boolean saveimage;
    public String filename;

    // startup-only options
    public boolean chipseqHistogramPainter = true;

    /* These constants correspond to the different input arrays.  They are used
       in WarpPaintable to store the type of input that created the paintable 
       (this must be set by the creator, currently done in RegionPanel).  This information
       is useful for removing a painter because it lets you figure out which field of
       the WarpOptions class should be modified */
    public static final int BINDINGSCAN = 1,
        GENES = 2,
        NCRNAS = 3,
        OTHERANNOTS = 4,
        AGILENTDATA = 5,
        RULERPROBES = 6,
        RULERINTERVALS = 7,
        BAYESRESULTS = 8,
        AGILENTLL = 9,
        MSP = 10,
        MOTIFS = 11,
        PEAKS = 12,
        SEQLETTERS = 13,
        GCCONTENT = 14,
        EXPRESSION = 15,
        CPG = 16,
        MOTIFSCANS = 17,
        REGEXMATCHER = 18,
        REGIONTRACKS = 19,
        PYRPURCONTENT = 20,
        CHIPSEQANALYSES = 21;
    
    
    public WarpOptions(String gname) {
        
    	this.loadOptions();
    	
        genome = gname;
        Genome g;
        try {
            g = Organism.findGenome(genome);
            species = g.getSpecies();
            
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(gname);
        }
        
        start = 10000;
        stop = 20000;        
        chrom = g.getChromList().get(0);
        gene = null;
        hash = true;
        chipseqHistogramPainter = true;
        genes = new ArrayList<String>();
        ncrnas = new ArrayList<String>();
        bindingScans = new ArrayList<BindingScan>();
        otherannots = new ArrayList<String>();
        agilentdata= new ArrayList<ExptNameVersion>();
        chiapetExpts = new HashMap<String,String>();
        chipseqExpts = new ArrayList<ChipSeqLocator>();
        pairedChipseqExpts = new ArrayList<ChipSeqLocator>();
        chiapetArcs = new ArrayList<ChipSeqLocator>();
        agilentll = new ArrayList<AnalysisNameVersion>();
        bayesresults = new ArrayList<AnalysisNameVersion>();
        msp = new ArrayList<AnalysisNameVersion>();
        motifs = new ArrayList<WeightMatrix>();
        motifscans = new ArrayList<WeightMatrixScan>();
        peakCallers = new ArrayList<ExptNameVersion>();
        exprExperiments = new ArrayList<Experiment>();
        regionTracks = new HashMap<String,String>();
        regexes = new HashMap<String,String>();
        chipseqAnalyses = new ArrayList<ChipSeqAnalysis>();
    }

    public WarpOptions() {
    	this.loadOptions();
    	
        start = -1;
        stop = -1;
        chrom = null;
        gene = null;
        hash = true;
        chipseqHistogramPainter = true;
        genes = new ArrayList<String>();
        ncrnas = new ArrayList<String>();
        bindingScans = new ArrayList<BindingScan>();
        otherannots = new ArrayList<String>();
        agilentdata= new ArrayList<ExptNameVersion>();
        chiapetExpts = new HashMap<String,String>();
        chipseqExpts = new ArrayList<ChipSeqLocator>();
        pairedChipseqExpts = new ArrayList<ChipSeqLocator>();
        chiapetArcs = new ArrayList<ChipSeqLocator>();
        agilentll = new ArrayList<AnalysisNameVersion>();
        bayesresults = new ArrayList<AnalysisNameVersion>();
        msp = new ArrayList<AnalysisNameVersion>();
        motifs = new ArrayList<WeightMatrix>();
        motifscans = new ArrayList<WeightMatrixScan>();
        peakCallers = new ArrayList<ExptNameVersion>();
        exprExperiments = new ArrayList<Experiment>();
        regionTracks = new HashMap<String,String>();
        regexes = new HashMap<String,String>();
        chipseqAnalyses = new ArrayList<ChipSeqAnalysis>();
    }

    /* adds options from this into union.  For lists, it generates the 
       union.  For other settings, this takes priority */
    public void mergeInto(WarpOptions union) throws IllegalArgumentException {
        if (union == null) {throw new NullPointerException("Must supply options to mergeInto");}
        if (species == null) {throw new NullPointerException("Tried to call mergeInto when species is null");}
        if (genome == null) {throw new NullPointerException("Tried to call mergeInto when gnome is null");}
        if (!(species.equals(union.species) &&
              genome.equals(union.genome))) {
            throw new IllegalArgumentException("Species and Genome version must match in WarpOptions.mergeInto");
        }
        if (chrom != null) {
            union.chrom = chrom;
            union.start = start;
            union.stop = stop;
        }
        if (gene != null) {
            union.gene = gene;
        }
        union.hash = hash;
        union.relative = relative;
        union.gccontent = gccontent;
        union.pyrpurcontent = pyrpurcontent;
        union.cpg = cpg;
        union.gene = gene;
        union.regexmatcher = regexmatcher;
        union.seqletters = seqletters;
        union.chipseqHistogramPainter = (chipseqHistogramPainter && union.chipseqHistogramPainter);
        mergeInto(bindingScans,union.bindingScans);
        mergeInto(genes,union.genes);
        mergeInto(ncrnas,union.ncrnas);
        mergeInto(otherannots,union.otherannots);
        mergeInto(agilentdata,union.agilentdata);
        mergeInto(chiapetExpts,union.chiapetExpts);
        mergeInto(chipseqExpts,union.chipseqExpts);
        mergeInto(pairedChipseqExpts,union.pairedChipseqExpts);
        mergeInto(chiapetArcs,union.chiapetArcs);
        mergeInto(bayesresults,union.bayesresults);
        mergeInto(agilentll,union.agilentll);
        mergeInto(msp,union.msp);
        mergeInto(motifs,union.motifs);
        mergeInto(motifscans,union.motifscans);
        mergeInto(peakCallers,union.peakCallers);
        mergeInto(exprExperiments,union.exprExperiments);
        mergeInto(regionTracks,union.regionTracks);
        mergeInto(regexes, union.regexes);
        mergeInto(chipseqAnalyses, union.chipseqAnalyses);
    }

    public void mergeInto(ArrayList source, ArrayList target) {
        for (int i = 0; i < source.size(); i++) {
            if (!target.contains(source.get(i))) {
                target.add(source.get(i));
            }
        }        
    }
    public void mergeInto(Map source, Map target) {
        for (Object k : source.keySet()) {
            if (!target.containsKey(k)) {
                target.put(k,source.get(k));
            }
        }
    }

    /* deletes options from this that are also present in other.
       Doesn't mess with chrom, start, stop, or gene 
     */
    public void differenceOf(WarpOptions other) {
        if (species != null &&
            other.species != null &&
            !other.species.equals(species)) {
            throw new IllegalArgumentException("Species must match in WarpOptions.mergeInto");
        }
        if (genome != null &&
            other.genome != null &&
            !other.genome.equals(genome)) {
            throw new IllegalArgumentException("Genome must match in WarpOptions.mergeInto");
        }
        if (other.position != null && !other.gene.equals("")) {
            position = other.position;
        }
        if (other.regionListFile != null) {
            regionListFile = other.regionListFile;
        }
        if (other.gene != null && !other.gene.equals("")) {
            gene = other.gene;
        }
        hash = hash && (!other.hash);        
        relative = relative || other.relative;
        gccontent = gccontent && (!other.gccontent);
        pyrpurcontent = pyrpurcontent && (!other.pyrpurcontent);
        cpg = cpg && (!other.cpg);
        regexmatcher = regexmatcher && (!other.regexmatcher);        
        seqletters = seqletters && (!other.seqletters);
        chipseqHistogramPainter = (chipseqHistogramPainter && other.chipseqHistogramPainter);
        differenceOf(bindingScans,other.bindingScans);
        differenceOf(genes,other.genes);
        differenceOf(ncrnas,other.ncrnas);
        differenceOf(otherannots,other.otherannots);
        differenceOf(agilentdata,other.agilentdata);
        differenceOf(chiapetExpts,other.chiapetExpts);
        differenceOf(chipseqExpts,other.chipseqExpts);
        differenceOf(pairedChipseqExpts,other.pairedChipseqExpts);
        differenceOf(chiapetArcs,other.chiapetArcs);
        differenceOf(bayesresults,other.bayesresults);
        differenceOf(agilentll,other.agilentll);
        differenceOf(msp,other.msp);
        differenceOf(motifs,other.motifs);
        differenceOf(motifscans,other.motifscans);
        differenceOf(peakCallers,other.peakCallers);
        differenceOf(exprExperiments, other.exprExperiments);
        differenceOf(regionTracks,other.regionTracks);
        differenceOf(regexes,other.regexes);
        differenceOf(chipseqAnalyses, other.chipseqAnalyses);
    }

    public void differenceOf(ArrayList removeFrom, ArrayList other) {
        removeFrom.removeAll(other);
    }
    public void differenceOf(Map removeFrom, Map other) {
        for (Object k : other.keySet()) {
            removeFrom.remove(k);
        }
    }

    public WarpOptions clone() {
        WarpOptions o = new WarpOptions();
        o.species = species;
        o.genome = genome;
        o.chrom = chrom;
        o.gene = gene;
        o.position = position;
        o.regionListFile = regionListFile;
        o.start = start;
        o.stop = stop;
        o.hash = hash;
        o.gccontent = gccontent;
        o.pyrpurcontent = pyrpurcontent;
        o.cpg = cpg;
        o.relative = relative;
        o.seqletters = seqletters;
        o.regexmatcher = regexmatcher;
        o.chipseqHistogramPainter = chipseqHistogramPainter;
        o.bindingScans = (ArrayList<BindingScan>)bindingScans.clone();
        o.genes = (ArrayList<String>) genes.clone();
        o.ncrnas = (ArrayList<String>) ncrnas.clone();
        o.otherannots = (ArrayList<String>) otherannots.clone();
        o.agilentdata = (ArrayList<ExptNameVersion>)agilentdata.clone();
        o.chiapetExpts = (HashMap<String,String>)chiapetExpts.clone();
        o.chipseqExpts = (ArrayList<ChipSeqLocator>)chipseqExpts.clone();
        o.pairedChipseqExpts = (ArrayList<ChipSeqLocator>)pairedChipseqExpts.clone();
        o.chiapetArcs = (ArrayList<ChipSeqLocator>)chiapetArcs.clone();
        o.bayesresults = (ArrayList<AnalysisNameVersion>)bayesresults.clone();
        o.agilentll = (ArrayList<AnalysisNameVersion>)agilentll.clone();
        o.msp = (ArrayList<AnalysisNameVersion>)msp.clone();
        o.motifs = (ArrayList<WeightMatrix>)motifs.clone();
        o.motifscans = (ArrayList<WeightMatrixScan>)motifscans.clone();
        o.peakCallers = (ArrayList<ExptNameVersion>)peakCallers.clone();
        o.exprExperiments = (ArrayList<Experiment>)exprExperiments.clone();
        o.regionTracks = (HashMap<String,String>)regionTracks.clone();
        o.regexes = (HashMap<String,String>)regexes.clone();
        o.chipseqAnalyses = (ArrayList<ChipSeqAnalysis>)chipseqAnalyses.clone();
        return o;
    }

    /* Fills in a WarpOptions from command line arguments */
    public static WarpOptions parseCL(String[] args) throws NotFoundException, SQLException, IOException {
    	/**
    	 * TODO add a restore defaults option and an import option so that if 
    	 * someone screws up the warp drive options that are maintained with 
    	 * java.util.prefs.Preferences reasonable values can be restored  
    	 */
        WarpOptions opts = new WarpOptions();

        try {        
            ResourceBundle res = ResourceBundle.getBundle("defaultgenome");
            opts.species = res.getString("species");
            opts.genome = res.getString("genome");
        } catch (MissingResourceException e) {
            // who cares, we're just getting defaults
        } catch (Exception e) {
            // ditto
        }
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                opts.species = args[++i];
                if (opts.species.indexOf(';') != -1) {
                    String[] pieces = opts.species.split(";");
                    opts.species = pieces[0];
                    opts.genome = pieces[1];
                }
            }  
            if (args[i].equals("--g")) {
                opts.genomeString = args[++i];
            }
            if (args[i].equals("--genome") || args[i].equals("--genomeversion")) {
                opts.genome = args[++i];
            }
            if (args[i].equals("--oldchipseq")) {
                opts.chipseqHistogramPainter = false;
                System.err.println("Will use old ChipSeq painters");
            }

        }
        ChipSeqLoader chipseqloader = null;
        WeightMatrixLoader wmloader = null;
        if (opts.genomeString==null){
            wmloader = new WeightMatrixLoader();
            chipseqloader = new ChipSeqLoader();
        }
        try {
            Genome genome = null; Organism organism = null;
            if (opts.species != null && opts.genome != null && opts.genomeString==null) {
                organism = new Organism(opts.species);
                genome = organism.getGenome(opts.genome);
            }

            System.err.println("Parsing args of length " + args.length);
            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("--chrom") || args[i].equals("--region")) {
                    if (args[i+1].matches(".+:.+\\-.+")) {
                        int colon = args[i+1].indexOf(':');
                        int dash = args[i+1].indexOf('-');
                        opts.chrom = args[i+1].substring(0, colon);
                        String startstring = args[i+1].substring(colon+1,dash);
                        String stopstring = args[i+1].substring(dash+1);
                        opts.start = stringToNum(startstring);
                        opts.stop = stringToNum(stopstring);
                        i++;
                    } else {
                        opts.chrom = args[i+1];
                    }
                }
                if (args[i].equals("--saveimage")) {
                    opts.saveimage = true;
                    opts.filename = args[++i];
                }
                if (args[i].equals("--gene")) {
                    opts.gene = args[++i];
                }
                if (args[i].equals("--hash")) { 
                    opts.hash = true;
                }
                if (args[i].equals("--gccontent")) {
                    opts.gccontent = true;
                }
                if (args[i].equals("--pyrpurcontent")) {
                    opts.pyrpurcontent = true;
                }
                if (args[i].equals("--cpg")) {
                    opts.cpg = true;
                }
                if (args[i].equals("--repeatmasked")) {
                    opts.otherannots.add("RepeatMasker");
                }
                if (args[i].equals("--repeatmasker")) {
                    opts.otherannots.add("RepeatMasker");
                }
                if (args[i].equals("--cpgislands")) {
                    opts.otherannots.add("CpGIslands");
                }
                if (args[i].equals("--relative")) {
                    opts.relative = true;
                }
                if (args[i].equals("--genes")) {
                    opts.genes.add(args[++i]);;
                }            
                if (args[i].equals("--chiapet")) {
                	String pieces[] = args[++i].split(";");
                    if (pieces.length == 1) {
                        opts.chiapetExpts.put(pieces[0],pieces[0]);
                    } else {
                        opts.chiapetExpts.put(pieces[0],pieces[1]);
                    }
                }
                if (args[i].equals("--chipseq")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 2) {
                        opts.chipseqExpts.add(new ChipSeqLocator(pieces[0], pieces[1]));
                    } else if (pieces.length >= 3) {
                        Set<String> repnames = new HashSet<String>();
                        for (int j = 1; j < pieces.length - 1; j++) {
                            repnames.add(pieces[j]);
                        }
                        opts.chipseqExpts.add(new ChipSeqLocator(pieces[0], repnames, pieces[pieces.length-1]));
                    } else {
                        System.err.println("Couldn't parse --chipseq " + args[i]);
                    }
                }
                if (args[i].equals("--pairedchipseq")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 2) {
                        opts.pairedChipseqExpts.add(new ChipSeqLocator(pieces[0], pieces[1]));
                    } else if (pieces.length == 3) {
                        opts.pairedChipseqExpts.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
                    } else {
                        System.err.println("Couldn't parse --pairedchipseq " + args[i]);
                    }
                }
                if (args[i].equals("--chiapetarc")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 2) {
                        opts.chiapetArcs.add(new ChipSeqLocator(pieces[0], pieces[1]));
                    } else if (pieces.length == 3) {
                        opts.chiapetArcs.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
                    } else {
                        System.err.println("Couldn't parse --chiapetarc " + args[i]);
                    }
                }
                if (args[i].equals("--chipseqanalysis")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 2) {
                        opts.chipseqAnalyses.add(ChipSeqAnalysis.get(chipseqloader, pieces[0], pieces[1]));
                    } else {
                        System.err.println("Couldn't parse --chipseqanalysis " + args[i]);
                    }                
                }

                if (args[i].equals("--agilent") || args[i].equals("--chipchip")) {                
                    System.err.println("Parsing AGILENT option");
                    System.err.println("args[i+1] = " + args[i+1]);
                    String pieces[] = args[++i].split(";");
                    ExptNameVersion env = null;
                    if (pieces.length == 2) {
                        env = new ExptNameVersion(pieces[0],pieces[1]);
                    } else if (pieces.length == 3) {
                        env = new ExptNameVersion(pieces[0],pieces[1], pieces[2]);
                    }
                    if (i < args.length - 2 && args[i + 1].equals("--label")) {
                        i += 2;
                        env.setLabel(args[i]);
                    }
                    opts.agilentdata.add(env);
                }
                if (args[i].equals("--ll") ||
                    args[i].equals("--mle")) {
                    String pieces[] = args[++i].split(";");
                    AnalysisNameVersion env = null;
                    env = new AnalysisNameVersion(pieces[0],pieces[1]);
                    if (i < args.length - 2 && args[i + 1].equals("--label")) {
                        i += 2;
                        env.setLabel(args[i]);
                    }
                    opts.agilentll.add(env);
                }            
                if (args[i].equals("--bayes")) {
                    String pieces[] = args[++i].split(";");
                    AnalysisNameVersion env = null;
                    env = new AnalysisNameVersion(pieces[0],pieces[1]);
                    if (i < args.length - 2 && args[i + 1].equals("--label")) {
                        i += 2;
                        env.setLabel(args[i]);
                    }
                    opts.bayesresults.add(env);
                }
                if (args[i].equals("--msp")) {
                    String pieces[] = args[++i].split(";");
                    AnalysisNameVersion env = null;
                    env = new AnalysisNameVersion(pieces[0],pieces[1]);
                    if (i < args.length - 2 && args[i + 1].equals("--label")) {
                        i += 2;
                        env.setLabel(args[i]);
                    }
                    opts.msp.add(env);             
                }
                if (args[i].equals("--sgdOther")) {
                    opts.otherannots.add("sgdOther");
                }
                if (args[i].equals("--regionList")) {
                    opts.regionListFile = args[++i];
                }
                if (args[i].equals("--otherannot")) {
                    opts.otherannots.add(args[++i]);
                }
                if (args[i].equals("--wmscan")) {
                    String[] pieces = args[++i].split(";");
                    Organism thisorg = organism;
                    if (pieces.length == 4) {
                        thisorg = new Organism(pieces[3]);
                    }                
                    WeightMatrix matrix = wmloader.query(thisorg.getDBID(),pieces[0],pieces[1]);
                    WeightMatrixScan scan = WeightMatrixScan.getScanForMatrix(matrix.dbid, pieces[2]);
                    opts.motifscans.add(scan);
                }
                if (args[i].equals("--wm")) {
                    String[] pieces = args[++i].split(";");
                    Organism thisorg = organism;
                    if (pieces.length == 3) {
                        thisorg = new Organism(pieces[2]);
                    }                
                    WeightMatrix matrix = wmloader.query(thisorg.getDBID(),pieces[0],pieces[1]);
                    opts.motifs.add(matrix);
                }
                if (args[i].equals("--regex")) {
                    opts.regexmatcher = true;
                    String[] pieces = args[++i].split(";");
                    if (pieces.length == 1) {
                        opts.regexes.put("regex" + i, pieces[0]);
                    } else if (pieces.length >= 2) {
                        opts.regexes.put(pieces[1],pieces[0]);
                    }
                }
                if (args[i].equals("--fileTrack")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 1) {
                        opts.regionTracks.put(pieces[0],pieces[0]);
                    } else {
                        opts.regionTracks.put(pieces[0],pieces[1]);
                    }
                }


            }
        } finally {
        	if (wmloader!=null)
        		wmloader.close();
        	if (chipseqloader!=null)
        		chipseqloader.close();
        }
        
        return opts;
    }


    /**
     * Import the options from an input stream
     * @param is the input stream from which to read the options
     * @throws IOException if reading from the specified output stream results in an IOException.
     * @throws InvalidPreferencesFormatException Data on input stream does not constitute a valid XML document with the mandated document type.
     */
    public void importOptions(InputStream is) throws IOException, WarpDriveException {
    	try {
			Preferences.importPreferences(is);
		} 
    	catch (InvalidPreferencesFormatException ipfex) {
    		throw new WarpDriveException(ipfex);
		}
    	this.loadOptions();
    	
    	//TODO the WarpOptionsPane needs to be updated after these options are
    	//imported
    }
    
    
    /**
     * Export the options to an output stream
     * @param os the output stream to which to write the options
     * @throws IOException if writing to the specified output stream results in an IOException.
     * @throws BackingStoreException if preference data cannot be read from backing store.
     */
    public void exportOptions(OutputStream os) throws IOException, WarpDriveException {
    	this.saveOptions();
    	try {
    		Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    		prefs.exportSubtree(os);
    	}
    	catch (BackingStoreException bsex) {
    		throw new WarpDriveException(bsex);
    	}
    }
    
    
    /**
     * Save these options to the system
     *
     */
    public void saveOptions() {
    	Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    	prefs.putInt(WINDOW_WIDTH, preferredWindowWidth);
    	prefs.putInt(WINDOW_HEIGHT, preferredWindowHeight);
    	prefs.putBoolean(WINDOW_IS_CENTERED, isWindowCentered);
    	prefs.putInt(WINDOW_TOP_LEFT_X, preferredWindowTopLeftX);
    	prefs.putInt(WINDOW_TOP_LEFT_Y, preferredWindowTopLeftY);
    }
    
    
    /**
     * Load the options from the system
     *
     */
    public void loadOptions() {
    	Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    	preferredWindowWidth = prefs.getInt(WINDOW_WIDTH, DEFAULT_WINDOW_WIDTH);
    	preferredWindowHeight = prefs.getInt(WINDOW_HEIGHT, DEFAULT_WINDOW_HEIGHT);
    	isWindowCentered = prefs.getBoolean(WINDOW_IS_CENTERED, DEFAULT_WINDOW_IS_CENTERED);
    	preferredWindowTopLeftX = prefs.getInt(WINDOW_TOP_LEFT_X, DEFAULT_TOP_LEFT_X);
    	preferredWindowTopLeftY = prefs.getInt(WINDOW_TOP_LEFT_Y, DEFAULT_TOP_LEFT_Y);
    }
    
    
    /**
     * Checks whether the specified width is within the allowable range for
     * the Warp Drive Main Frame Width
     * @param testWidth the width to test
     */
    public boolean checkPreferredWindowWidth(int testWidth) {
    	return ((testWidth <= MAX_WINDOW_WIDTH) && (testWidth >= MIN_WINDOW_WIDTH));
    }
    
    
    /**
     * Checks whether the specified height is within the allowable range for
     * the Warp Drive Main Frame Height
     * @param testHeight the height to test
     */
    public boolean checkPreferredWindowHeight(int testHeight) {
    	return ((testHeight <= MAX_WINDOW_HEIGHT) && (testHeight >= MIN_WINDOW_HEIGHT));	
    }
    
    
    /**
     * Checks whether the specified X coordinate is within the allowable range
     * for the location of the top left corner of the Warp Drive Main Frame
     * @param testX the x-coord to test
     */
    public boolean checkPreferredTopLeftX(int testX) {
    	return ((testX <= MAX_TOP_LEFT_X) && (testX >= MIN_TOP_LEFT_X));
    }
    
    
    /**
     * Checks whether the specified Y coordinate is within the allowable range
     * for the location of the top left corner of the Warp Drive Main Frame
     * @param the y-coord to test
     */
    public boolean checkPreferredTopLeftY(int testY) {
    	return ((testY <= MAX_TOP_LEFT_Y) && (testY >= MIN_TOP_LEFT_Y));
    }
 
    
    /**
     * Returns the preferred width for the Warp Drive Main Frame
     * @return the preferred width for the Warp Drive Main Frame
     */
    public int getPreferredWindowWidth() {
    	return preferredWindowWidth;
    }
    
    
    /**
     * Sets the preferred width for the Warp Drive Main Frame
     * @param newWidth the preferred width for the Warp Drive Main Frame
     */
    public boolean setPreferredWindowWidth(int newWidth) {
    	if (this.checkPreferredWindowWidth(newWidth)) {
    		this.preferredWindowWidth = newWidth;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
    
    
    /**
     * Returns the preferred height for the Warp Drive Main Frame
     * @return the preferred height for the Warp Drive Main Frame
     */
    public int getPreferredWindowHeight() {
    	return preferredWindowHeight;
    }
    
    
    /**
     * Sets the preferred height for the Warp Drive Main Frame
     * @param newHeight the preferred height for the Warp Drive Main Frame
     */
    public boolean setPreferredWindowHeight(int newHeight) {
    	if (this.checkPreferredWindowHeight(newHeight)) {
    		this.preferredWindowHeight = newHeight;
    		return true;
    	}
    	else {
    		return false;
    	}
    	
    }

    
    /**
     * Returns whether the Warp Drive Main Frame should be centered on the screen
     * @return true if the Warp Drive Main Frame should be centered on the screen
     */
    public boolean isWindowCentered() {
    	return isWindowCentered;
    }
    
    
    /**
     * Sets whether the Warp Drive Main Frame should be centered on the screen
     * @param isCentered true causes the Warp Drive Main Frame to be centered on the screen
     */
    public void setWindowCentered(boolean isCentered) {
    	this.isWindowCentered = isCentered;
    }
    
    
    /**
     * Returns the x-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     * @return the x-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     */
    public int getPreferredTopLeftX() {
    	return preferredWindowTopLeftX;
    }
    
    
    /**
     * Sets the x-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     * @param newX the x-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     */
    public boolean setPreferredTopLeftX(int newX) {
    	if (this.checkPreferredTopLeftX(newX)) {
    		this.preferredWindowTopLeftX = newX;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
    
   
    /**
     * Returns the y-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     * @return the y-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     */
    public int getPreferredTopLeftY() {
    	return preferredWindowTopLeftY;
    }
    
    
    /**
     * Sets the y-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     * @param newY the y-coord of the top left corner of the preferred location for 
     * the Warp Drive Main Frame
     */
    public boolean setPreferredTopLeftY(int newY) {
    	if (this.checkPreferredTopLeftY(newY)) {
    		this.preferredWindowTopLeftY = newY;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
 
    
    
    
    public static int stringToNum(String s) {
        if (s.matches(".*[kK]$")) {
            return 1000 * Integer.parseInt(s.substring(0,s.length()-1));
        }
        if (s.matches(".*[mM]$")) {
            return 1000000 * Integer.parseInt(s.substring(0,s.length()-1));
        } 
        return Integer.parseInt(s);
    }
}
