package edu.mit.csail.cgs.tools.utils;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.sql.SQLException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.ewok.verbs.MotifScanResultsGenerator;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;

/**
 * <code>Args</code> is a utility class for parsing command line arguments.  It can parse 
 * different types of values (eg, strings, integers, Genomes, WeightMatrixScans).
 *
 * <code>Args</code> provides some internal caching of previously
 * parsed objects based on the String[] object such that it returns the same Genome object
 * no matter how many times parseGenome is called.
 *
 * @author <a href="mailto:arolfe@mit.edu">Alex Rolfe</a>
 * @version 1.0
 */
public class Args {
    private static Map<String[], Organism> orgs = new HashMap<String[], Organism>();
    private static Map<String[], Genome> genomes = new HashMap<String[], Genome>();
    private static Map<String[], Set<String>> flags = new HashMap<String[],Set<String>>();
    private static Map<String[], Set<String>> arguments = new HashMap<String[],Set<String>>();
    private static Map<String, String> readdbName2Id = new HashMap<String, String>();
   
    /**
     * Parses all arguments. Similar to parseFlags, but with no restrictions on 
     * the argument not taking a value <br>
     * Returns all strings preceded by "<tt>--</tt>"
     * @param args The command line options of the form <tt>--foo</tt>
     * @return
     */
    public static Set<String> parseArgs(String args[]) {
    	if (arguments.containsKey(args)) {
            return arguments.get(args);
        }

        HashSet<String> output = new HashSet<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].matches("^--.*")){
                output.add(args[i].substring(2));
            }
        }
        arguments.put(args,output);
        return output;    	
    }
    
    /** parses flags.  These are command line options of the form 
     * --foo
     * followed by another option (eg, --foo --bar quux) or
     * the end of the command line. They take no value after the name of the argument
     *
     * @returns the Set of flags present in args[].  The Strings returned do not include the leading --
     */
    public static Set<String> parseFlags(String args[]) {
        if (flags.containsKey(args)) {
            return flags.get(args);
        }

        HashSet<String> output = new HashSet<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].matches("^--.*") &&
                ((i == args.length - 1) ||
                 args[i+1].matches("^--.*"))) {
                output.add(args[i].substring(2));
            }
        }
        flags.put(args,output);
        return output;
    }

    /** Parses the integer value of the argument named by <code>key</code> from the specified command line.   
     *  If no value is present, returns <code>defaultValue</code>.  If the key is present multiple times 
     *  on the command line, the first instance is returned.  
     *  Example:
     *  parseInteger(args,"foo",10); where args={"--minimum","1.3", "--foo","50"} returns 50.
     */
    public static int parseInteger(String args[], String key, int defaultValue) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                return Integer.parseInt(args[++i]);
            }
        }
        return defaultValue;
    }
    
    /**
     * Parses all the integers of the arguments that are named by <tt>key</tt>
     * and returns them as a <tt>Collection</tt> of <tt>Integers</tt> <br>
     * Example: parseIntegers(args, "foo"); 
     * where args ={"--foo", "3", "--min", "2.5", "--foo", "4"}  returns [3, 4].
     * @param args arguments of the command line
     * @param key the argument named by <tt>key</tt>
     * @return
     */
    public static Collection<Integer> parseIntegers(String args[], String key) {
    	ArrayList<Integer> output = new ArrayList<Integer>();
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                output.add(new Integer(args[++i]));
            }
        }
        return output;
    }
   
    /**
     * Parses all the doubles of the arguments that are named by <tt>key</tt>
     * and returns them as a <tt>Collection</tt> of <tt>Doubles</tt> <br>
     * Example: parseDoubles(args, "foo"); 
     * where args ={"--foo", "3.2", "--min", "2.5", "--foo", "4.3"}  returns [3.2, 4.3].
     * @param args arguments of the command line
     * @param key the argument named by <tt>key</tt>
     * @return
     */
    public static Collection<Double> parseDoubles(String args[], String key) {
    	ArrayList<Double> output = new ArrayList<Double>();
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                output.add(new Double(args[++i]));
            }
        }
        return output;
    }
    
    /** Parses a long from the specified command line.
     * @see edu.mit.csail.cgs.tools.utils.Args.parseInteger
     */
    public static long parseLong(String args[], String key, long defaultValue) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                return Long.parseLong(args[++i]);
            }
        }
        return defaultValue;
    }
    /** Parses a double from the specified command line.
     * @see edu.mit.csail.cgs.tools.utils.Args.parseInteger
     */
    public static double parseDouble(String args[], String key, double defaultValue) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                return Double.parseDouble(args[++i]);
            }
        }
        return defaultValue;
    }
    
    /** Parses a float from the specified command line.
     * @see edu.mit.csail.cgs.tools.utils.Args.parseInteger
     */
    public static float parseFloat(String args[], String key, float defaultValue) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                return Float.parseFloat(args[++i]);
            }
        }
        return defaultValue;
    }
    /** Parses a string from the specified command line.
     * @see edu.mit.csail.cgs.tools.utils.Args.parseInteger
     */
    public static String parseString(String args[], String key, String defaultValue) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                return args[++i].trim();
            }
        }
        return defaultValue;
    }
    
    /** Parses all strings of the argument by the name <tt>key<tt>.
     * @see edu.mit.csail.cgs.tools.utils.Args.parseIntegers
     */
    public static Collection<String> parseStrings(String args[], String key) {
        ArrayList<String> output = new ArrayList<String>();
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                output.add(args[++i]);
            }
        }
        return output;
    }
    /** Parses a filename from the command line.  The file name is either 
     * immediately preceded by <tt>--file</tt> or is any argument(s) that come(s)
     * after <tt>--</tt> and then follows the end of the command line. <br>  
     * For example:
     * <tt>--max 10 --min 3 -- foo.txt bar.txt baz.txt</tt>
     */
    public static List<String> parseFile(String args[]) {
        ArrayList<String> output = new ArrayList<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--file")) {
                output.add(args[++i]);
            }
            if (args[i].equals("--")) {
                for (int j = i + 1; j < args.length; j++) {
                    output.add(args[j]);
                }
                break;
            }
        }
        return output;
    }
    
    /** Parses a list of files of the argument by the name <tt>key</tt> 
     * from the command line and returns file handles
     * 
     */
    public static List<File> parseFileHandles(String args[], String key) {
    	if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
    	ArrayList<File> output = new ArrayList<File>();
    	for (int i = 0; i < args.length; i++) {
    		if (args[i].equals(key)) {
            	output.add(new File(args[++i]));
            }
        }
    	return(output);
    }
    /** This method returns a list containing all values preceded by the specified key.
     * for example:
     *   <tt>--input foo.txt --input bar.txt --input baz.txt</tt>
     * would return the list "foo.txt","bar.txt","baz.txt" if key="input"
     */
    public static List<String> parseList(String args[], String key) {
        if (!key.matches("^\\-\\-.*")) {
            key = "--" + key;
        }
        ArrayList<String> output = new ArrayList<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                output.add(args[++i]);
            }
        }
        return output;
    }

    /** parses <tt>--species "Mus musculus;mm8"</tt> into a Species and Genome
     *  @see edu.mit.csail.cgs.datasets.species.Organism
     *  @see edu.mit.csail.cgs.datasets.species.Genome
     */
    public static Pair<Organism,Genome> parseGenome(String args[]) throws NotFoundException {
        if (orgs.containsKey(args) && genomes.containsKey(args)) {
            return new Pair<Organism,Genome>(orgs.get(args),
                                             genomes.get(args));
        }

        String speciesname = null, genomename = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                String[] pieces = args[++i].split(";");
                speciesname = pieces[0];
                genomename = pieces[1];
            }
        }
        if(speciesname==null && genomename==null){
        	return null;
        }
        Organism org = new Organism(speciesname);
        Genome genome = org.getGenome(genomename);
        orgs.put(args,org);
        genomes.put(args,genome);
        return new Pair<Organism,Genome>(org,genome);
    }
    
    /** returns a list of <tt>ExptNameVersions</tt> from the <tt>--expt</tt> parameters. <br>
     *  <tt>--expt</tt> takes 2 or 3 semicolon separated values: name, version, and an optional replicate 
     */

    public static List<ExptNameVersion> parseENV(String args[]) {
        return parseENV(args,"expt");
    }
    
    /**
     * returns a list of <tt>ExptNameVersions</tt> from the parameters named by
     * <tt>argname</tt> <br>
     * <tt>argname</tt> takes 2 or 3 semicolon separated values: name, version, and an optional replicate
     * @see edu.mit.csail.cgs.datasets.chipchip.ExptNameVersion
     */
    public static List<ExptNameVersion> parseENV(String args[], String argname) {
        argname = "--" + argname;
        List<ExptNameVersion> output = new ArrayList<ExptNameVersion>();
        String exptname = null, exptversion = null, exptrep = null;        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(argname)) {
                String[] pieces = args[++i].split(";");
                exptname = pieces[0];
                exptversion = pieces[1];
                if (pieces.length >= 3) {
                    exptrep = pieces[2];
                }
                output.add(new ExptNameVersion(exptname,exptversion,exptrep));
            }
        }
        return output;
    }
    /**
     * parses ChipSeqLocators from the <tt>--chipseq</tt> parameters.  Takes
     * either "name;alignment" or "name;replicate;alignment"
     */
    public static List<ChipSeqLocator> parseChipSeq(String args[]) {
        return parseChipSeq(args,"chipseq");
    }

    /**
     * parses ChipSeqLocators from the <tt>argname</tt> parameters.  <br>Takes only "name;replicate;alignment".<br>
     * A readdb meta file is loaded from --readdb parameter to look up the alignment id of the readdb name.
     * @see edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator
     */
    public static List<ChipSeqLocator> parseChipSeq(String args[], String argname) {
    		if (readdbName2Id.isEmpty()) {
    			for (int i = 0; i < args.length; i++) {
    	            if (args[i].equals("--readdb")) {
    	            		ArrayList<String> lines = CommonUtils.readTextFile(args[i+1]);
    	            		for (String line: lines) {
    	            			String[] fs = line.trim().split("\t");
    	            			readdbName2Id.put(fs[1], fs[4]);
    	            		}
    	            }
    			}
    		}
    		
        argname = "--" + argname;
        ArrayList<ChipSeqLocator> output = new ArrayList<ChipSeqLocator>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(argname)) {
            		String readdbName = args[++i].trim();
//            		if (argname.equals("--chipseq") || argname.equals("--pairedchipseq") || argname.equals("--chiapetarc")) {
//            			String f[] = readdbName.split(";");
//            			output.add(new ChipSeqLocator(f[0], f[1], f[2]));
//            		}
//            		else 
            		if (readdbName2Id.containsKey(readdbName)) {
        				output.add(new ChipSeqLocator(readdbName2Id.get(readdbName)));
        			}
                else{
                    throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + args[i]);
                }
            }
        }
        return output;
    }

    public static Collection<ChipSeqAnalysis> parseChipSeqAnalyses(String args[], String argname) throws NotFoundException {
        Collection<String> bases = parseStrings(args,argname);
        Collection<ChipSeqAnalysis> out = new ArrayList<ChipSeqAnalysis>();
        ChipSeqLoader loader = null;
        try {
            for (String base : bases) {
                String pieces[] = base == null ? null : base.split(";");
                if (pieces != null && pieces.length != 2 ) {
                    throw new RuntimeException("Invalid string for ChipSeqAnalysis " + base);
                }
                ChipSeqAnalysis a = null;
                loader = new ChipSeqLoader(false);
                if (pieces != null) {
                    a = ChipSeqAnalysis.get(loader,pieces[0],pieces[1]);
                }
                out.add(a);
            }
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        } catch (Exception e) {
            throw new RuntimeException(e.toString(),e);
        } finally {
            if (loader != null) {
                loader.close();        
            }
            loader = null;
        }
        return out;
    }

    public static ChipSeqAnalysis parseChipSeqAnalysis(String args[], String argname) {
        String base = parseString(args,argname,null);
        String aname = parseString(args,"analysisname",null);
        String aversion = parseString(args,"analysisversion",null);
        String pieces[] = base == null ? null : base.split(";");
        if (pieces != null && pieces.length != 2 ) {
            throw new RuntimeException("Invalid string for ChipSeqAnalysis " + base);
        }
        ChipSeqAnalysis a = null;
        ChipSeqLoader loader = null;
        try {
            loader = new ChipSeqLoader(false);
            if (pieces != null) {
                try {
                    a = ChipSeqAnalysis.get(loader,pieces[0],pieces[1]);
                } catch (NotFoundException e) {
                    System.err.println("Couldn't find analysis from " + pieces[0] + " and " + pieces[1]);
                }
            }
            if (a == null && aname != null && aversion != null) {
                try {
                    a = ChipSeqAnalysis.get(loader,aname,aversion);
                } catch (NotFoundException e) {}
            }
            if (a == null) {
                throw new RuntimeException("Couldn't parse or find a ChipSeqAnalysis from " + base + " or " + aname +","+aversion);
            }            
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        } catch (Exception e) {
            throw new RuntimeException(e.toString(),e);
        } finally {
            if (loader != null) {
                loader.close();        
            }
            loader = null;
        }

        return a;
    }

    /**
     * parses ChipChipLocators from the <tt>chipchip</tt> parameters.
     * @see edu.mit.csail.cgs.datasets.locators.ChipChipLocator
     */
    public static List<ChipChipLocator> parseChipChip(Genome gen, String args[]) {
        return parseChipChip(gen, args,"chipchip");
    }
    
    /**
     * parses ChipChipLocators from the <tt>argname</tt> parameters.
     * @see edu.mit.csail.cgs.datasets.locators.ChipChipLocator
     */
    public static List<ChipChipLocator> parseChipChip(Genome gen, String args[], String argname) {
        argname = "--" + argname;
        ArrayList<ChipChipLocator> output = new ArrayList<ChipChipLocator>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(argname)) {
                String[] pieces = args[++i].split(";");
                if (pieces.length == 2) {
                    output.add(new ChipChipLocator(gen, pieces[0], pieces[1]));
                } else if (pieces.length == 3) {
                    output.add(new ChipChipLocator(gen, pieces[0], pieces[1], pieces[2]));
                } else {
                    System.err.println("Couldn't parse a ChipChipLocator from " + args[i]);
                }
            }
        }
        return output;
    }
    /** parses <tt>--scan</tt> or <tt>wmscan</tt> options into a list of WeightMatrixScans.  --scan takes 3 or 4 semicolon separated values:
     *   - matrix name
     *   - matrix version
     *   - scan name
     *   - organism for the matrix if not the same as what was specified in --species
     *   @see edu.mit.csail.cgs.datasets.motifs.WeightMatrix
     *   @see edu.mit.csail.cgs.datasets.motifs.WeightMatrixLoader
     *   @see edu.mit.csail.cgs.datasets.motifs.WeightMatrixScan
     */
    public static List<WeightMatrixScan> parseWMScans(String args[]) throws NotFoundException {
        Organism org = parseGenome(args).getFirst();
        List<WeightMatrixScan> output = new ArrayList<WeightMatrixScan>();
        WeightMatrixLoader wmloader = new WeightMatrixLoader();
        for (int i = 0; i < args.length; i++) {    
            if (args[i].equals("--scan") || args[i].equals("--wmscan")) {
                String[] pieces = args[++i].split(";");
                Organism thisorg = org;
                if (pieces.length == 4) {
                    thisorg = new Organism(pieces[3]);
                }

                WeightMatrix matrix = wmloader.query(thisorg.getDBID(),pieces[0],pieces[1]);
                WeightMatrixScan scan = WeightMatrixScan.getScanForMatrix(matrix.dbid, pieces[2]);
                output.add(scan);
            }
        }
        return output;
    }
    /** Similar to parseWMScans, but returns a MotifScanResultGenerator for each scan
     */
    public static List<MotifScanResultsGenerator> parseWMGenerators(String args[]) throws NotFoundException {
        ArrayList<MotifScanResultsGenerator> out = new ArrayList<MotifScanResultsGenerator>();
        for (WeightMatrixScan scan : parseWMScans(args)) {
            out.add(new MotifScanResultsGenerator(scan));
        }
        return out;
    }
    private static List<Pattern> makePatterns(Collection<String> strings) {
        List<Pattern> out = new ArrayList<Pattern>();
        for (String s : strings) {
            out.add(Pattern.compile(s));
        }
        return out;
    }
    private static boolean matchesAny(String s, List<Pattern> patterns) {
        for (Pattern p : patterns) {
            Matcher m = p.matcher(s);
            if (m.find()) {
                return true;
            }
        }
        return false;
    }
    public static Collection<WeightMatrix> filterMatrices(Collection<String> accepts,
                                                          Collection<String> rejects,
                                                          Collection<String> acceptvers,
                                                          Collection<String> rejectvers,
                                                          Collection<String> accepttypes,
                                                          Collection<String> rejecttypes,
                                                          Collection<WeightMatrix> matrices) {
        ArrayList<WeightMatrix> out = new ArrayList<WeightMatrix>();
        List<Pattern> acceptp = makePatterns(accepts);
        List<Pattern> rejectp = makePatterns(rejects);

        List<Pattern> acceptversp = makePatterns(acceptvers);
        List<Pattern> rejectversp = makePatterns(rejectvers);

        List<Pattern> accepttypep = makePatterns(accepttypes);
        List<Pattern> rejecttypep = makePatterns(rejecttypes);
        
        for (WeightMatrix wm : matrices) {
            if (matchesAny(wm.name, rejectp) ||
                matchesAny(wm.version, rejectversp) ||
                matchesAny(wm.type, rejecttypep)) {
                continue;
            }
            if ((matchesAny(wm.name, acceptp) ||
                 matchesAny(wm.version, acceptversp) ||
                 matchesAny(wm.type, accepttypep)) ||
                (accepts.size() == 0 && acceptvers.size() == 0 && accepttypes.size() == 0)) {
                out.add(wm);
            }
        }
        return out;
    }
    public static Collection<WeightMatrix> parseWeightMatrices(String args[]) throws NotFoundException {
        Collection<String> awm = parseStrings(args,"acceptwm");
        Collection<String> rwm = parseStrings(args,"rejectwm");

        Collection<String> awmv = parseStrings(args,"acceptwmver");
        Collection<String> rwmv = parseStrings(args,"rejectwmver");

        Collection<String> awmt = parseStrings(args,"acceptwmtype");
        Collection<String> rwmt = parseStrings(args,"rejectwmtype");

        if (awm.size() > 0 ) { System.err.println("Acceping wmnames " + awm);}
        if (awmv.size() > 0 ) { System.err.println("Acceping wmvers " + awmv);}
        if (awmt.size() > 0 ) { System.err.println("Acceping wmtypes " + awmt);}
        if (rwm.size() > 0 ) { System.err.println("Rejecting wmnames " + rwm);}
        if (rwmv.size() > 0 ) { System.err.println("Rejecting wmvers " + rwmv);}
        if (rwmt.size() > 0 ) { System.err.println("Rejecting wmtypes " + rwmt);}


        Collection<WeightMatrix> out = new ArrayList<WeightMatrix>();
        if (awm.size() > 0 || rwm.size() > 0 || awmv.size() > 0 ||
            rwmv.size() > 0 || awmt.size() > 0 || rwmt.size() > 0) {
            ArrayList<WeightMatrix> matrices = new ArrayList<WeightMatrix>();
            matrices.addAll(WeightMatrix.getAllWeightMatrices());
            out.addAll(filterMatrices(awm, rwm, awmv, rwmv, awmt, rwmt, matrices));
        }

        Collection<String> namevers = parseStrings(args,"wm");
        if (namevers.size() > 0) {
            WeightMatrixLoader loader = new WeightMatrixLoader();
            for (String nv : namevers) {
                String pieces[] = nv.split(";");                     
                for (WeightMatrix m : loader.query(pieces[0], pieces[1], null)) {
                    if (m.getName().equals(pieces[0]) && m.getVersion().equals(pieces[1])) {
                        out.add(m);
                    }
                }                
            }
            loader.close();
        } else if (awm.size() == 0 && rwm.size() == 0 && awmv.size() == 0 &&
                   rwmv.size() == 0 && awmt.size() == 0 && rwmt.size() == 0) {
            return WeightMatrix.getAllWeightMatrices();
        }
        System.err.println("parseWeightMatrices returning " + out.size());

        return out;
    }

    /** regularization computes ratios as ratio = (ip + alpha)/(wce + alpha) <br>
     * Argument should appear as: <tt>--regularize alpha</tt>
     */
    public static int parseRegularize(String args[]) {
        int alpha = 0;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--regularize")) {
                alpha = Integer.parseInt(args[++i]);
            }
        }
        return alpha;
    }
    /** returns a list of RefGeneGenerator.  If --genes is specified, uses the value as the type rather than
     *  the default type of refGene.  The type is the name of the table that the RefGeneGenerator will query.
     */
    public static List<RefGeneGenerator> parseGenes (String args[]) throws NotFoundException {
        ArrayList<RefGeneGenerator> output = new ArrayList<RefGeneGenerator>();
        boolean flipstrands = parseFlags(args).contains("flipgenestrands") || parseFlags(args).contains("flipgenestrand");
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--genes")) {
                output.add(new RefGeneGenerator(parseGenome(args).getLast(),args[++i]));
            }
        }
        if (flipstrands) {
            for (RefGeneGenerator rg : output) {
                rg.setFlipStrand(true);
            }
        }
        return output;
    }

    /** Takes a <tt>key</tt> that specifies the name of the command line option.  
     * For example, if key is <tt>quux</tt>, then looks for <tt>--quux</tt>.  
     * The value after each <tt>--quux</tt> is parsed as specifying a filename 
     * that should be opened and read.  
     * Each line is parsed as a region and those regions are returned, sorted.
     * @see edu.mit.csail.cgs.datasets.species.Genome
     * @see edu.mit.csail.cgs.datasets.general.Region
     */
    public static List<Region> readLocations(String args[], String key) throws IOException, NotFoundException {
        Genome genome = parseGenome(args).getLast();
        ArrayList<Region> output = null;
        if (!key.matches("^\\-\\-")) {
            key = "--" + key;
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                String fname = args[++i];
                output = new ArrayList<Region>();
                BufferedReader reader;
                if (fname.equals("-")) {
                    reader = new BufferedReader(new InputStreamReader(System.in));
                } else {
                    reader = new BufferedReader(new FileReader(fname));
                }
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.matches("^#.*")) {
                        continue;
                    }
                    Region r = Region.fromString(genome, line);
                    if (r != null) {
                        output.add(r);
                    }
                }
                Collections.sort(output);
                break;
            }            
        }
        return output;        
    }
    /** Parses the regions from the command line as specified by <tt>--region</tt>.
     *  The command line must also contain a <tt>--species</tt> option somewhere
     */
    public static List<Region> parseRegions(String args[]) throws NotFoundException {
        Genome genome = parseGenome(args).getLast();
        ArrayList<Region> regions = new ArrayList<Region>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--region")) {
                regions.add(Region.fromString(genome, args[++i]));
            }
        }
        return regions;
    }
    /** Parses the <tt>--region</tt> options from the command line.  If none are specified,
     *   returns regions corresponding to the chromosomes in the genome specified
     *   on the command line.
    */
    public static List<Region> parseRegionsOrDefault(String args[]) throws NotFoundException {
        Genome genome = parseGenome(args).getLast();
        ArrayList<Region> regions = new ArrayList<Region>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--region")) {
                regions.add(Region.fromString(genome, args[++i]));
            }
        }
        if (regions.size() == 0) {
            ChromRegionIterator chroms = new ChromRegionIterator(genome);                
            while (chroms.hasNext()) {
                regions.add(chroms.next());
            }                
        }
        return regions;
    }

}// end of Args class