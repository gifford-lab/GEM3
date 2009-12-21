package edu.mit.csail.cgs.tools.motifs;

import java.io.*;
import java.util.*;
import java.sql.*;
import java.text.ParseException;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import edu.mit.csail.cgs.utils.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;
import edu.mit.csail.cgs.ewok.verbs.motifs.PerBaseMotifMatch;

/** Scans a genome for a previously loaded weight matrix.  Can read either from FASTA files or from the database.  Can put results
 *  into the database or print to the screen.  Can also load scan results from a file
 *
 * To load previously scanned results from a file:
 *   java edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner --wmspecies 'Saccharomyces cerevisiae' --wmname HSF1 --wmversion MacIsaac06 --cutoff 60% --scanname 'for S288C comparison' --loadfile HSF1.txt
 * To perform the scan and save results to the DB:
 *   java edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner --species 'Saccharomyces cerevisiae;sigma_broad_v1' --wmname HSF1 --wmversion MacIsaac06 --cutoff 60% --scanname 'for S288C comparison' 
 *   java edu.mit.csail.cgs.toos.motifs.WeightMatrixScanner --wmspecies 'Saccharomyces cerevisiae' --wmname HSF1 --wmversion MacIsaac06 --cutoff 60% --species "Homo sapiens;hg18" --scanname 'for S288C comparison' 
 *
 * wmspecies is the species associated with the weight matrix
 * species is the species;genome that you want to scan
 *
 * ALEX TODO: is the score cutoff always relative to frequency or log-likelihood, or does it depend on how the WM is stored?
 */

public class WeightMatrixScanner {
    static int inserted = 0;
    static java.sql.Connection cxn;
    static java.sql.Connection core;
    static PreparedStatement getScan, insertScan, getscanid, insertScanned, getScanned, getScannedGenome, insertScannedGenome, getString, insertHit;


    public static void main(String args[]) {
        String wmspecies = null;
        String wmname = null, wmversion = null, scanname = null;
        String species = null, genomeversion = null;
        String cutoff = null;
        float cutoffscore = 10000000;
        boolean scan = true;
        String loadfile = null;
        ArrayList<String> fastafiles = new ArrayList<String>();
        ArrayList<Region> regions = new ArrayList<Region>();
        boolean print = false;
        WMConsumer consumer = null;
        int scanid = -1;

        try {
            cxn =DatabaseFactory.getConnection("annotations");
            core =DatabaseFactory.getConnection("core");
            setup();
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
            System.exit(1);
        } catch (SQLException ex) {
            ex.printStackTrace();
            System.exit(1);
        }

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) { 
                species = args[++i];
                if (species.indexOf(';') != -1) {
                    String[] pieces = species.split(";");
                    species = pieces[0];
                    genomeversion = pieces[1];
                }
            }
            if (args[i].equals("--genomeversion")) {
                genomeversion = args[++i];
            }
            if (args[i].equals("--wmspecies")) {
                wmspecies = args[++i];
            }
            if (args[i].equals("--wmname")) {
                wmname = args[++i];
                if (wmname.indexOf(';') != -1) {
                    String[] pieces = wmname.split(";");
                    wmname = pieces[0];
                    wmversion = pieces[1];
                }
            }
            if (args[i].equals("--wmversion")) {
                wmversion = args[++i];
            }
            if (args[i].equals("--cutoff")) {
                cutoff = args[++i];
            }
            if (args[i].equals("--scanname")) {
                scanname = args[++i];
            }
            if (args[i].equals("--loadfile")) {
                loadfile = args[++i];
            }
            if (args[i].equals("--fasta")) {
                fastafiles.add(args[++i]);
            }
            if (args[i].equals("--print")) {
                print = true;
            }
        }
        if (species == null) {
            System.err.println("Must supply a --species"); System.exit(1);
        }
        if (wmname == null) {
            System.err.println("Must supply a --wmname"); System.exit(1);
        }
        if (wmversion == null) {
            System.err.println("Must supply a --wmversion"); System.exit(1);
        }
        if (wmspecies ==null) {
            System.err.println("Using " + species + " as --wmspecies");
            wmspecies = species;
        }        
        if (!print) {            
            if (genomeversion == null) {
                System.err.println("Must supply a --genomeversion"); System.exit(1);
            }
            if (scanname == null) {
                System.err.println("Must supply a --scanname"); System.exit(1);
            }
        }
        Organism organism = null;
        try {
            organism = new Organism(wmspecies);
        } catch (NotFoundException ex) {
            System.err.println("Couldn't find species " + species);
            System.exit(1);
        }
        Genome genome = null;
        if (genomeversion != null) {
            try {
                genome = new Genome(species,genomeversion);
            } catch (NotFoundException ex) {
                System.err.println("Couldn't find genome " + species + ", " + genomeversion);
                System.exit(1);
            }
        }
        WeightMatrix matrix = null;
        try {

            int wmid = WeightMatrix.getWeightMatrixID(organism.getDBID(),
                    wmname,
                    wmversion);
            matrix = WeightMatrix.getWeightMatrix(wmid);
        } catch (NotFoundException ex) {
            System.err.println("Can't find WM " + wmname + "," + wmversion);
            System.exit(1);
        }
        float maxscore = (float)matrix.getMaxScore();
        if (cutoff != null && cutoff.matches(".*%.*")) {  // handle percent rather than absolute cutoffs
            String c = cutoff.replaceAll("\\s*%.*","");
            System.err.println("Parsing " + c);
            System.err.println("maxscore is " + maxscore);
            cutoffscore = Float.parseFloat(c) / (float)100.0 * maxscore;
        } else if (cutoff == null) {
            cutoffscore = (float).9 * maxscore;
        } else  {
            cutoffscore = Float.parseFloat(cutoff);
        }

        if (print) {
            consumer = new PrintConsumer(genome);
            scanid = -1;
        } else {
            try {
                scanid = getScanID(matrix.dbid,scanname,cutoffscore);
                System.err.println("SCAN ID is " + scanid);
                insertHit = cxn.prepareStatement("insert into wms_hits(scan,chromosome,startpos,stoppos,strand,score) " +
                        " values (" + scanid + ",?,?,?,?,?)");
                consumer = new StoreConsumer(cxn,
                        insertHit);
                getScannedGenome.setInt(1,scanid);
                getScannedGenome.setInt(2,genome.getDBID());
                ResultSet rs = getScannedGenome.executeQuery();
                rs.next();
                if (rs.getInt(1) == 0) {
                    insertScannedGenome.setInt(1,scanid);
                    insertScannedGenome.setInt(2,genome.getDBID());
                    insertScannedGenome.execute();
                }
                rs.close();
                getScannedGenome.close();

            } catch (SQLException ex) {
                ex.printStackTrace();
                System.exit(1);
            }
        }

        // look for chrom arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--chrom")) {
                String c = args[++i].trim();
                Region region = null;
                String mp = "^\\s*([\\w\\d]+):\\s*([,\\d]+[mMkK]?)\\s*-\\s*([,\\d]+[mMkK]?)\\s*";
                Pattern factorPattern = Pattern.compile(mp);
                Matcher m = factorPattern.matcher(c);
                if(m.find() && m.groupCount() == 3) {
                    String chromStr = m.group(1);
                    String startStr = m.group(2);
                    String endStr = m.group(3);
                    if(chromStr.startsWith("chr")) { chromStr = chromStr.substring(3, chromStr.length()); }
                    startStr = startStr.replaceAll(",", "");
                    endStr = endStr.replaceAll(",", "");
                    region = new Region(genome,
                            chromStr,
                            Region.stringToNum(startStr),
                            Region.stringToNum(endStr));
                }
                if (region != null) {
                    regions.add(region);
                } else {
                    try {
                        int start = 1;
                        int end = genome.getChromLength(c);
                        regions.add(new Region(genome,
                                c,
                                start,
                                end));
                    } catch (IllegalArgumentException ex) {
                        // no such chrom
                    } 
                }
            }
        }
        // if no regions specified, then do whole genome
        // don't do this if we're scanning a fasta file (output goes to stdout)
        // or loading a file of previously scanned results
        if (regions.size() == 0 && fastafiles.size() == 0) {
            System.err.println("No --chrom specified.  Scanning the whole genome");
            Map<String,Integer> cmap = genome.getChromLengthMap();
            for (String c : cmap.keySet()) {
                regions.add(new Region(genome,
                        c,
                        1,
                        cmap.get(c)));
            }
        }

        if (!print) {
            storeRegionList(genome,scanid,regions);
        }

        /* load file thing here */
        if (loadfile != null) {
            try {
                loadFile(genome,loadfile,cutoffscore,consumer);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            try {
                cxn.commit();
            } catch (SQLException ex) {
                ex.printStackTrace();
            }
            System.exit(0);
        }

        if (fastafiles.size() != 0) {
            regions.clear();
            for (String fastafile : fastafiles) {
                System.err.println("Scanning file " + fastafile);
                scanFasta(genome,
                        matrix,
                        consumer,
                        cutoffscore,
                        fastafile,
                        regions);

            }                
            if (!print) {
                storeRegionList(genome,scanid,regions);
            }
        } else {
            try {
                scanFromDB(genome,
                        matrix,
                        consumer,
                        cutoffscore,
                        regions);
            } catch (NotFoundException ex) {
                ex.printStackTrace();
            }
        }
        try {
            cxn.commit();
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    /* sets up SQL statements.  Must be called before anything else is done */
    public static void setup() throws SQLException {
        getScan = cxn.prepareStatement("select id from weightmatrixscan where weightmatrix = ? and name = ?");            
        insertScan = cxn.prepareStatement("insert into weightmatrixscan(id,weightmatrix,name,cutoff) values " +
        "(weightmatrixscan_id.nextval,?,?,?)");
        getscanid = cxn.prepareStatement("select weightmatrixscan_id.currval from dual");

        insertScanned = cxn.prepareStatement("insert into wms_scanned_regions(scan,chromosome,startpos,stoppos) values(?,?,?,?)");
        getScanned = cxn.prepareStatement("select count(*) from wms_scanned_regions where scan = ? and chromosome = ? and startpos = ? and stoppos = ?");
        getScannedGenome = cxn.prepareStatement("select count(*) from wms_scanned_genomes where scan = ? and genome = ?");
        insertScannedGenome = cxn.prepareStatement("insert into wms_scanned_genomes (scan,genome) values (?,?)");
        getString = core.prepareStatement("select upper(substr(sequence,?,?)) from chromsequence where id = ?");
        cxn.setAutoCommit(false);
    }

    /* retrieves the DBID of the scan with the specified name and cutoff for a weight matrix */
    public static int getScanID(int wmid,
            String scanname,
            float cutoffscore) throws SQLException {
        int scanid;
        getScan.setInt(1,wmid);
        getScan.setString(2,scanname);
        ResultSet rs = getScan.executeQuery();
        if (rs.next()) {
            scanid = rs.getInt(1);
            rs.close();
        } else {
            insertScan.setInt(1,wmid);
            insertScan.setString(2,scanname);
            insertScan.setFloat(3,cutoffscore);
            insertScan.execute();
            insertScan.close();
            rs = getscanid.executeQuery();
            if (rs.next()) {
                scanid = rs.getInt(1);
            } else {
                throw new DatabaseException("Couldn't get id for wms we just inserted");
            }
            rs.close();
        }
        return scanid;
    }

    /* adds a record to the wms_scanned_regions table.
       This needs to be done for all weight matrix scans that are storing 
       results in the database */
    public static boolean addScannedRegion(int scanid, int chromid, int start, int stop) throws SQLException {
        boolean newregion = false;
        getScanned.setInt(1,scanid);
        getScanned.setInt(2,chromid);
        getScanned.setInt(3,start);
        getScanned.setInt(4,stop);
        ResultSet rs = getScanned.executeQuery();
        rs.next();
        if (rs.getInt(1) == 0) {                    
            insertScanned.setInt(1,scanid);
            insertScanned.setInt(2,chromid);
            insertScanned.setInt(3,start);
            insertScanned.setInt(4,stop);
            insertScanned.execute();
            newregion = true;
        }           
        rs.close();
        if (inserted++ % 1000 == 0) {
            cxn.commit();
        }
        return newregion;
    }

    /* stores a WM hit to the database */
    public static void addHit(WMHit hit) throws SQLException {
        insertHit.setInt(1,hit.scanid);
        insertHit.setInt(2,hit.chromid);
        insertHit.setInt(3,hit.start);
        insertHit.setInt(4,hit.end);
        insertHit.setString(5,hit.strand);
        insertHit.setFloat(6,hit.score);
        try {
            insertHit.execute();
        } catch (SQLException ex) {
            if (!ex.getMessage().trim().matches(".*ORA.00001.*")) {
                throw ex;
            }
        }
    }

    public static void storeRegionList(Genome genome,
            int scanid,
            List<Region> regions) {
        for (Region r : regions) {                
            try {
                addScannedRegion(scanid, genome.getChromID(r.getChrom()),r.getStart(),r.getEnd());
            } catch (SQLException ex) {
                ex.printStackTrace();
                System.exit(1);
            }
        }        
    }

    public static void loadFile (Genome genome,
            String fname,
            float cutoffscore,
            WMConsumer consumer) throws FileNotFoundException, IOException {
        try {
            File resultsfile = new File(fname);
            if (!resultsfile.exists()) {
                System.err.println(fname + " doesn't exist");
                System.exit(1);
            }
            BufferedReader reader = new BufferedReader(new FileReader(resultsfile));
            String line;
            ArrayList<WMHit> hits = new ArrayList<WMHit>();
            while ((line = reader.readLine()) != null) {
                String pieces[] = line.split("\t");                
                if (pieces.length >= 5) {
                    float score = Float.parseFloat(pieces[4]);
                    if (score >= cutoffscore) {
                        hits.add(new WMHit(-1,
                                genome.getChromID(pieces[0]),
                                Integer.parseInt(pieces[1]),
                                Integer.parseInt(pieces[2]),
                                pieces[3],
                                score));
                    }
                }
            }
            consumer.consume(hits);                            
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /* scans a FASTA file for a weight matrix. */
    public static List<WMHit> scanFasta(Genome genome,
            WeightMatrix matrix,
            WMConsumer consumer,
            float cutoffscore,
            String fastafile,
            List<Region> regions) {
        try {
            File file = new File(fastafile);
            FASTAStream stream = new FASTAStream(file);
            Pattern p = Pattern.compile("(.*):(\\d*)\\-(\\d*)");
            while (stream.hasNext()) {
                Pair<String,String> pair = stream.next();
                String name = pair.getFirst();
                String seq = pair.getLast();
                pair = null;
                char[] aschars = seq.toCharArray();

                int chromid = -1;
                int offset = 0, end = -1;
                if(name.startsWith("chr")) { name = name.substring(3, name.length()); }
                Matcher m = p.matcher(name);
                if (m.matches()) {
                    name = m.group(0);
                    offset = Integer.parseInt(m.group(1));
                    end = Integer.parseInt(m.group(2));
                }
                try {
                    if (genome != null) {
                        chromid = genome.getChromID(name);
                    }
                } catch (NullPointerException e) {
                    chromid = -1;
                }
                if (chromid != -1) {
                    if (end == -1) {
                        regions.add(new Region(genome,
                                               name,
                                               0,aschars.length));
                    } else {
                        regions.add(new Region(genome,
                                               name,
                                               offset,end));
                    }                    
                }

                for (int i = 0; i < aschars.length; i++) {
                    if (aschars[i] == 'a') {aschars[i]='A';}
                    if (aschars[i] == 'c') {aschars[i]='C';}
                    if (aschars[i] == 'g') {aschars[i]='G';}
                    if (aschars[i] == 't') {aschars[i]='T';}
                }
                List<WMHit> hits = scanSequence(matrix,
                        cutoffscore,
                        aschars);

                if (chromid > 0) {
                    for (WMHit hit : hits) {
                        hit.start += offset;
                        hit.chromid = chromid;
                    }
                }
                consumer.consume(hits);
            }            
            stream.close();
        } catch (FileNotFoundException ex) {
            System.err.println(ex.toString());
            ex.printStackTrace();
        } catch (Exception ex) {
            System.err.println(ex.toString());
            ex.printStackTrace();
        }
        return new ArrayList<WMHit>();
    }

    /* Scans a list of regions for a weight matrix using the given cutoff. */
    public static void scanFromDB(Genome genome,
            WeightMatrix matrix,
            WMConsumer consumer,
            float scorecutoff,
            ArrayList<Region> regions) throws NotFoundException {
        try {
            for (Region region : regions) {
                int rstart = region.getStart();
                int chunksize = 8000000;
                int chromid = region.getGenome().getChromID(region.getChrom());
                int length = matrix.length();

                ArrayList<WMHit> results = new ArrayList<WMHit>();
                // work over the target region in pieces
                while (rstart < region.getEnd()) {
                    int rend = rstart + chunksize;
                    if (rend > region.getEnd()) {
                        rend = region.getEnd();
                    }
                    System.err.println("Working on " + rstart + " to " + rend);
                    if (rend - rstart < length) {break;}
                    getString.setInt(1,rstart);
                    getString.setInt(2,rend - rstart + 1);
                    getString.setInt(3,chromid);
                    ResultSet rs = getString.executeQuery();            
                    if (!rs.next()) {
                        throw new DatabaseException("Couldn't get desired sequence for " + chromid + " starting at " +
                                rstart + " of length " + (rend - rstart + 1));
                    }
                    char[] bytes = rs.getString(1).toCharArray();
                    rs.close();
                    for (WMHit hit : scanSequence(matrix,
                            scorecutoff,
                            bytes)) {
                        hit.chromid = chromid;
                        hit.start += rstart;
                        hit.end += rstart;
                        results.add(hit);
                    }
                    rstart = rend - length + 1;
                }
                consumer.consume(results);
            }
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }

    /* Scans a list of regions for a weight matrix and returns the best hit in each Region (no cutoff). */
    public static ArrayList<WMHit> scanRegionsReturnBest(Genome genome,
            WeightMatrix matrix,
            ArrayList<Region> regions) throws NotFoundException {
    	ArrayList<WMHit> results = new ArrayList<WMHit>();
        try {
            for (Region region : regions) {
                int rstart = region.getStart();
                int rend = region.getEnd();
                int chromid = region.getGenome().getChromID(region.getChrom());
                int length = matrix.length();
                
                //This method is meant for getting the best hit in small regions; I am not working the region in chunks
                if (rend - rstart < length) {break;}
                getString.setInt(1,rstart);
                getString.setInt(2,rend - rstart + 1);
                getString.setInt(3,chromid);
                ResultSet rs = getString.executeQuery();            
                if (!rs.next()) {
                    throw new DatabaseException("Couldn't get desired sequence for " + chromid + " starting at " +
                            rstart + " of length " + (rend - rstart + 1));
                }
                char[] bytes = rs.getString(1).toCharArray();
                rs.close();
                PerBaseMotifMatch matcher = new PerBaseMotifMatch(matrix);
                Double [] scoreScape =matcher.execute(bytes);
                char [] hitStrands = matcher.getHitStrands();
                int maxI=-1; Double maxScore=-10000000.0;
                for(int i=0; i<scoreScape.length; i++){
                	if(scoreScape[i]>maxScore){maxScore=scoreScape[i]; maxI=i;}
                	//I'm introducing a slight bias towards the center of the sequence for my own selfish reasons
                	if(scoreScape[i]==maxScore && Math.abs(i-(region.getWidth()/2))<Math.abs(maxI-(region.getWidth()/2))){maxScore=scoreScape[i]; maxI=i;}
                }
                WMHit hit =new WMHit(-1, chromid, maxI+rstart, maxI+rstart+matrix.length(), String.format("%c", hitStrands[maxI]), maxScore.floatValue());
                results.add(hit);                
            }
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }return(results);
    }
    
    /* returns a list of WMHits.  Since this doesn't know the chromosome or
       scanid, it just fills those in with -1 for someone else to fix later.
     */
    public static List<WMHit> scanSequence(WeightMatrix matrix,
            float scorecutoff,
            char[] sequence) {
        ArrayList<WMHit> results = new ArrayList<WMHit>();
        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i < sequence.length - length; i++) {
            float score = 0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            // save the hit
            if (score >= scorecutoff) {
                results.add(new WMHit(-1,-1,i,i + length - 1,"+",score));
            }
        }
        /* now reverse complement the sequence and scan that.
           basically the same as before but with slightly different
           computation of the hit coordinates */
        SequenceUtils.reverseComplement(sequence);
        for (int i = 0; i < sequence.length - length; i++) {
            float score = 0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            if (score >= scorecutoff) {
                results.add(new WMHit(-1,-1,sequence.length - i - length,sequence.length - i -1,"-",score));
            }
        }
        return results;
    }

    /* returns an array of scores that describe how well the matrix
       matches each base of the input sequence.  Indices in the output
       array correspond ot the score of the match *starting* at that
       position.  This means that the last n entries are meaningless 
       (where n = matrix length) */
    public static Double[] scanSequence(WeightMatrix matrix,
            char[] sequence) {
        return (new PerBaseMotifMatch(matrix)).execute(sequence);
    }

}

abstract class WMConsumer implements Sink<WMHit> { 
    public void consume(Collection<WMHit> hits) { 
        init();
        consume(hits.iterator());
        finish();
    }
}

class PrintConsumer extends WMConsumer {
    
    private Genome genome;
    
    public PrintConsumer(Genome g) {
        this.genome = g;
    }
    
    public void consume(WMHit hit) { 
        String chromname;
        try {
            chromname = genome.getChromName(hit.chromid);
        } catch (NullPointerException e) {
            chromname = "unknown";
        }

        System.out.println(chromname + "\t" +
                hit.start + "\t" + 
                hit.end + "\t" + 
                hit.score + "\t" +
                hit.strand);        
    }
    
    public void consume(Iterator<WMHit> hits) { 
        while(hits.hasNext()) { consume(hits.next()); }
    }

    public void finish() {
    }

    public void init() {
    }
}

class StoreConsumer extends WMConsumer {
    
    private int done;
    private java.sql.Connection cxn;
    private PreparedStatement insertHit;

    public StoreConsumer(java.sql.Connection c,
            PreparedStatement i) {
        cxn = c;
        insertHit = i;
        done = 0;
    }
    
    public void consume(WMHit hit) { 
        if (hit.chromid < 1) {
            throw new RuntimeException("chromid == -1");
        }
        try {
            insertHit.setInt(1,hit.chromid);
            insertHit.setInt(2,hit.start);
            insertHit.setInt(3,hit.end);
            insertHit.setString(4,hit.strand);
            insertHit.setFloat(5,hit.score);
            insertHit.execute();
            if (done++ % 1000 == 0) {
                cxn.commit();
            }
        } catch (SQLException ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

    public void consume(Iterator<WMHit> hits) {        
        while(hits.hasNext()) { consume(hits.next()); }
    }

    public void finish() {
        try {
            cxn.commit();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    public void init() {
        done = 0;
    }
}
