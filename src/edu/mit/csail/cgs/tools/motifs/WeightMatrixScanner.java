package edu.mit.csail.cgs.tools.motifs;

import java.io.*;
import java.util.*;
import java.util.Date;
import java.sql.*;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.motifs.PerBaseMotifMatch;
import edu.mit.csail.cgs.tools.utils.Args;

/** Scans a genome for a previously loaded weight matrix.  Can read either from FASTA files or from the database.  Can put results
 * into the database or print to the screen.  Can also load scan results from a file.  All matrices are converted to log-odds
 * before scanning.
 *
 * --species "$MM;mm9"
 * --scanname "90%"  [required only when storing to db]
 *
 * --wm, --acceptwm, --rejectwm, --acceptwmver, --rejectwmver, --acceptwmtype, --rejectwmtype used to specify set of matrices
 * [--print]  print results rather than storing to db
 * [--loadfile foo.txt]  load results from file rather than doing a new scan
 * [--cutoff .9] as a fraction of maximum log-odds score
 */

public class WeightMatrixScanner {
    private int inserted = 0;
    private java.sql.Connection cxn;
    private java.sql.Connection core;
    private PreparedStatement getScan, insertScan, getscanid, insertScanned, getScanned, getScannedGenome, insertScannedGenome, insertHit;
    private int madeupChromosomeID = -1;
    private Map<Integer,String> madeupChromMap = new HashMap<Integer,String>();
    private Genome genome;
    private List<String> fastafiles;
    private List<Region> regions;
    private List<WeightMatrix> matrices;
    private double cutoff;
    private String loadfile, scanname;
    private boolean print;
    private WMConsumer consumer;
    private String outfile;

    public static void main(String args[]) throws Exception {
        WeightMatrixScanner scanner = new WeightMatrixScanner();
        scanner.parseArgs(args);
        scanner.run();
        scanner.close();
    }
    public WeightMatrixScanner() throws SQLException, UnknownRoleException {
        cxn =DatabaseFactory.getConnection("annotations");
        core =DatabaseFactory.getConnection("core");
        setup();
    }
    public void close() throws SQLException {
        getScan.close();
        insertScan.close();
        getscanid.close();
        insertScanned.close();
        getScanned.close();
        getScannedGenome.close();
        insertScannedGenome.close();
        if (insertHit != null) {
            insertHit.close();
        }
        core.close();
        cxn.close();
        consumer.close();
    }
    public void parseArgs(String args[]) throws Exception {
        boolean scan = true;
        print = false;
        consumer = null;
        loadfile = null;
        genome = Args.parseGenome(args).cdr();
        cutoff = Args.parseDouble(args,"cutoff",.8);
        scanname = Args.parseString(args,"scanname",null);
        loadfile = Args.parseString(args,"loadfile",null);
        fastafiles = new ArrayList<String>();
        fastafiles.addAll(Args.parseStrings(args,"fasta"));
        print = Args.parseFlags(args).contains("print");
        regions = Args.parseRegionsOrDefault(args);
        outfile = Args.parseString(args, "outfile", "");

        if (!print) {            
            if (scanname == null) {
                System.err.println("Must supply a --scanname"); System.exit(1);
            }
        }

        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());

        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }                
        matrices = new ArrayList<WeightMatrix>();
        matrices.addAll(Args.parseWeightMatrices(args));
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }        
    }
    public void run() throws Exception {
        /* load file thing here */
        if (loadfile != null) {
            loadFile();
        } else {
            scanMatrices();
        }
        cxn.commit();
    }
    public void loadFile() {
        try {
            loadFile(genome,loadfile,(float)cutoff,consumer);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        try {
            cxn.commit();
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    } 
    public void scanMatrices() throws SQLException, FileNotFoundException {
    	DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    	if (print) {
    		if (!outfile.equals("")) {
    			consumer = new PrintConsumer(genome, null, outfile);
    		}
    	}
    	int count = 0;
        for (WeightMatrix matrix : matrices) {
            int scanid;
            float cutoffscore =(float) (matrix.getMaxScore() * cutoff);

            if (print) {
                scanid = -1;
                if (outfile.equals("")) {
                	consumer = new PrintConsumer(genome, matrix);
                } else {
                	consumer.setMatrix(matrix);
                }
            } else {
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
            }
        
            if (fastafiles.size() != 0) {
                regions.clear();
                for (String fastafile : fastafiles) {
                    scanFasta(genome,
                              matrix,
                              consumer,
                              cutoffscore,
                              fastafile,
                              regions);
                    
                }                
            } else {
                try {
                    scanFromDB(genome,
                               matrix,
                               consumer,
                               (float)cutoff,
                               regions);
                } catch (NotFoundException ex) {
                    ex.printStackTrace();
                }
            }
            if (!print) {
                storeRegionList(genome,scanid,regions);
            }
            count++;
            if (count % (matrices.size()/100)==0) {
            	System.err.println(count+" "+dfm.format(new Date()));
            }
        }
    }
    
    /* sets up SQL statements.  Must be called before anything else is done */
    public void setup() throws SQLException {
        getScan = cxn.prepareStatement("select id from weightmatrixscan where weightmatrix = ? and name = ?");            
        insertScan = cxn.prepareStatement("insert into weightmatrixscan(id,weightmatrix,name,cutoff) values " +
                                          "(weightmatrixscan_id.nextval,?,?,?)");
        getscanid = cxn.prepareStatement("select weightmatrixscan_id.currval from dual");

        insertScanned = cxn.prepareStatement("insert into wms_scanned_regions(scan,chromosome,startpos,stoppos) values(?,?,?,?)");
        getScanned = cxn.prepareStatement("select count(*) from wms_scanned_regions where scan = ? and chromosome = ? and startpos = ? and stoppos = ?");
        getScannedGenome = cxn.prepareStatement("select count(*) from wms_scanned_genomes where scan = ? and genome = ?");
        insertScannedGenome = cxn.prepareStatement("insert into wms_scanned_genomes (scan,genome) values (?,?)");
        cxn.setAutoCommit(false);
    }

    /* retrieves the DBID of the scan with the specified name and cutoff for a weight matrix */
    public int getScanID(int wmid,
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
    public boolean addScannedRegion(int scanid, int chromid, int start, int stop) throws SQLException {
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
    public void addHit(WMHit hit) throws SQLException {
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

    public void storeRegionList(Genome genome,
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

    public void loadFile (Genome genome,
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
            List<WMHit> hits = new ArrayList<WMHit>();
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
    public List<WMHit> scanFasta(Genome genome,
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
                    try {
                    	Region tmp = Region.fromString(genome, name);
                    	name = tmp.getChrom();
                    	offset = tmp.getStart();
                    	end = tmp.getEnd();
						//name = m.group(0);
						//offset = Integer.parseInt(m.group(1));
						//end = Integer.parseInt(m.group(2));
					} catch (Exception e) {
						System.err.println(name);
					}
                }
                try {
                    if (genome != null) {
                        chromid = genome.getChromID(name);
                    }
                } catch (NullPointerException e) {
                    chromid = madeupChromosomeID--;
                    madeupChromMap.put(chromid, name);                    
                }
                if (chromid >= 0) {
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
                } else {
                    for (WMHit hit : hits) {
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
            //System.err.println(ex.toString());
            //ex.printStackTrace();
        }
        return new ArrayList<WMHit>();
    }

    /* Scans a list of regions for a weight matrix using the given cutoff. */
    public void scanFromDB(Genome genome,
                           WeightMatrix matrix,
                           WMConsumer consumer,
                           float scorecutoff,
                           List<Region> regions) throws NotFoundException {
        SequenceGenerator seqgen = new SequenceGenerator();
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
                char[] bytes = seqgen.execute(new Region(region.getGenome(), region.getChrom(), rstart,rend)).toCharArray();
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
    }

    /* Scans a list of regions for a weight matrix and returns the best hit in each Region (no cutoff). */
    public static ArrayList<WMHit> scanRegionsReturnBest(Genome genome,
                                                         WeightMatrix matrix,
                                                         ArrayList<Region> regions) throws NotFoundException {
        SequenceGenerator seqgen = new SequenceGenerator();
    	ArrayList<WMHit> results = new ArrayList<WMHit>();
        for (Region region : regions) {
            int rstart = region.getStart();
            int rend = region.getEnd();
            int chromid = region.getGenome().getChromID(region.getChrom());
            int length = matrix.length();
                
            //This method is meant for getting the best hit in small regions; I am not working the region in chunks
            if (rend - rstart < length) {break;}
            char[] bytes = seqgen.execute(region).toCharArray();
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
        return results;
    }
    
    /**
     * Returns the best hit (by score) or null if no hits meet the scorecutoff
     */
    public static WMHit scanSequenceBestHit(WeightMatrix matrix,
                                            float scorecutoff,
                                            char[] sequence) {
        List<WMHit> hits = scanSequence(matrix,
                                        scorecutoff,
                                        sequence);
        if (hits.size() == 0) {
            return null;
        } else {
            WMHit output = hits.get(0);
            for (int i = 1; i < hits.size(); i++) {
                WMHit o = hits.get(i);
                if (o.getScore() > output.getScore()) {
                    output = o;
                }
            }
            return output;
        }
    }
    
    /* returns a list of WMHits.  Since this doesn't know the chromosome or
       scanid, it just fills those in with -1 for someone else to fix later.
    */
    public static List<WMHit> scanSequence(WeightMatrix matrix,
                                           float scorecutoff,
                                           char[] sequence) {
        ArrayList<WMHit> results = new ArrayList<WMHit>();
        float[] scoreleft = new float[matrix.matrix.length + 1];
        scoreleft[scoreleft.length - 1] = 0;
        for (int i = scoreleft.length - 2; i >= 0; i--) {
            double maxval = Double.NEGATIVE_INFINITY;
            if (matrix.matrix[i]['A'] > maxval) {
                maxval = matrix.matrix[i]['A'];
            }
            if (matrix.matrix[i]['C'] > maxval) {
                maxval = matrix.matrix[i]['C'];
            }
            if (matrix.matrix[i]['T'] > maxval) {
                maxval = matrix.matrix[i]['T'];
            }
            if (matrix.matrix[i]['G'] > maxval) {
                maxval = matrix.matrix[i]['G'];
            }
            scoreleft[i] = scoreleft[i+1] + (float)maxval;
        }

        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i < sequence.length - length; i++) {
            float score = 0;
            for (int j = 0; j < length && scoreleft[j] + score > scorecutoff; j++) {
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
            for (int j = 0; j < length && scoreleft[j] + score > scorecutoff; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            if (score >= scorecutoff) {
                results.add(new WMHit(-1,-1,sequence.length - i - length,sequence.length - i -1,"-",score));
            }
        }
        SequenceUtils.reverseComplement(sequence); // put the sequence back the way we found it

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


abstract class WMConsumer implements Sink<WMHit> { 
    public void consume(Collection<WMHit> hits) { 
        init();
        consume(hits.iterator());
        finish();
    }
    
    public abstract void close();
    public abstract void setMatrix(WeightMatrix m);
}

class PrintConsumer extends WMConsumer {
    
    private Genome genome;
    private WeightMatrix matrix;
    private PrintStream out;
    private boolean fileout;
    
    public PrintConsumer(Genome g, WeightMatrix m) {
        this.genome = g;
        this.matrix = m;
        this.out = System.out;
        this.fileout = false;
    }
    
    public PrintConsumer(Genome g, WeightMatrix m, String outfile) throws FileNotFoundException {
    	this.genome = g;
    	this.matrix = m;
    	this.out = new PrintStream(outfile);
    	this.fileout = true;
    }
    
    public void setMatrix(WeightMatrix m) {
    	this.matrix = m;
    }
    
    public void consume(WMHit hit) { 
        String chromname;
        try {
            chromname = genome.getChromName(hit.chromid);
        } catch (NullPointerException e) {
            if (hit.chromid < 0) {
                chromname = madeupChromMap.get(hit.chromid);
            } else {
                chromname = "unknown";
            }
        }

        out.println(matrix.toString() + "\t" + 
                           chromname + "\t" +
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
    
    public void close() {
    	if (fileout) {
    		out.flush();
    		out.close();
    	}
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
    
    public void close() {
    	
    }
    
    public void setMatrix(WeightMatrix m) {
    	
    }
}

}

