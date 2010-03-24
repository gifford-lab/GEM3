package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.util.regex.*;
import java.util.logging.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.parsing.FASTAStream;
import edu.mit.csail.cgs.ewok.verbs.io.*;
import edu.mit.csail.cgs.tools.utils.Args;


/**
 * Base ChipSeq importer.  Programs derived from GenericImporter
 * run in five phases:
 *  - read input file of read sequences and write to chipseqreads
 *  - build a map from read name to read dbid
 *  - read hits from input file(s)
 *  - sort hits by chromosome,start,stop
 *  - insert hits into chipseqhits
 *
 * Programs derived from GenericImporter will take a lot of memory to run 
 * (8GB for 9.5m reads and 30m hits) because of the map and sorting.
 * This approach is fast with our current database because inserts
 * to chipseqhits are fast if the hits are sorted by the same
 * fields as the B-tree index on that table.
 *
 * Subclasses of GenericImporter need a main() method and 
 * will parse input files.  GenericImporter provides the database interaction.
 */

public class GenericImporter {

    public static final char plus = '+', minus='-', other=' ';
    private Genome genome;
    private Organism species;
    private String paramsfname;
    private String exptname, replicate, alignmentname, factor, cells, condition;
    private int exptID, alignID, readlength;
    private PreparedStatement insertRead, insertHit, getReadID;
    private java.sql.Connection cxn;
	private Logger logger;
    private Map<String,Long> readIDs;
    private Map<Long,Integer> readCounts;
    private List<ImportHit> hits;
    private int rowsWritten;
    private boolean delayReadIDs;

    public GenericImporter() throws NotFoundException, SQLException {
        cxn = DatabaseFactory.getConnection("chipseq");
        DatabaseMetaData dbmd = cxn.getMetaData();    	
        logger = Logger.getLogger("edu.mit.csail.cgs.tools.chipseq.NovoAlignImporter");
        cxn.setAutoCommit(false);
        hits = new ArrayList<ImportHit>();
        readIDs = new HashMap<String,Long>();
        readCounts = new HashMap<Long,Integer>();
        rowsWritten = 0;
        delayReadIDs = false;
    }

    /**
     * Must be called to pull command line args:
     * --genome
     * --alignment
     * --expt
     * --readlength 
     * --factor
     * --cells
     * --condition
     */
    public void parseArgs(String args[]) throws NotFoundException {
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        species = pair.car();
        genome = pair.cdr();
        alignmentname = Args.parseString(args,"alignment",null);
        paramsfname = Args.parseString(args,"paramsfile","alignment.params");
        String exprep = Args.parseString(args,"expt",null);
        if (exprep == null) {
            throw new IllegalArgumentException("Must supply --expt 'name;replicate'");
        }
        String pieces[] = exprep.split(";");
        if (pieces.length != 2) {
            throw new IllegalArgumentException("Must supply --expt 'name;replicate'");
        }
        exptname = pieces[0];
        replicate = pieces[1];
        readlength = Args.parseInteger(args,"readlength",36);
        factor = Args.parseString(args,"factor",null);
        cells = Args.parseString(args,"cells",null);
        condition = Args.parseString(args,"condition",null);
    }

    /**
     * Must be called after parseArgs.  Creates the experiment and alignment, if necessary
     */
    public void prepare() throws SQLException, IOException {
        ChipSeqLoader chipseq = new ChipSeqLoader();
        ChipSeqExpt experiment=null;
        ChipSeqAlignment alignment;
        /* get or create the experiment */
        try {
            experiment = chipseq.loadExperiment(exptname, replicate);
        } catch (NotFoundException e) {
            logger.log(Level.INFO,"Creating a new experiment for " + exptname + "," + replicate);
            MetadataLoader core = new MetadataLoader();
            PreparedStatement insert = ChipSeqExpt.createInsert(cxn);
            insert.setString(1, exptname);
            insert.setString(2, replicate);
            insert.setInt(3, species.getDBID());
            insert.setInt(4, readlength);
            insert.setInt(5, core.getCells(cells).getDBID());
            insert.setInt(6, core.getCondition(condition).getDBID());
            insert.setInt(7, core.getFactor(factor).getDBID());
            insert.execute();
            try {
                experiment = chipseq.loadExperiment(exptname, replicate);
            } catch (NotFoundException e2) {
                /* failed again means the insert failed.  you lose */
                throw new DatabaseException("Couldn't create " + exptname + "," + replicate);
            }
        }
        
        
        alignment = null;
        try {
            alignment = chipseq.loadAlignment(experiment, alignmentname, genome);
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
            //			e.printStackTrace();
		}

        if (alignment == null) {
            try {
                logger.log(Level.INFO,"Creating a new alignment for " + exptname + "," + replicate + "," + alignmentname);
                PreparedStatement insert = ChipSeqAlignment.createInsertStatement(cxn);
                insert.setInt(1, experiment.getDBID());
                insert.setString(2, alignmentname);
                insert.setInt(3, genome.getDBID());
                insert.executeQuery();
                alignment = chipseq.loadAlignment(experiment, alignmentname, genome);
                cxn.commit();
                File f = null;
                if (paramsfname != null) {
                    f = new File(paramsfname);
                }
                if (f != null && f.exists()) {
                    System.err.println("Reading alignment parameters from " + f);
                    chipseq.addAlignmentParameters(alignment, f);

                }
            } catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
                System.err.println("Couldn't add alignment parameters");
                e.printStackTrace();
            }

        }
        if (alignment == null) {
            throw new DatabaseException("Couldn't create alignment " + alignmentname + " for " + experiment);
        }
        
        
        exptID = experiment.getDBID();
        alignID = alignment.getDBID();
        cxn.commit();
        insertRead = cxn.prepareStatement(String.format("insert /*+ append */ into chipseqreads (id, expt, name, sequence) values " +
					"(%s, %d, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqread_id"), exptID));
        getReadID = cxn.prepareStatement(Sequence.getLastSQLStatement(cxn,
                                                                      "chipseqread_id"));
        insertHit = cxn.prepareStatement(String.format("insert /*+ append */into chipseqhits " +
                                                       "(read, expt, alignment, chromosome, startpos, stoppos, strand, weight) values" +
                                                       " (?, %d, %d, ?, ?, ?, ?, ?)", exptID, alignID));
        cxn.commit();
        chipseq.close();
    }

    public Connection getConnection() {return cxn;}
    /** call when done
     */
    public void close() throws SQLException {
        cxn.commit();
        insertRead.close();
        insertHit.close();
        getReadID.close();
        cxn.close();
    }
    public void commit() throws SQLException {
        cxn.commit();
    }
    /**
     * delayReadIDs is true if the table of read name to read dbid
     * mapping will be filled in only once, by calling
     * getReadIDs().  IfdelayReadIDs is false, then the
     * table is filled in with each call to addRead.
     * A subclass that first loads a fasta file of reads
     * and then loads hits should set delayReadIDs(true)
     * to achieve better performance.
     */
    public boolean delayReadIDs() {return delayReadIDs;}
    public void delayReadIDs(boolean b) {delayReadIDs = b;}


    /** Add a read to the database
     */
    public void addRead(String readname, String readsequence) throws SQLException {
        insertRead.setString(1,readname);
        insertRead.setString(2,readsequence);
        insertRead.execute();
        if (!delayReadIDs) {
            ResultSet rs = getReadID.executeQuery();
            rs.next();
            readIDs.put(readname, rs.getLong(1));
            rs.close();
        }
        if (rowsWritten++ % 100000 == 0) {
            cxn.commit();
            logger.log(Level.INFO,"Committed through read " + readname);
        }
    }
    /**
     * Load read ids from the database
     */
    public void getReadIDs() throws SQLException {
        Statement stmt = cxn.createStatement();
        ResultSet rs = stmt.executeQuery("select name, id from chipseqreads where expt = " +
                                         exptID);
        while (rs.next()) {
            readIDs.put(rs.getString(1),
                        rs.getLong(2));
        }
        rs.close();
        stmt.close();        
    }
    /** Add a hit to the buffer of hits but doesn't add it yet to the database.
     * You must have first called addRead() on the read or called getReadIDs().
     */
    public void addHit(String readname,
                       String chrom,
                       int startpos,
                       int stoppos,
                       char strand,
                       float weight) throws NotFoundException {
        long readid = readIDs.get(readname);
        int chromid = genome.getChromID(Genome.fixChrom(chrom));
        if (!readCounts.containsKey(readid)) {
            readCounts.put(readid,0);
        }
        readCounts.put(readid,readCounts.get(readid) + 1);
        hits.add(new ImportHit(readid, chromid, startpos, stoppos, strand, weight));
        
    }
    /** Assigns a weight to a hit.  Identity function by default.  Override
     * to implement your own weighting scheme; you can call assignWeightByCount() to
     * have weight = 1 / readcount
     */
    public float assignWeight(ImportHit h) {
        return h.weight;
    }
    public float assignWeightByCount(ImportHit h) {
        return (float)(1.0 / readCounts.get(h.readid));
    }
    /**
     * Flushes the buffered hits to the database.  Calls assignWeight() on each hit
     * before it is written to the database.
     */
    public void flushHits() throws SQLException {
        Collections.sort(hits);
        for (ImportHit h : hits) {
            h.weight = assignWeight(h);
            insertHit.setLong(1,h.readid);
            insertHit.setInt(2,h.chromid);
            insertHit.setInt(3,h.startpos);
            insertHit.setInt(4,h.stoppos);
            insertHit.setString(5,Character.toString(h.strand));
            insertHit.setFloat(6,h.weight);
            insertHit.execute();
            if (rowsWritten++ % 100000 == 0) {
                logger.log(Level.INFO,String.format("COMMITING up to hit row %d, chrom %d, startpos %d, read %d",rowsWritten,h.chromid, h.startpos,h.readid));
                cxn.commit();
            }
        }
        cxn.commit();        
        hits.clear();
    }

    public Logger getLogger() {return logger;}
    public int getReadCount(String readid) {
        return readCounts.get(readIDs.get(readid));
    }
}

/** record class to hold the hits that have been buffered but not flushed
 */
class ImportHit implements Comparable<ImportHit> {

    public long readid;
    public int chromid, startpos, stoppos;
    public char strand;
    public float weight;
    
    public ImportHit(long rid, int cid, int start, int stop, char s, float w) {
        readid = rid;
        chromid = cid;
        startpos = start;
        stoppos = stop;
        strand = s;
        weight = w;
    }

    public boolean equals(Object o) {
        if (o instanceof ImportHit) {
            ImportHit other = (ImportHit)o;
            return readid == other.readid &&
                chromid == other.chromid &&
                startpos == other.startpos &&
                stoppos == other.stoppos &&
                strand == other.strand;
        } else {
            return false;
        }
    }
    public int compareTo(ImportHit o) {
        if (chromid == o.chromid) {
            if (startpos == o.startpos) {
                return stoppos - o.stoppos;
            } else {
                return startpos - o.startpos;
            }
        } else {
            return chromid - o.chromid;
        }
    }

}