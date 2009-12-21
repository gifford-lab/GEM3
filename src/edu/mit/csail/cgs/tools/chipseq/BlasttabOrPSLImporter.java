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
 * <code>BlasttabOrPSLImporter</code> loads ChipSeq reads that have been aligned to the
 * genome; the alignments must be in either PSL or Blast8 (tabular blast output) format.
 * The FASTA file gives the sequence and name for each read.  The alignment file must
 * use the same names and present the reads in the same order.  
 *   The factor, cells, and condition for the experiment will be created in the Core
 * schema if they don't already exist. 
 *
 *
 * Usage:
 *  BlasttabOrPSLImporter --species "$SC;SGDV1" --expt "Sc Gcn4, sigma, YPD;7/2/08" \
 *   --alignemnt "PSL > 32 matches" --readlength 25
 *   --fasta "reads.fasta" --blat "reads.against_SGDv1.psl"
 *   --factor "Gcn4" --cells "Sigma" --condition "YPD"
 *
 * Optional params
 *  --minPercentIdentity .80
 *  --minHitLength 30
 *  --maxMismatchPercent 20
 *  --skipPast readname
 *  --paramsfile alignment.params
 *
 * skipPast lets you skip part of the input.  This is useful because the import commits periodically
 * and prints the last read loaded.  You can set skipPast to pick up the import at that point.
 *
 *
 *  The blat file can be either in PSL or BLAST8 format.
 */

public class BlasttabOrPSLImporter {

    private String plus = "+", minus="-", other="";
    private Genome genome;
    private Organism species;
    private String fastafname, alignmentfname, paramsfname;
    private String exptname, replicate, alignmentname, factor, cells, condition;
    private double minPercentIdentity, minHitLength, maxMismatchPercent, maxTargetGapBases;
    private String skipPast;

    private int exptID, alignID, readlength;
    private PreparedStatement insertRead, insertHit, getReadID;
    private java.sql.Connection cxn;
	private Logger logger;

    public BlasttabOrPSLImporter() throws NotFoundException, SQLException {
        cxn = DatabaseFactory.getConnection("chipseq");		
        logger = Logger.getLogger("edu.mit.csail.cgs.tools.chipseq.BlasttabOrPSLImporter");
        cxn.setAutoCommit(false);
    }

    public void parseArgs(String args[]) throws NotFoundException {
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        species = pair.car();
        genome = pair.cdr();
        fastafname = Args.parseString(args,"fasta",null);
        alignmentfname = Args.parseString(args,"blat",null);
        alignmentname = Args.parseString(args,"alignment",null);
        paramsfname = Args.parseString(args,"paramsfile",null);
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
        
        minPercentIdentity = Args.parseDouble(args,"minPercentIdentity",.60);
        minHitLength = Args.parseDouble(args,"minHitLength",0);
        maxMismatchPercent = Args.parseDouble(args,"maxMismatchPercent",10);
        maxTargetGapBases = Args.parseDouble(args,"maxTargetGapBases",4);

        skipPast = Args.parseString(args,"skipPast",null);
    }

    /**
     * Gets the experiment and alignment IDs, creating them if necessary.
     * Creates the preparedStatements that will insert the reads and hits
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
            alignment = chipseq.loadAlignment(experiment, alignmentname);
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
                alignment = chipseq.loadAlignment(experiment, alignmentname);
                cxn.commit();
                File f = null;
                if (paramsfname != null) {
                    f = new File(paramsfname);
                }
                if (f == null || (!f.exists())) {
                    File fastafile = new File(fastafname);
                    System.err.println("Fastafile is " + fastafile + " and parent is " +
                                       fastafile.getParentFile());
                    File parent = fastafile.getParentFile();
                    f = new File((parent != null ? (parent.toString() + fastafile.separator) : "") +
                                 "alignment.params");
                    System.err.println("Using default alignment parameters file name " + f);
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
        insertRead = cxn.prepareStatement(String.format(
					"insert into chipseqreads (id, expt, name, sequence) values " +
					"(%s, %d, ?, ?)", Sequence.getInsertSQL(cxn, "chipseqread_id"), exptID));
        getReadID = cxn.prepareStatement(Sequence.getLastSQLStatement(cxn,
                                                                      "chipseqread_id"));
        insertHit = cxn.prepareStatement(String.format(
					"insert into chipseqhits " +
					"(read, expt, alignment, chromosome, startpos, stoppos, strand) values" +
                    " (%s, %d, %d, ?, ?, ?, ?)", Sequence.getLastSQL(cxn, "chipseqread_id"),exptID, alignID));
        cxn.commit();
        chipseq.close();
    }

    public void parseAndInsert() throws SQLException, CGSException, IOException {
        BLATResultFileInputStream blatstream = new BLATResultFileInputStream(alignmentfname,
                                                                             genome);
        FASTAStream fastastream = new FASTAStream(new File(fastafname));
        BLATAlignmentHitRegion saved = null;
        /* keep track of how many reads have been handled so we can commit periodically */
        int readcount = 0;
        int heldfor = 0;
        boolean found = false;
        if (skipPast == null) {
            found = true;
        }

        /* This reads through the two files in parallel.  There may
           be reads with no hits, but there can't be hits with no read.
           The <code>saved</code> hit is there for when we reach
           the next read in the alignments file (ie, a hit for the next sequence
           in the fatsa file).  That hit gets saved until we move
           to it in the fasta file.
             We expect this to happen at the end of each read's section in the
           alignments file.  However, there's a warning that gets printed 
           when <code>saved</code> isn't used by the next read.  This
           happens when the next read doesn't have any hits, but it might
           also signal a problem with the input files.
        */

        while (fastastream.hasNext()) {
            Pair<String,String> pair = fastastream.next();
            String readname = pair.car();
            if (!found && skipPast != null) {
                if (skipPast.equals(readname)) {
                    found = true;
                    continue;
                } else {
                    continue;
                }
            }
            if (pair.cdr() == null || pair.cdr().length() < 10 ) {
                System.err.println("no sequence for " + pair.car() + ":" + pair.cdr());
                continue;
            }

            insertRead.setString(1, pair.car());
            insertRead.setString(2, pair.cdr());
            insertRead.execute();

            if (saved != null && saved.getName().equals(readname)) {
                insert(saved);
                saved = null;
                heldfor = 0;
            }
            if (saved == null) {
                while (blatstream.hasNext()) {
                    BLATAlignmentHitRegion hitregion;
                    try {
                        hitregion = blatstream.next();
                    } catch (Exception e) {
                        System.err.println(e.toString());
                        continue;
                    }
                    if (hitregion == null) {
                        continue;
                    }


                    if (hitregion.getName().equals(readname)) {
                        insert(hitregion);
                        saved = null;
                    } else {
                        saved = hitregion;
                        break;
                    }
                }
            } else {
                heldfor++;
                if (heldfor > 1) {
                    //                    logger.log(Level.INFO,"Holding onto hit for " + saved.getName() + " at " + readname);
                }
            }
            if (readcount++ % 10000 == 0) {
                logger.log(Level.INFO,"COMMITING up to " + readname);
                cxn.commit();
            }

        }       
        cxn.commit();
    }

    /**
     * Dumps one hit into the database.
     * This is where the filtering for minHitLength, minPercentIdentity, etc
     * takes place.  This must be called after the Read is inserted
     * because it uses chipseqread_id.currval
     */
    private void insert(BLATAlignmentHitRegion hit) throws SQLException {
        if (hit.getPercentIdentity() < minPercentIdentity) {
            return;
        }
        if (hit.getWidth() < minHitLength) {
            return;
        }
        if (hit instanceof PSLHitRegion) {
            if (((PSLHitRegion)hit).getTBaseInsert() > maxTargetGapBases) {
                return;
            }           
        }

        double mmp = (100.0 * hit.getNumMismatches()) / hit.getWidth();
        if (mmp > maxMismatchPercent) {
            return;
        }

        insertHit.setInt(1, genome.getChromID(Genome.fixChrom(hit.getChrom())));
        insertHit.setInt(2, hit.getStart());
        insertHit.setInt(3, hit.getEnd());
        char c = hit.getStrand();
        insertHit.setString(4, c == '+' ? plus : (c == '-' ? minus : other));
        try {
            insertHit.execute();
        } catch (SQLException e) {
            if (e.getErrorCode() == 1) {
                // ignore it.  duplicate row in blat output
            } else {
                throw e;
            }
        }
    }

    public static void main(String args[]) throws Exception {
        BlasttabOrPSLImporter importer = new BlasttabOrPSLImporter();
        importer.parseArgs(args);
        importer.prepare();
        importer.parseAndInsert();
    }
}
