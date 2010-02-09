package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.util.regex.*;
import java.util.logging.*;
import java.io.*;
import java.sql.*;
import oracle.jdbc.driver.OraclePreparedStatement;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.ewok.verbs.io.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Import results from NovoAlign, Bowtie, or Maq (or any other program with tab
 * delimited output that includes the read sequence, chromosome, startpos, and strand
 * 
 * You will need to give java a large (8GB) heap to run this program or any other program
 * extended from GenericImporter.
 *
 * NovoAlignImporter --species "$SC;SGDV1" --expt "Sc Gcn4, sigma, YPD;7/2/08" \
 *   --alignemnt "PSL > 32 matches" [--novo]
 *   --factor "Gcn4" --cells "Sigma" --condition "YPD" -- results.novoalign
 *
 * Optional params
 *  --paramsfile alignment.params
 *
 *  --maq --fasta reads.fasta
 *  --bowtie --fasta reads.fasta
 *  --eland --fasta reads.fasta
 *  --skipfasta   (reads already loaded, so only does positions)
 *
 *
 *
 */


public class NovoAlignImporter extends GenericImporter{
    private List<String> files;
    private String fastafile, skipPastRead;
    private int readNameColumn, readSequenceColumn, hitChromColumn, hitPositionColumn, hitStrandColumn;
    private boolean novo, maq, bowtie,eland, skipfasta;

	/**
	 * @throws NotFoundException
	 * @throws SQLException
	 */
    public NovoAlignImporter() throws NotFoundException, SQLException {
        super();
    }
        
    /**
     * @param args
     * @throws NotFoundException
     */
    public void parseArgs(String args[]) throws NotFoundException {
        super.parseArgs(args);
        files = Args.parseFile(args);
        fastafile = Args.parseString(args,"fasta",null);
        skipPastRead = Args.parseString(args,"skippastread",null);
        Set<String> flags = Args.parseFlags(args);
        novo = false; maq = false; bowtie = false; eland=false;
        skipfasta = flags.contains("skipfasta");
        if (flags.contains("maq")) {
            maq = true;
            readNameColumn = 0;
            readSequenceColumn = 14;
            hitChromColumn = 1;
            hitPositionColumn = 2;
            hitStrandColumn = 3;
        } else if (flags.contains("bowtie")) {
            bowtie = true;
            readNameColumn = 0;
            readSequenceColumn = 4;
            hitChromColumn = 2;
            hitPositionColumn = 3;
            hitStrandColumn = 1;
        } else if (flags.contains("eland")) {
            eland = true;
            readNameColumn = 0;
            readSequenceColumn = 1;
            hitChromColumn = 6;
            hitPositionColumn = 7;
            hitStrandColumn = 8;
        } else {
            novo = true;
            readNameColumn = 0;
            readSequenceColumn = 2;
            hitChromColumn = 7;
            hitPositionColumn = 8;
            hitStrandColumn = 9;
        }
    }
          
    public void parseAndInsert() throws SQLException, CGSException, IOException {
        boolean adding = skipPastRead == null;
        if (fastafile != null) {
            delayReadIDs(true);
            FASTAStream fastastream = new FASTAStream(new File(fastafile));
            while (fastastream.hasNext()) {
                Pair<String,String> pair = fastastream.next();
                if (adding) {
                    addRead(pair.car(), pair.cdr());
                } else if (pair.car().equals(skipPastRead)) {
                    adding = true;
                }
            }
            fastastream.close();
            getReadIDs();
        }
        if (skipfasta) {
            getReadIDs();
        }

        for (String fname : files) {
            BufferedReader reader = new BufferedReader(new FileReader(new File(fname)));
            String line=null;
            while ((line = reader.readLine()) != null) {
                try {
                    if (line.charAt(0) == '#') {continue;}                    
                    String pieces[] = line.split("\\t");
                    char s = pieces[hitStrandColumn].charAt(0);       
                    if (fastafile == null) {
                        addRead(pieces[readNameColumn], pieces[readSequenceColumn]);
                    }
                    if (novo) {
                        if (pieces.length<=6 || pieces[4].equals("QC") || pieces[4].equals("NM") || pieces[4].equals("QL")) {
                            continue;  // no hits to parse
                        }
                        s = (s == 'F' ? plus : minus);
                    }else if (eland) {
                        if (pieces.length<=8 || pieces[2].equals("QC") || pieces[2].equals("NM") || pieces[2].equals("QL")) {
                            continue;  // no hits to parse
                        }
                        s = (s == 'F' ? plus : minus);
                        pieces[hitChromColumn].replaceAll(".fa","");
                    }

                    addHit(pieces[readNameColumn],
                           pieces[hitChromColumn].replaceAll("^>",""),
                           Integer.parseInt(pieces[hitPositionColumn]) - 1,
                           Integer.parseInt(pieces[hitPositionColumn]) + pieces[readSequenceColumn].length() - 2,
                           s,
                           (float)1.0);
                } catch (Exception e) {
                    getLogger().log(Level.SEVERE,"LINE IS " + line);
                    e.printStackTrace();
                    if (e instanceof SQLException) {
                        throw (SQLException)e;
                    }
                    if (e instanceof IOException) {
                        throw (IOException)e;
                    }
                }
            }
            flushHits();
            commit();
            reader.close();
            getLogger().log(Level.INFO,"Finished with file " + fname);
        }
    }
    public float assignWeight(ImportHit h) {
        return assignWeightByCount(h);
    }
    public static void main(String args[]) throws Exception {
        NovoAlignImporter importer = new NovoAlignImporter();
        importer.parseArgs(args);
        importer.prepare();
        importer.parseAndInsert();
        importer.close();
    }
}

   
class NovoAlignHit {
    public String strand;
    public int chrom, position, length;

}