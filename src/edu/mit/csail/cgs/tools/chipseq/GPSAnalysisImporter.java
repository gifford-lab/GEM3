package edu.mit.csail.cgs.tools.chipseq;

import java.io.IOException;
import java.sql.SQLException;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAnalysisResult;

/**
 * See AnalysisImporter docs.  Command line options are the same; the only difference
 * is that GPSAnalysisImporter parses the GPS native output format.
 */

public class GPSAnalysisImporter extends AnalysisImporter {

    /* oracle complains about underflow if we don't limit the pvalues.  the actual 
       min value is somewhere between E-100 and E-200, but I didn't bother tracking 
       it down more closely since I don't think the difference really matters
    */
    public final static double minpval = Math.pow(10,-100);

    public static void main(String args[]) throws NotFoundException, SQLException, DatabaseException, IOException {
        GPSAnalysisImporter importer = new GPSAnalysisImporter();
        importer.parseArgs(args);
        importer.run(System.in);
        importer.close();
    }
    public ChipSeqAnalysisResult parseLine(String line) {
        String pieces[] = line.split("\\t");
        if (pieces[0].equals("Position")) {
            return null;
        }
        String chrompos[] = pieces[0].split(":");
        int pos = Integer.parseInt(chrompos[1]);
        if (pieces[2].equals("NA")) {
            pieces[2] = "1";
        }

        double ip = Double.parseDouble(pieces[1]);
        double wce = Math.max(Double.parseDouble(pieces[2]), .1);
        double pval = Math.max(Math.pow(10, -1 * Double.parseDouble(pieces[4])),
                               minpval);
        return new ChipSeqAnalysisResult(getGenome(),
                                         chrompos[0],
                                         pos,
                                         pos+1,
                                         pos,
                                         ip,
                                         wce,
                                         0.0,
                                         0.0,
                                         pval,
                                         ip/wce);
    }


}