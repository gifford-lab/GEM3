package edu.mit.csail.cgs.tools.chipseq;

import java.util.*;
import java.io.*;
import java.sql.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.datasets.chipseq.*;

public class AnalysisImporter {

    public static void main(String args[]) throws NotFoundException, SQLException, DatabaseException, IOException {
        String name = Args.parseString(args,"name",null);
        String version = Args.parseString(args,"version",null);
        String program = Args.parseString(args,"program",null);
        Genome g = Args.parseGenome(args).cdr();
        ChipSeqLoader loader = new ChipSeqLoader();
        ChipSeqAnalysis analysis = new ChipSeqAnalysis(name,version,program);

        String paramsfname = Args.parseString(args,"paramfile",null);
        if (paramsfname != null) {
            analysis.setParameters(ChipSeqLoader.readParameters(new BufferedReader(new FileReader(paramsfname))));
        }
        Set<ChipSeqAlignment> fg = new HashSet<ChipSeqAlignment>();
        Set<ChipSeqAlignment> bg = new HashSet<ChipSeqAlignment>();
        for (String s : Args.parseStrings(args,"foreground")) {
            String pieces[] = s.split(";");
            if (pieces.length == 2) {
                fg.addAll(loader.loadAlignments(new ChipSeqLocator(pieces[0],pieces[1])));
            } else if (pieces.length == 3) {
                fg.addAll(loader.loadAlignments(new ChipSeqLocator(pieces[0],pieces[1],pieces[2])));
            } else {
                System.err.println("Bad alignment spec: " + s);
            }
        }
        for (String s : Args.parseStrings(args,"background")) {
            String pieces[] = s.split(";");
            if (pieces.length == 2) {
                bg.addAll(loader.loadAlignments(new ChipSeqLocator(pieces[0],pieces[1])));
            } else if (pieces.length == 3) {
                bg.addAll(loader.loadAlignments(new ChipSeqLocator(pieces[0],pieces[1],pieces[2])));
            } else {
                System.err.println("Bad alignment spec: " + s);
            }
        }
        analysis.setInputs(fg,bg);

        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");
            analysis.addResult(new ChipSeqAnalysisResult(g,
                                                         pieces[0],
                                                         Integer.parseInt(pieces[1]),
                                                         Integer.parseInt(pieces[2]),
                                                         pieces[3].length() > 0 ? Integer.parseInt(pieces[3]) : null,
                                                         pieces[4].length() > 0 ? Double.parseDouble(pieces[4]) : null,
                                                         pieces[5].length() > 0 ? Double.parseDouble(pieces[5]) : null,
                                                         pieces[6].length() > 0 ? Double.parseDouble(pieces[6]) : null,
                                                         pieces[7].length() > 0 ? Double.parseDouble(pieces[7]) : null,
                                                         pieces[8].length() > 0 ? Double.parseDouble(pieces[8]) : null,
                                                         pieces[9].length() > 0 ? Double.parseDouble(pieces[9]) : null));
        }
        analysis.store();

    }
    
}
