package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;
import java.util.*;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.tools.motifs.WeightMatrixScanner;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;


/**
 * "Train" an HMM without access to ChIPSeq data for the transcription factor.  Instead,
 * analyze the DNAseq data and the motif to make a guess at the HMM parameters
 *
 * java edu.mit.csail.cgs.projects.dnaseq.FakeHMMTrain --species "$HS;hg19" --wm "CTCF;JASPAR 11/09 MA0139.1" --dnaseq "Crawford GM12878 DNaseSeq GM12878 against Input;statistical 1/11/11" --region "7:17m-18m" --cutoff .7 --minfold 1.5
 *
 * --region 1:3-5 (can give multiple)
 * --cutoff .7  motif cutoff as multiple of maximum LL score
 * --minfold 1.5 minimum fold change for hypersensitive region
 * --modelfile hmm.model  save hmm params to this file
 *
 * The analysis specified by --dnaseq isn't actually used; the
 * code just uses the same input alignments as that analysis.
 */

public class FakeHMMTrain extends HMMTrain {

    private double motifSampleFraction;
    public FakeHMMTrain() throws IOException, ClientException, SQLException {
        super();
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, IOException {
        super.parseArgs(args);
        motifSampleFraction = Args.parseDouble(args,"sample",.3);
    }

    protected boolean hasBinding(int position) {
        return (Math.random() < motifSampleFraction);
    }
    public static void main(String args[]) throws Exception {
        FakeHMMTrain train = new FakeHMMTrain();
        train.parseArgs(args);
        train.train();
        train.printModel();
        train.saveModel();
    }

}