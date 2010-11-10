package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import java.text.DecimalFormat;

import libsvm.*;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.*;
//import edu.mit.csail.cgs.utils.probability.Binomial;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;

public class SVMCombinatorial extends CombinatorialEnrichment {

    private double trainfrac = .3;
    private double[] trainY, testY;
    private svm_node[][] trainX, testX; 
    private double[] matrixMaxScores;
    private List<String> trainKeys, testKeys;

    private svm_parameter param;        // set by parse_command_line
    private svm_problem prob;       // set by read_problem
    private svm_model model;                                                         
    private int traini, testi;
    private PrintWriter saveCalls;

    public SVMCombinatorial() {
        super();
        trainKeys = new ArrayList<String>();
        testKeys = new ArrayList<String>();
    }
    public void parseArgs(String args[]) throws Exception {
        super.parseArgs(args);
        trainfrac = Args.parseDouble(args,"trainfrac",trainfrac);
        String calls = Args.parseString(args,"savecalls",null);
        if (calls == null) {
            saveCalls = new PrintWriter("/dev/null");
        } else {
            saveCalls = new PrintWriter(calls);
        }


    }
    private void fillsvm(Map<String, WMHit[]> hits, double val) {
        for (String s : hits.keySet()) {
            WMHit[] list = hits.get(s);
            if (traini < trainY.length &&
                (Math.random() < trainfrac || testi >= testY.length)) {
                trainY[traini] = val;
                for (int j = 0; j < list.length; j++) {
                    trainX[traini][j] = new svm_node();
                    trainX[traini][j].index = j;
                    trainX[traini][j].value = list[j] == null ? 0 : list[j].getScore()/matrixMaxScores[j];
                }
                trainKeys.add(s);
                traini++;
            } else {
                testY[testi] = val;
                for (int j = 0; j < list.length; j++) {
                    testX[testi][j] = new svm_node();
                    testX[testi][j].index = j;
                    testX[testi][j].value = list[j] == null ? 0 : list[j].getScore()/matrixMaxScores[j];
                }
                testKeys.add(s);
                testi++;
            }
        }
    }
    public void setupSVM() {
        int totalexamples = fghits.size() + bghits.size();
        int trainsize = (int)(totalexamples * trainfrac);
        int testsize = totalexamples - trainsize;
        trainY = new double[trainsize];
        testY = new double[testsize];
        trainX = new svm_node[trainsize][matrices.size()];
        testX = new svm_node[testsize][matrices.size()];
        traini = 0;
        testi = 0;
        matrixMaxScores = new double[matrices.size()];
        for (int i = 0; i < matrices.size(); i++) {
            matrixMaxScores[i] = matrices.get(i).getMaxScore();
        }

        fillsvm(fghits,1.0);
        System.err.println("Used " + traini + " of " + trainY.length + " from fg dataset for training");
        fillsvm(bghits,-1.0);
        
		param = new svm_parameter();
		param.svm_type = svm_parameter.C_SVC;
		param.kernel_type = svm_parameter.LINEAR;
		param.degree = 1;
		param.gamma = 0;	// 1/num_features
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = new int[0];
		param.weight = new double[0];

        prob = new svm_problem();
        prob.l = trainY.length;
        prob.x = trainX;
        prob.y = trainY;        
    }
    public void trainSVM() {
        model = svm.svm_train(prob,param);
    }
    public void testSVM() {
        int pospos = 0, posneg = 0, negpos = 0, negneg = 0;
        for (int i = 0; i < testY.length; i++) {
            double pred = svm.svm_predict(model, testX[i]);
            if (testY[i] > 0) {
                if (pred > 0) {
                    saveCalls.println(testKeys.get(i) + " ++");
                    pospos++;
                } else {
                    saveCalls.println(testKeys.get(i) + " +-");
                    posneg++;
                }
            } else {
                if (pred > 0) {
                    saveCalls.println(testKeys.get(i) + " -+");
                    negpos++;
                } else {
                    saveCalls.println(testKeys.get(i) + " --");
                    negneg++;
                }
            }
        }
        System.out.println(String.format("++ %d, +- %d, -+ %d, -- %d", pospos, posneg, negpos, negneg));
    }
    public void report() {
        
    }
    public static void main(String args[]) throws Exception {
        SVMCombinatorial ce = new SVMCombinatorial();
        ce.parseArgs(args);
        System.err.println("Masking and saving");
        ce.maskSequence();
        ce.saveSequences();
        System.err.println("Doing weight matrix scanning");
        ce.doScans();
        System.err.println("Translating to SVM Format");
        ce.setupSVM();
        System.err.println("Training");
        ce.trainSVM();
        System.err.println("Testing");
        ce.testSVM();
        ce.report();
    }




}