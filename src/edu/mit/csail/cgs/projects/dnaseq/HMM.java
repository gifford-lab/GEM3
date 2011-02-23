package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;

public class HMM {

    public int numStates;
    public double[][] transitions;
    public double[] initialProbabilities;
    public HMMState[] states;

    public HMM (int nStates) {
        numStates = nStates;
        transitions = new double[nStates][nStates];
        initialProbabilities = new double[nStates];
        states = new HMMState[nStates];
    }
    public HMM (String fname) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = reader.readLine();
        numStates = Integer.parseInt(line);

        transitions = new double[numStates][numStates];
        initialProbabilities = new double[numStates];
        states = new HMMState[numStates];
        line = reader.readLine();
        String pieces[] = line.split("\\t");
        for (int i = 0; i < numStates; i++) {
            initialProbabilities[i] = Double.parseDouble(pieces[i]);
        }

        for (int i = 0; i < numStates; i++) {
            states[i] = HMMState.deserialize(reader.readLine());
        }

        for (int i = 0; i < numStates; i++) {
            line = reader.readLine();
            pieces = line.split("\\t");
            double sum = 0;
            for (int j = 0; j < numStates; j++) {
                transitions[i][j] = Double.parseDouble(pieces[j]);
                sum += transitions[i][j];
            }
            for (int j = 0; j < numStates; j++) {
                transitions[i][j] /= sum;
            }
        }
        reader.close();        
    }
    public void toLogProbabilities() {
        for (int i = 0; i < numStates; i++) {
            initialProbabilities[i] = Math.log(initialProbabilities[i]);
            for (int j = 0; j < numStates; j++) {
                transitions[i][j] = Math.log(transitions[i][j]);
            }
        }
    }


}