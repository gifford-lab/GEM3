package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.NominalPrediction;
import weka.classifiers.evaluation.Prediction;
import weka.classifiers.meta.AdaBoostM1;
import weka.classifiers.trees.RandomForest;
import weka.core.Instances;

public class WekaUtil {
    public static BufferedReader readDataFile(String filename) {
        BufferedReader inputReader = null;
        
        try {
            inputReader = new BufferedReader(new FileReader(filename));
        } catch (FileNotFoundException ex) {
            System.err.println("File not found: " + filename);
        }
        
        return inputReader;
    }
    
    public static Evaluation simpleClassify(Classifier model, Instances trainingSet, Instances testingSet) throws Exception {
        Evaluation validation = new Evaluation(trainingSet);
        
        model.buildClassifier(trainingSet);
        
//        System.out.println(model.toString());
//        for (int i = 0; i < testingSet.numInstances(); i++) {
//        	double[] pred = model.distributionForInstance(testingSet.instance(i));
//        	System.out.print(pred[0]);
//        	System.out.print("/");
//        	System.out.print(pred[1]);
//        	System.out.print(" ");
//        }
//        System.out.println();
        
        validation.evaluateModel(model, testingSet);
        
        return validation;
    }
    
    public static double calculateAccuracy(ArrayList<Prediction> predictions) {
        double correct = 0;
        
        for (int i = 0; i < predictions.size(); i++) {
            NominalPrediction np = (NominalPrediction) predictions.get(i);
            if (np.predicted() == np.actual()) {
                correct++;
            }
        }
        
        return 100 * correct / predictions.size();
    }
    
    public static Instances[][] crossValidationSplit(Instances data, int numberOfFolds) {
        Instances[][] split = new Instances[2][numberOfFolds];
        
        for (int i = 0; i < numberOfFolds; i++) {
            split[0][i] = data.trainCV(numberOfFolds, i);
            split[1][i] = data.testCV(numberOfFolds, i);
        }
        
        return split;
    }
    
    public static void main(String[] args) throws Exception {
        // I've commented the code as best I can, at the moment.
        // Comments are denoted by "//" at the beginning of the line.
        
        BufferedReader datafile = readDataFile("/Users/yguo/data/projects/codes/tools/weka-3-8-1/data/iris.arff");
        
        Instances data = new Instances(datafile);
        data.setClassIndex(data.numAttributes() - 1);
        
        // Choose a type of validation split
        Instances[][] split = crossValidationSplit(data, 10);
        
        // Separate split into training and testing arrays
        Instances[] trainingSplits = split[0];
        Instances[] testingSplits  = split[1];
        
        // Choose a set of classifiers
        Classifier[] models = {     //new J48(),
//                                    new PART(),
//                                    new DecisionTable(),
//                                    new OneR(),
//                                    new DecisionStump(),
                                    new RandomForest(), 
                                    new AdaBoostM1()};
        
        // Run for each classifier model
        for(int j = 0; j < models.length; j++) {

            // Collect every group of predictions for current model in a FastVector
            ArrayList<Prediction> predictions = new ArrayList<Prediction>();
            
            // For each training-testing split pair, train and test the classifier
            for(int i = 0; i < trainingSplits.length; i++) {
                Evaluation validation = simpleClassify(models[j], trainingSplits[i], testingSplits[i]);
                predictions.addAll(validation.predictions());
                
                // Uncomment to see the summary for each training-testing pair.
                // System.out.println(models[j].toString());
            }
            
            // Calculate overall accuracy of current classifier on all splits
            double accuracy = calculateAccuracy(predictions);
            
            // Print current classifier's name and accuracy in a complicated, but nice-looking way.
            System.out.println(models[j].getClass().getSimpleName() + ": " + String.format("%.2f%%", accuracy) + "\n=====================");
        }
        
    }
}
