/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.Random;

public class BNSearch {

	private GraphStepper stepper;
	private BNModelScore scorer;
	public BN network;
	private Random rand;
	
	private boolean walkBack;
	
	public BNSearch(BN bn) { 
		this(bn, new SimpleGraphStepper(), new MDLGraphScore(1.0));
	}
	
	public BNSearch(BN bn, GraphStepper st, BNModelScore sc) {
		rand = new Random();
		this.network = bn;
		stepper = st;
		scorer = sc;
		network.learnCPDs();
		walkBack = true;
	}
	
	public void searchStep(int i) {
		double oldScore = network.logLikelihood() - scorer.graphScore(network);
		GraphStep step = stepper.step(network.graph);

		step.forward(network.graph);
		network.learnCPDs();
		double newScore = network.logLikelihood() - scorer.graphScore(network);
		
		double ratio = Math.exp(newScore - oldScore);
		//double ratio = newScore / oldScore;
		
		boolean accepted = true;
		
		if(ratio < 1.0) {
			double p =rand.nextDouble();
			
			if(walkBack && p <= 1.0-ratio) { 
				step.reverse(network.graph);
				network.learnCPDs();
				//System.out.println("\tRejected.");
				accepted = false;
			}
		}

		if(accepted) { 
			System.out.println(String.format("Step %d: %s", i, step.toString()));
			if(ratio < 1.0) { 
				System.out.println(String.format("\tRejection Ratio: %.4f", ratio));
			}
		}
	}
	
	public void search(int steps) { 
		for(int i = 0; i < steps; i++) {
			searchStep(i);
		}
	}
}
