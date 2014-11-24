package edu.mit.csail.cgs.ewok.verbs.motifs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class WeightMatrixBestHitMapper implements Mapper<Region,WeightMatrixHit> {

	private WeightMatrixScorer scorer;
	
	public WeightMatrixBestHitMapper(WeightMatrixScorer s) { 
		scorer = s;
	}
	
	public WeightMatrixHit execute(Region a) {

		WeightMatrixScoreProfile prof = scorer.execute(a);
		int bestHit = prof.getMaxIndex();
		char bestStrand = prof.getHigherScoreStrand(bestHit);
		double score = prof.getHigherScore(bestHit);
		
		WeightMatrix m = prof.getMatrix();
	
		WeightMatrixHit hit = new WeightMatrixHit(a.getGenome(), a.getChrom(),
				bestHit, bestHit + m.matrix.length-1, score, bestStrand, m);
		return hit;
	}

	
}
