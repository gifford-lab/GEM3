package edu.mit.csail.cgs.deepseq.features;

import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KmerGroup;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class ComponentFeature extends Feature  implements Comparable<ComponentFeature>{
	// these should be set at the beginning of BindingMixure
	private static ArrayList<String> conditionNames;
	private static int numConditions=0;
	// these are set later
	private static double non_specific_ratio[];
	private static int sortingCondition=0;	// the condition to compare p-values
	
	protected Point position;
	private char strand='*';
	// isJointEvent first set using mixing prob, but later updated 
	// using inter-event distance (500) only considering significant events
	protected boolean isJointEvent = false;
	protected double mfold;
	protected double[] conditionBeta;
	protected boolean[] condSignificance;
	protected double[] condSumResponsibility;
	protected double totalSumResponsibility=0;
//	protected double logKL_plus[];
//	protected double logKL_minus[];
	protected double ipCtrl_logKL_plus[];
	protected double ipCtrl_logKL_minus[];
	protected double shapeDeviation[];
	protected double unScaledControlCounts[];
	protected double expectedCounts[];
	protected double p_values[];
	protected double p_values_wo_ctrl[];
	protected double q_value_log10[];
	protected Point EM_position;		//  EM result
	protected double alpha;
	protected double noiseFraction;
	public char getStrand(){return strand;}
	public void setStrand(char s){strand=s;}
	public double getNoiseFraction() {
		return noiseFraction;
	}
	public void setNoiseFraction(double noiseFraction) {	// store the noise fraction of the region
		this.noiseFraction = noiseFraction;
	}
	
	protected KmerGroup kmerGroup;
	public KmerGroup getKmerGroup() { return kmerGroup; }
	protected char kmerStrand;
	public char getKmerStrand(){return kmerStrand;}
	protected double enrichedKmerHGPLog10=99;
	public double getEnrichedKmerHGPLog10() {
		return enrichedKmerHGPLog10;
	}
	public void setEnrichedKmerHGPLog10(double enrichedKmerHGP) {
		this.enrichedKmerHGPLog10 = enrichedKmerHGP;
	}
	private String boundSequence;						// the aligned sequence string flanking kmer underlying this event position
	public String getBoundSequence(){return boundSequence;}
	public void setBoundSequence(String boundSequence){this.boundSequence = boundSequence;}
	public void flipBoundSequence(){boundSequence = SequenceUtils.reverseComplement(boundSequence);}
	
	public ComponentFeature(BindingComponent b){
		super(null);
		position = b.getLocation();
		strand = b.getStrand();
		coords = position.expand(1);
		if(b.getMixProb()!=1)
			isJointEvent = true;
		
		condSumResponsibility = new double[numConditions];
		for(int c = 0; c < numConditions; c++) { condSumResponsibility[c] = b.getSumResponsibility(c); }
		
		totalSumResponsibility = 0;
		for(int c = 0; c < numConditions; c++) { totalSumResponsibility  += condSumResponsibility[c]; }	
		
		conditionBeta = new double[numConditions];
		for(int c=0; c<numConditions; c++){ conditionBeta[c]=b.getConditionBeta(c); }
		
		condSignificance = new boolean[numConditions];
		p_values = new double[numConditions];
		p_values_wo_ctrl = new double[numConditions];
		q_value_log10 = new double[numConditions];
		EM_position = b.getEMPosition();
		noiseFraction = b.getNoiseFraction();
		alpha = b.getAlpha();
		kmerGroup = b.getKmerGroup();
		kmerStrand = b.getKmerStrand();
		ipCtrl_logKL_plus = new double[numConditions];
		ipCtrl_logKL_minus = new double[numConditions];
		boundSequence = b.getBoundSequence();
	}
	
	//Accessors 
	public double getAlpha(){return(alpha);}
	public Point getPosition() { return position;}
	public Point getEMPosition() { return EM_position;}
	public boolean isJointEvent() {return isJointEvent;}
	public void setJointEvent(boolean isJointEvent) {this.isJointEvent = isJointEvent;}
	public double[] getCondBetas() { return conditionBeta; }
	public boolean[] getCondSignificance() { return condSignificance; }
	public double getTotalEventStrength(){return(totalSumResponsibility);}
	public double getCondSumResponsibility(int cond){return condSumResponsibility[cond];}
	public double[] getUnscaledControlCounts(){return(unScaledControlCounts);}
	public double getExpectedCount(int cond){return(expectedCounts[cond]);}
	public double getScaledControlCounts(int cond){
		return unScaledControlCounts[cond]*non_specific_ratio[cond];
	}
	public double getEventReadCounts(int cond){
		return condSumResponsibility[cond]; 
	}
	public double get_mfold() { return mfold; }
	public void set_mfold(double mf) { mfold = mf; }
	public double getShapeDeviation(int cond) {
		return shapeDeviation[cond];
	}
	public double getAvgShapeDeviation() {
		double sum=0;
		for(int c=0; c<numConditions; c++){
			sum += shapeDeviation[c];
		}
		return sum/numConditions;
	}
	public void setShapeDeviation(double[] shapeDeviation) {
		this.shapeDeviation = shapeDeviation;
	}	
//	public void setProfileLogKL(double[] logKL_plus, double[] logKL_minus){
//		this.logKL_plus = logKL_plus;
//		this.logKL_minus = logKL_minus;
//	}

	public void setIpCtrlLogKL(int cond, double logKL_plus, double logKL_minus){
		ipCtrl_logKL_plus[cond] = logKL_plus;
		ipCtrl_logKL_minus[cond] = logKL_minus;
	}
	/** Set control read count for the specified condition */
	public void setControlReadCounts(double controlCount, int cond) {
		if (unScaledControlCounts == null){
			unScaledControlCounts = new double[numConditions];
		}
		unScaledControlCounts[cond] = controlCount;
	}
	
	/** Set the local lambda (expected IP count) for Poisson test */
	public void setExpectedCounts(double expectedCount, int cond) {
		if (expectedCounts == null){
			expectedCounts = new double[numConditions];
		}
		expectedCounts[cond] = expectedCount;
	}
	
	public void setCondSignificance(int cond, boolean val) { condSignificance[cond] = val; }

	// overall average logKL
//	public double getAverageLogKL(){
//		double sum=0;
//		for(int c=0; c<numConditions; c++){
//			sum += logKL_plus[c] + logKL_minus[c];
//		}
//		return sum/(2*numConditions);
//	}
	// overall average logKL for control
	public double getAverageIpCtrlLogKL(){
		double sum=0;
		for(int c=0; c<numConditions; c++){
			sum += ipCtrl_logKL_plus[c] + ipCtrl_logKL_minus[c];
		}
		return sum/(2*numConditions);
	}	
	/**
	 * the average logKL and symmetric score (of +/- strands) for given condition
	 * @param cond
	 * @return average logKL and symmetric score
	 */
//	public Pair<Double, Double> getLogKL(int cond){
//		return new Pair<Double, Double>((logKL_plus[cond] + logKL_minus[cond])/2,
//			Math.abs(logKL_plus[cond] - logKL_minus[cond])/(logKL_plus[cond] + logKL_minus[cond]));
//	}
	
	public static void setSortingCondition(int cond){
		sortingCondition = cond;
	}

	public boolean onSameChrom(ComponentFeature f){
		return position.getChrom().equalsIgnoreCase(f.getPosition().getChrom());
	}
	
	/** Comparable default method, sort by location*/
	public int compareTo(ComponentFeature f) {
		return compareByLocation(f);
	}
	
	// The following methods can be called from creating a new Comparator.
	// See benjaminiHochbergCorrection() in BindingMixture for example use
	/** 
	 * sort componentFeatures in decreasing binding strength
	 */
	public int compareByTotalResponsibility(ComponentFeature f) {
		if(totalSumResponsibility>f.getTotalEventStrength()){return(-1);}
		else if(totalSumResponsibility<f.getTotalEventStrength()){return(1);}
		else return(0);
	}
	/** 
	 * sort componentFeatures first by having a k-mer group match, then by decreasing binding strength
	 */
	public int compareByTotalResponsibilityWithKmerMatch(ComponentFeature f) {
		if(kmerGroup!=null && f.getKmerGroup()==null){return(-1);}
		else if(kmerGroup==null && f.getKmerGroup()!=null){return(1);}
		else return compareByTotalResponsibility(f);
	}
	public int compareByCondResponsibility(ComponentFeature f) {
		if(condSumResponsibility[sortingCondition]>f.getCondSumResponsibility(sortingCondition)){return(-1);}
		else if(condSumResponsibility[sortingCondition]<f.getCondSumResponsibility(sortingCondition)){return(1);}
		else return(0);
	}
//	public int compareBylogKL(ComponentFeature f) {
//		double diff = getAverageLogKL() - f.getAverageLogKL();
//		return diff==0?0:(diff<0)?-1:1; // smaller logKL, more similar distribution
//	}	
	public int compareByLocation(ComponentFeature f) {
		return getPeak().compareTo(f.getPeak());
	}
	/** Compare p-value for expt with control <br>
	 *  Also calculate possion p-value based on IP only<br>
	 *  When comparing, use the larger p-value */
	public int compareByPValue(ComponentFeature f) {		
		if(getPValue(sortingCondition)<f.getPValue(sortingCondition)){return(-1);}
		else if(getPValue(sortingCondition)>f.getPValue(sortingCondition)){return(1);}
		else return(compareByCondResponsibility(f)); //P-values equal, compare by responsibility now
	}
	public int compareByPValue_wo_ctrl(ComponentFeature f) {
		if(p_values_wo_ctrl[sortingCondition]<f.getPValue_wo_ctrl(sortingCondition)){return(-1);}
		else if(p_values_wo_ctrl[sortingCondition]>f.getPValue_wo_ctrl(sortingCondition)){return(1);}
		else return(compareByCondResponsibility(f)); //P-values equal, compare by responsibility now
	}
	public int compareByQValue(ComponentFeature f) {
		if(q_value_log10[sortingCondition]<f.getQValueLog10(sortingCondition)){return(-1);}
		else if(q_value_log10[sortingCondition]>f.getQValueLog10(sortingCondition)){return(1);}
		else return(compareByCondResponsibility(f)); //Q-values equal, compare by responsibility now
	}
	public int compareByMfold(ComponentFeature f) {
		if(mfold < f.get_mfold())      { return -1; }
		else if(mfold > f.get_mfold()) { return  1; }
		else                           { return  0; }
	}
	public int compareByAvgQValue(ComponentFeature f) {
		double q_avg = 0;
		double q_avg2 = 0;
		for (int i=0;i<numConditions;i++){
			q_avg+=q_value_log10[i];
			q_avg2+=f.getQValueLog10(i);
		}
		double diff = q_avg-q_avg2;
		return diff==0?0:(diff<0)?-1:1;
	}	
	//Print the sequence
	public String toSequence(int win) {
		String seq;
		SequenceGenerator seqgen = new SequenceGenerator();
		Region peakWin=null;
	
		int start = position.getLocation()-((int)(win/2));
		if(start<1){start=1;}
		int end = position.getLocation()+((int)win/2)-1;
		if(end>coords.getGenome().getChromLength(coords.getChrom())){end =coords.getGenome().getChromLength(coords.getChrom())-1;} 
		peakWin = new Region(coords.getGenome(), coords.getChrom(), start, end);
		
		String seqName = new String(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth());
		String sq = seqgen.execute(peakWin);
		seq = seqName +"\n"+ sq+"\n";
		return seq;
	}
	public double getExptCtrlRatio(int condition){
		double ratio = 200;
		if (unScaledControlCounts[condition]!=0)
			ratio = getTotalEventStrength()*conditionBeta[condition]/unScaledControlCounts[condition];
		return ratio;
	}
	
	public double getAvgExptCtrlRatio(){
		int ratios=0;
		for(int c=0; c<numConditions; c++){
			ratios+=getExptCtrlRatio(c);
		}
		return ratios/numConditions;
	}
	/** Return p-value for expt with control for the specified condition<br>
	 *  the larger value of binomial p-value and possion p-value based on IP only */
	public double getPValue(int cond){
		return Math.max(p_values[cond], p_values_wo_ctrl[cond]);
	}
	public double getPValue_w_ctrl(int cond){
		return p_values[cond];
	}
	public double getPValue_wo_ctrl(int cond){
		return p_values_wo_ctrl[cond];
	}	
	public double getQValueLog10(int cond){
		return q_value_log10[cond];
	}
	public void setPValue_w_ctrl(double pValue, int cond){
		p_values[cond] = pValue;
	}
	public void setPValue_wo_ctrl(double pValue, int cond){
		p_values_wo_ctrl[cond] = pValue;
	}
	//Corrected_Q_value
	public void setQValueLog10(double qValue, int cond){
		q_value_log10[cond] = qValue;
	}

	//Print the feature
	//each field should match header String
	public String toString() {
		StringBuilder result = new StringBuilder();
		String strandStr = "";
		if (strand=='-')
			strandStr = ":-";
		else if (strand=='+')
			strandStr = ":+";
		result.append(position.getLocationString()).append(strandStr).append("\t");
		result.append(String.format("%7.1f\t", totalSumResponsibility));
			
        for(int c=0; c<numConditions; c++){
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		result.append(String.format("%d\t", condSignificance[c] ? 1 : 0 ));
        		result.append(String.format("%7.1f\t", getEventReadCounts(c) ));
        	}
        	if(unScaledControlCounts!=null){
        		double fold = 0;
        		if (getScaledControlCounts(c)==0)
        			fold = 9999.9;
        		else
        			fold = getEventReadCounts(c)/getScaledControlCounts(c);
        		result.append(String.format("%7.1f\t", getScaledControlCounts(c)))
        			  .append(String.format("%7.1f\t", fold));
        	}
        	else
        		result.append("NaN\t").append("NaN\t");
        
        	double q_lg = getQValueLog10(c);
    		if (q_lg==Double.POSITIVE_INFINITY)
    			q_lg= 999;
        	result.append(String.format("%7.2f\t", q_lg));    
        	double p_lg = -Math.log10(getPValue(c));
        	if(unScaledControlCounts==null)
        		p_lg = -Math.log10(getPValue(c));
        	if (p_lg==Double.POSITIVE_INFINITY)
    			p_lg= 999;
        	result.append(String.format("%7.2f\t", p_lg));
        	
    		result.append(String.format("%7.2f\t", getShapeDeviation(c)));
    		
        	if(unScaledControlCounts!=null)
        		result.append(String.format("%7.2f\t", getAverageIpCtrlLogKL()));
        	else
        		result.append("NaN\t");
        }

        result.append(isJointEvent?1:0).append("\t");
         
		String gene = nearestGene == null ? "NONE" : nearestGene.getName();
		result.append(gene).append("\t");
		result.append(distToGene).append("\t");
		
        if (annotations != null) {
            for (Region r : annotations) {
                result.append(r.toString() + ",");
            }
            result.append("\t");
        }
        result.append(String.format("%7.1f\t", alpha));
        result.append(EM_position.getLocationString());
        result.append("\n");

		return result.toString();
	}

	//GFF3 (fill in later)
	public String toGFF(){
		return("");
	}
	
	//generate Header String, each field should match toString() output
	public String headString(){
		StringBuilder header = new StringBuilder();
		
		header.append("Position\t")
			  .append("   IP\t");
        
        for(int c=0; c<numConditions; c++){
        	String name = numConditions==1?"":conditionNames.get(c)+"_";
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		header.append(name+"Present\t");
        		header.append(name+"IP\t");
        	}
        	header.append(name+"Control\t")
        		  .append(name+"   Fold\t")
        	      .append(name+"Q_-lg10\t")
  	      		  .append(name+"P_-lg10\t")
  	      		  .append("IPvsEMP\t")
  	      		  .append(name+"IPvsCTR\t");
        }
        
        header.append("Joint\t");
		header.append("NearestGene\t").append("Distance\t");
		
        if (annotations != null) {
            header.append("Annotations").append("\t");
        }
        header.append("Alpha\t").append("EM_Position");
        header.append("\n");
        return header.toString();
	}
	
	//generate Header String, each field should match toString() output
	// for GPS release v1
	public String headString_v1(){
		StringBuilder header = new StringBuilder();
		
		header.append("Position\t")
			  .append("     IP\t");        
        
        for(int c=0; c<numConditions; c++){
        	String name = numConditions==1?"":conditionNames.get(c)+"_";
        	if (numConditions!=1) {		// if single condition, IP is same as total
        		header.append(name+"Present\t");
        		header.append(name+"IP\t");
        	}
        	header.append(name+"Control\t")
  		  		  .append(name+(numConditions==1?"   ":"")+"Fold\t")
        	      .append(name+"Expectd\t")
        	      .append(name+"Q_-lg10\t")
  	      		  .append(name+"P_-lg10\t")
  	      		  .append(name+"P_poiss\t")
  	      		  .append(name+"IPvsEMP\t")
  	      		  .append(name+"  Noise\t");
        	if (c<numConditions-1)
        		header.append("\t");
        }
        if (boundSequence!=null){
	        header.append("KmerGroup  \tMotifId\tKG_hgp\tStrand"); 
        }
        header.append("\n");
        return header.toString();
	}
	//Print the feature
	//each field should match header String
	public String toString_v1() {
		StringBuilder result = new StringBuilder();
		String strandStr = "";
		if (strand=='-')
			strandStr = ":-";
		else if (strand=='+')
			strandStr = ":+";
		result.append(position.getLocationString()).append(strandStr).append("\t");

		result.append(String.format("%7.1f\t", totalSumResponsibility));
      
        for(int c=0; c<numConditions; c++){
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		result.append(String.format("%d\t", condSignificance[c] ? 1 : 0 ));
        		result.append(String.format("%7.1f\t", getEventReadCounts(c) ));
        	}
        	if(unScaledControlCounts!=null){
        		double fold = getEventReadCounts(c)/getScaledControlCounts(c);
        		if (fold>9999.9)
        			fold = 9999.9;
        		result.append(String.format("%7.1f\t", getScaledControlCounts(c)))
        			  .append(String.format("%7.1f\t", fold));
        	} else {
        		result.append("NaN\t").append("NaN\t");                
            }  
        	if(expectedCounts!=null){
        		result.append(String.format("%7.1f\t", getExpectedCount(c)));
        	}else {
        		result.append("NaN\t");                
            }
        	double q_lg = getQValueLog10(c);
    		if (q_lg==Double.POSITIVE_INFINITY)
    			q_lg= 999;
        	result.append(String.format("%7.2f\t", q_lg));   
        	if(unScaledControlCounts!=null){
	        	double p_lg1 = -Math.log10(getPValue_w_ctrl(c));
	        	if (p_lg1==Double.POSITIVE_INFINITY)
	    			p_lg1= 999;
	        	result.append(String.format("%7.2f\t", p_lg1)); 
        	}
        	else
        		result.append("NaN\t");
        	
        	double p_lg2 = -Math.log10(getPValue_wo_ctrl(c));
        	if (p_lg2==Double.POSITIVE_INFINITY)
        		p_lg2= 999;
        	result.append(String.format("%7.2f\t", p_lg2)); 

    		result.append(String.format("%7.2f\t", getShapeDeviation(c))); 
    		result.append(String.format("%7.2f\t", getNoiseFraction()));
        }
        if (boundSequence!=null){
	        if (kmerGroup!=null)
	        	result.append(String.format("%s_%d/%d\t%d\t%.2f\t%s\t%s", kmerGroup.getBestKmer().getKmerString(),kmerGroup.getGroupHitCount(), kmerGroup.getGroupNegHitCount(), kmerGroup.getClusterId(), kmerGroup.getHgp(), kmerStrand, boundSequence));
	        else
	        	result.append(CommonUtils.padding(8, '-')).append("\t-1\t0.00\t*\t").append(boundSequence);
	        result.append("\n");
        }
        else
        	CommonUtils.replaceEnd(result, '\n');
		return result.toString();
	}
	
	//Print the feature in BED format
	public String toBED() {
		StringBuilder bed = new StringBuilder();
		bed.append("chr"+position.getChrom()).append("\t");
		Region r = position.expand(100);
		bed.append(r.getStart()).append("\t");
		bed.append(r.getEnd()).append("\t");
		bed.append(position.getLocationString()).append("\t");
		bed.append(String.format("%7.1f\t", totalSumResponsibility));
		if (strand=='-')
			bed.append("\t-");
		else if (strand=='+')
			bed.append("\t+");
        bed.append("\n");
		return bed.toString();
	}
	
	/** 
	 * Output ENCODE narrowPeak format<br>
	 * starting from GEM release v2.0<br>
	 * The score is based on q-value (corrected from the larger p-value from binomial p-value (IP vs Control) and Poisson p-value (IP only)
	 * @return
	 */
	public String toNarrowPeak(int score){
		double q_lg = getQValueLog10(0);
		if (q_lg==Double.POSITIVE_INFINITY)
			q_lg= 999;   
		double p_lg1=0;
		if(unScaledControlCounts!=null){
        	p_lg1 = -Math.log10(getPValue_w_ctrl(0));
        	if (p_lg1==Double.POSITIVE_INFINITY)
    			p_lg1= 999;
    	}
    	double p_lg2 = -Math.log10(getPValue_wo_ctrl(0));
    	if (p_lg2==Double.POSITIVE_INFINITY)
    		p_lg2= 999;
		String strandStr = ".";
		if (strand=='-')
			strandStr = "-";
		else if (strand=='+')
			strandStr = "+";

		return "chr"+position.getChrom()+"\t"+(position.getLocation()-100)+"\t"+(position.getLocation()+101)+"\t"+
    				position.toString()+"\t"+score+"\t"+strandStr+"\t"+String.format("%.1f",totalSumResponsibility)+"\t"+
    				String.format("%.2f", p_lg1>p_lg2?p_lg2:p_lg1)+"\t"+String.format("%.2f", q_lg)+"\t"+100+"\n";
	}
	
	//generate Header String for BED file, each field should match BED() output
	// for GPS release v1
	public String headString_BED(){
		StringBuilder header = new StringBuilder();
		
		header.append("Chrom\t")
			  .append("Coord-5\t")
			  .append("Coord+5\t")
			  .append("Coordinate\t")
			  .append("     IP\t")
			  .append("Strand\t");        
        
        for(int c=0; c<numConditions; c++){
        	String name = numConditions==1?"":conditionNames.get(c)+"_";
        	if (numConditions!=1) {		// if single condition, IP is same as total
        		header.append(name+"Present\t");
        		header.append(name+"IP\t");
        	}
        	header.append(name+"Control\t")
  		  		  .append(name+"IP/Ctrl\t")
        	      .append(name+"Q_-lg10\t")
  	      		  .append(name+"P_-lg10\t")
  	      		  .append("  Shape");
        	if (c<numConditions-1)
        		header.append("\t");
        }
        header.append("\n");
        return header.toString();
	}
	public void releaseMemory(){
		boundSequence = null;
		conditionBeta = null;
//		logKL_plus = null;
//		logKL_minus = null;
		ipCtrl_logKL_plus = null;
		ipCtrl_logKL_minus = null;
		shapeDeviation = null;
		kmerGroup = null;		
	}
	
	public static void setNon_specific_ratio(double[] non_specific_ratio) {
		ComponentFeature.non_specific_ratio = non_specific_ratio;
	}

	public static void setConditionNames(ArrayList<String> conditionNames) {
		ComponentFeature.conditionNames = conditionNames;
		numConditions = conditionNames.size();
		non_specific_ratio = new double[numConditions];
		for (int i=0;i<numConditions;i++)
			non_specific_ratio[i] = 1.0;
	}
	

}
