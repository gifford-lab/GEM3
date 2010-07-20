package edu.mit.csail.cgs.deepseq.features;

import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class ComponentFeature extends Feature  implements Comparable<ComponentFeature>{
	private static double non_specific_ratio[];
	private static int sortingCondition;	// the condition to compare p-values
	private static ArrayList<String> conditionNames;
	protected Point position;
	protected double mixingProb;
	protected double mfold;
	protected double[] conditionBeta;
	protected boolean[] condSignificance;
	protected int numConditions=0;
	protected double[] condSumResponsibility;
	protected double totalSumResponsibility=0;
//	protected double[][] readProfile_plus;
//	protected double[][] readProfile_minus;
	protected double logKL_plus[];
	protected double logKL_minus[];
	protected double logKL_ctrl_plus[];
	protected double logKL_ctrl_minus[];
	protected double unScaledControlCounts[];
	protected double p_values[];
	protected double p_values_wo_ctrl[];
	protected double q_value_log10[];
	protected double shape_z_weighted[];
	protected double shape_z_scores[];
	// the sum of ranking by control based p-value and peak shape parameter
	// each condition will have its contribution to the ranking sum
	protected int rank_sum=0; 			
	protected Point EM_position;		//  EM result
	protected double alpha;
	protected boolean use_internal_em_train;
	
	
	public ComponentFeature(BindingComponent b, boolean use_internal_em){
		this(b);
		this.use_internal_em_train = use_internal_em;
	}

	public ComponentFeature(BindingComponent b){
		super(null);
		position = b.getLocation();
		coords = position.expand(1);
		mixingProb = b.getMixProb();
		numConditions = b.getNumConditions();
		
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
		shape_z_weighted = new double[numConditions];
		shape_z_scores = new double[numConditions];
		EM_position = b.getEMPosition();
		alpha = b.getAlpha();
		logKL_ctrl_plus = new double[numConditions];
		logKL_ctrl_minus = new double[numConditions];
	}
	
	//Accessors 
	public double getMixProb(){return(mixingProb);}
	public double getAlpha(){return(alpha);}
	public Point getPosition() { return position;}
	public Point getEMPosition() { return EM_position;}

	public double[] getCondBetas() { return conditionBeta; }
	public boolean[] getCondSignificance() { return condSignificance; }
	public double getTotalSumResponsibility(){return(totalSumResponsibility);}
	public double getCondSumResponsibility(int cond){return condSumResponsibility[cond];}
	public double[] getUnscaledControlCounts(){return(unScaledControlCounts);}
	public double getScaledControlCounts(int cond){
		return unScaledControlCounts[cond]*non_specific_ratio[cond];
	}
	public double getEventReadCounts(int cond){
		if(use_internal_em_train) { return totalSumResponsibility*conditionBeta[cond]; }
		else                      { return condSumResponsibility[cond];                }
	}
	public double get_mfold() { return mfold; }
	public void set_mfold(double mf) { mfold = mf; }
	
//	public double getReadProfile(int cond, int index, char strand){
//		double result=0;
//		if (strand=='+')
//			result=readProfile_plus[cond][index];
//		else if (strand=='-')
//			result=readProfile_minus[cond][index];
//		return result;
//	}
//	public void calc_logKL_Divergence(double[] model, double[]smootherKernel){
//		logKL_plus = new double[numConditions];
//		logKL_minus = new double[numConditions];
//		for(int c=0; c<numConditions; c++){
//			logKL_plus[c] = StatUtil.log_KL_Divergence(model, StatUtil.symmetricKernelSmoother(readProfile_plus[c],smootherKernel));
//			logKL_minus[c] = StatUtil.log_KL_Divergence(model, StatUtil.symmetricKernelSmoother(readProfile_minus[c],smootherKernel));
//		}
//	}
	public void setProfileLogKL(double[] logKL_plus, double[] logKL_minus){
		this.logKL_plus = logKL_plus;
		this.logKL_minus = logKL_minus;
	}
//	public void clearReadProfiles(){
//		readProfile_plus = null;
//		readProfile_minus = null;
//	}
	public void setCtrlProfileLogKL(int cond, double logKL_plus, double logKL_minus){
		logKL_ctrl_plus[cond] = logKL_plus;
		logKL_ctrl_minus[cond] = logKL_minus;
	}
	public void setControlReadCounts(double controlCount, int cond) {
		if (unScaledControlCounts == null){
			unScaledControlCounts = new double[numConditions];
		}
		unScaledControlCounts[cond] = controlCount;
	}
	
	public void setCondSignificance(int cond, boolean val) { condSignificance[cond] = val; }

	// overall average logKL
	public double getAverageLogKL(){
		double sum=0;
		for(int c=0; c<numConditions; c++){
			sum += logKL_plus[c] + logKL_minus[c];
		}
		return sum/(2*numConditions);
	}
	// overall average logKL for control
	public double getAverageCtrlLogKL(){
		double sum=0;
		for(int c=0; c<numConditions; c++){
			sum += logKL_ctrl_plus[c] + logKL_ctrl_minus[c];
		}
		return sum/(2*numConditions);
	}	
	/**
	 * the average logKL and symmetric score (of +/- strands) for given condition
	 * @param cond
	 * @return average logKL and symmetric score
	 */
	public Pair<Double, Double> getLogKL(int cond){
		return new Pair<Double, Double>((logKL_plus[cond] + logKL_minus[cond])/2,
			Math.abs(logKL_plus[cond] - logKL_minus[cond])/(logKL_plus[cond] + logKL_minus[cond]));
	}
	
	public static void setSortingCondition(int cond){
		sortingCondition = cond;
	}
	
	//Comparable default method
	public int compareTo(ComponentFeature f) {
//		return compareByAvgQValue(f);
		return compareByLocation(f);
	}
	
	// The following method can be called from creating a new Comparator.
	// See benjaminiHochbergCorrection() in BindingMixture for example use
	public int compareByRankSum(ComponentFeature f) {
		double diff = getRank_sum()-f.getRank_sum();
		return diff==0?0:(diff<0)?-1:1;
	}
	public int compareByTotalResponsibility(ComponentFeature f) {
		if(totalSumResponsibility>f.getTotalSumResponsibility()){return(-1);}
		else if(totalSumResponsibility<f.getTotalSumResponsibility()){return(1);}
		else return(0);
	}
	public int compareBylogKL(ComponentFeature f) {
		double diff = getAverageLogKL() - f.getAverageLogKL();
		return diff==0?0:(diff<0)?-1:1; // smaller logKL, more similar distribution
	}	
	public int compareByPeakShape(ComponentFeature f) {
		if(shape_z_weighted[sortingCondition]<f.getShapeParameter(sortingCondition)){return(-1);}
		else if(shape_z_weighted[sortingCondition]>f.getShapeParameter(sortingCondition)){return(1);}
		else return(0);
	}
	public int compareByLocation(ComponentFeature f) {
		return getPeak().compareTo(f.getPeak());
	}
	public int compareByPValue(ComponentFeature f) {
		if(p_values[sortingCondition]<f.getPValue(sortingCondition)){return(-1);}
		else if(p_values[sortingCondition]>f.getPValue(sortingCondition)){return(1);}
		else return(0);
	}
	public int compareByPValue_wo_ctrl(ComponentFeature f) {
		if(p_values_wo_ctrl[sortingCondition]<f.getPValue_wo_ctrl(sortingCondition)){return(-1);}
		else if(p_values_wo_ctrl[sortingCondition]>f.getPValue_wo_ctrl(sortingCondition)){return(1);}
		else return(0);
	}
	public int compareByQValue(ComponentFeature f) {
		if(q_value_log10[sortingCondition]<f.getQValueLog10(sortingCondition)){return(-1);}
		else if(q_value_log10[sortingCondition]>f.getQValueLog10(sortingCondition)){return(1);}
		else return(0);
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
			ratio = getTotalSumResponsibility()*conditionBeta[condition]/unScaledControlCounts[condition];
		return ratio;
	}
	
	public double getAvgExptCtrlRatio(){
		int ratios=0;
		for(int c=0; c<numConditions; c++){
			ratios+=getExptCtrlRatio(c);
		}
		return ratios/numConditions;
	}
	public double getPValue(int cond){
		return p_values[cond];
	}
	public double getPValue_wo_ctrl(int cond){
		return p_values_wo_ctrl[cond];
	}	
	public double getQValueLog10(int cond){
		return q_value_log10[cond];
	}
	public void setPValue(double pValue, int cond){
		p_values[cond] = pValue;
	}
	public void setPValue_wo_ctrl(double pValue, int cond){
		p_values_wo_ctrl[cond] = pValue;
	}
	//Corrected_Q_value
	public void setQValueLog10(double qValue, int cond){
		q_value_log10[cond] = qValue;
	}
	public double getShapeParameter(int cond){
		return shape_z_weighted[cond];
	}
	public void setShapeParameter(double pValue, int cond){
		shape_z_weighted[cond] = pValue;
	}
	public double getShapeZScore(int cond){
		return shape_z_scores[cond];
	}
	public void setShapeZScore(double zscore, int cond){
		shape_z_scores[cond] = zscore;
	}
	//Print the feature
	//each field should match header String
	public String toString() {
		StringBuilder result = new StringBuilder();
		
		result.append(position.getLocationString()).append("\t");
//		result.append(rank_sum).append("\t");
		result.append(String.format("%.1f\t", totalSumResponsibility));
//		result.append(EM_position.getLocationString()).append("\t");
		result.append(String.format("%.2f\t", getAverageLogKL()));
				
//        for(int c=0; c<numConditions; c++){
//        	Pair<Double, Double> logKL = getLogKL(c);
//        	if (numConditions>1){
//	        	result.append(String.format("%.4f\t", conditionBeta[c]));
//	        	result.append(String.format("%.2f\t", logKL.car()));
//        	}
//        	result.append(String.format("%.4f\t", getShapeZScore(c)));
//        	result.append(String.format("%.6f\t", getShapeParameter(c)));
//        	result.append(String.format("%.2f\t", logKL.cdr()));
//        }

//        if (!(unScaledControlCounts==null)){
//	        for(int c=0; c<numConditions; c++){
//	        	if (numConditions!=1)	// if single condition, IP is same as total
//	        		result.append(String.format("%.1f\t", getEventReadCounts(c) ));
//	        	result.append(String.format("%.1f\t", getScaledControlCounts(c)));
//	        	result.append(String.format("%.2f\t", getQValueLog10(c)));
//	        	result.append(String.format("%.2f\t", -Math.log10(getPValue(c))));
//	        }
//        }
//        else{
//	        for(int c=0; c<numConditions; c++){
//	        	result.append(String.format("%.1f\t", getEventReadCounts(c) ));
//	        }
//	        result.append("NA\t");
//        }
       
        for(int c=0; c<numConditions; c++){
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		result.append(String.format("%c\t", condSignificance[c] ? 'T' : 'F' ));
        		result.append(String.format("%.1f\t", getEventReadCounts(c) ));
        	}
        	if(unScaledControlCounts!=null)
        		result.append(String.format("%.1f\t", getScaledControlCounts(c)));
        	else
        		result.append("NA\t");
        
        	result.append(String.format("%.2f\t", getQValueLog10(c)));
        	
        	if(unScaledControlCounts!=null)
        		result.append(String.format("%.2f\t", -Math.log10(getPValue(c))));
        	else
        		result.append(String.format("%.2f\t", -Math.log10(getPValue_wo_ctrl(c))));
        }

        result.append(mixingProb==1?1:0).append("\t");
//        result.append(String.format("%.4f\t", mixingProb));
         
		String gene = nearestGene == null ? "NONE" : nearestGene.getName();
		result.append(gene).append("\t");
		result.append(distToGene).append("\t");
		
        if (annotations != null) {
            for (Region r : annotations) {
                result.append(r.toString() + ",");
            }
            result.append("\t");
        }
        result.append(String.format("%.1f\t", alpha));
        result.append(String.format("%.2f\t", getAverageCtrlLogKL()));
        result.append(EM_position.getLocationString() + "\t");
        result.append("\n");

		return result.toString();
	}
	//GFF3 (fill in later)
	public String toGFF(){
		return("");
	}
	
	//generate Header String, each field should match toString() output
	public String headString(){
		StringBuilder header = new StringBuilder("%");
		
		header.append("Position\t")
//			  .append("Rank_Sum\t")
			  .append("IpStrength\t")
//			  .append("EM_Posi\t")
			  .append("Shape\t");
				
//        for(int c=0; c<numConditions; c++){
//        	String name = numConditions==1?"":conditionNames.get(c)+"_";
//        	if (numConditions>1)
//        		header.append(name+"Proportion\t")
//        	      	  .append(name+"Shape\t");
//    	    header.append("Shape_Z\t");
//    	    header.append("Shape_Param\t");
//    	    header.append(name+"ShapeAsymmetry\t");
//        }

//        if (!(unScaledControlCounts==null)){
//	        for(int c=0; c<numConditions; c++){
//	        	String name = numConditions==1?"":conditionNames.get(c)+"_";
//	        	if (numConditions!=1)		// if single condition, IP is same as total
//	        		header.append(name+"IpStrength\t");
//	        	header.append(name+"CtrlStrength\t")
//      	      		  .append(name+"Q_value_log10\t")
//      	      		  .append(name+"P_value_log10\t");
//	        }
//        }
//        else{
//	        for(int c=0; c<numConditions; c++){
//	        	String name = numConditions==1?"":conditionNames.get(c)+"_";
//	        	header.append(name+"IpStrength\t");
//	        }
//	        header.append("Control\t");
//        }
//        
        
        for(int c=0; c<numConditions; c++){
        	String name = numConditions==1?"":conditionNames.get(c)+"_";
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		header.append(name+"Significant\t");
        		header.append(name+"IpStrength\t");
        	}
        	header.append(name+"CtrlStrength\t")
        	      .append(name+"Q_value_log10\t")
  	      		  .append(name+"P_value_log10\t");
        }
        
        header.append("UnaryEvent\t");
		header.append("NearestGene\t").append("Distance\t");
		
        if (annotations != null) {
            header.append("Annotations").append("\t");
        }
        header.append("Alpha");
        header.append("\t");
        header.append("ControlLogKL");
        header.append("\t");
        header.append("EM_Position");
        header.append("\t");
        header.append("\n");
        return header.toString();
	}
	
	//generate Header String, each field should match toString() output
	// for GPS release v1
	public String headString_v1(){
		StringBuilder header = new StringBuilder("");
		
		header.append("Position\t")
			  .append("IpStrength\t");       
        
        for(int c=0; c<numConditions; c++){
        	String name = numConditions==1?"":conditionNames.get(c)+"_";
        	if (numConditions!=1) {		// if single condition, IP is same as total
        		header.append(name+"Significant\t");
        		header.append(name+"IpStrength\t");
        	}
        	header.append(name+"CtrlStrength\t")
        	      .append(name+"Q_value_log10\t")
  	      		  .append(name+"P_value_log10\t");
        }
        header.append("\n");
        return header.toString();
	}
	//Print the feature
	// for GPS release v1
	//each field should match header String
	public String toString_v1() {
		StringBuilder result = new StringBuilder();
		
		result.append(position.getLocationString()).append("\t");
		result.append(String.format("%.1f\t", totalSumResponsibility));
       
        for(int c=0; c<numConditions; c++){
        	if (numConditions!=1) {	// if single condition, IP is same as total
        		result.append(String.format("%c\t", condSignificance[c] ? 'T' : 'F' ));
        		result.append(String.format("%.1f\t", getEventReadCounts(c) ));
        	}
        	if(unScaledControlCounts!=null)
        		result.append(String.format("%.1f\t", getScaledControlCounts(c)));
        	else
        		result.append("NA\t");
        
        	result.append(String.format("%.2f\t", getQValueLog10(c)));
        	
        	if(unScaledControlCounts!=null)
        		result.append(String.format("%.2f\t", -Math.log10(getPValue(c))));
        	else
        		result.append(String.format("%.2f\t", -Math.log10(getPValue_wo_ctrl(c))));
        }

        result.append("\n");
		return result.toString();
	}
	public static void setNon_specific_ratio(double[] non_specific_ratio) {
		ComponentFeature.non_specific_ratio = non_specific_ratio;
	}

	public static void setConditionNames(ArrayList<String> conditionNames) {
		ComponentFeature.conditionNames = conditionNames;
	}

	public int getRank_sum() {
		return rank_sum;
	}

	public void addRank_sum(int rank_sum) {
		this.rank_sum += rank_sum;
	}
}
