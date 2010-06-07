package edu.mit.csail.cgs.deepseq;

import edu.mit.csail.cgs.datasets.general.Region;

public abstract class BackgroundModel {
	protected int modelType; //-1 for genome-wide, 0 for current region, positive integer for local window size  
	protected double logConfidence;
	protected double confThreshold;
	protected double totalReads, readLength, regionLength, mappableRegion, binWidth, binStep;
	protected int countThreshold=0;
	protected char strand='.';
	protected boolean useThisExpt=true;
	protected double scaling=1;
	
	public BackgroundModel(int mtype, double lc, double r, double l, double rl, double mr, double bw, double bo){this(mtype,lc,r,l,rl,mr,bw,bo,'.',1, true);}
	public BackgroundModel(int mtype, double lc, double r, double l, double rl, double mr, double bw, double bo, char str, double sc, boolean ute){
		modelType = mtype;
		logConfidence=lc;
		confThreshold = Math.pow(10,logConfidence);
		totalReads=r;
		readLength=l;
		regionLength=rl;
		mappableRegion=mr;
		binWidth=bw;
		binStep=bo;	
		strand=str;
		scaling=sc;	
		useThisExpt=ute;
	}
	
	//Accessors
	public char getStrand(){return strand;}
	public boolean isGenomeWide(){return modelType==-1 ? true : false;}
	public int getThreshold(){return countThreshold;}
	
	//Required
	public abstract boolean passesThreshold(int count);
	public abstract boolean underThreshold(int count);
	protected abstract int calcCountThreshold();
	
	//Update the threshold (depends on the type of model... local, etc)
	public void updateModel(Region currReg, int currOffset, double [] thisExptHitCounts, double [] otherExptHitCounts){
		double [] stackedHitCounts = useThisExpt ? thisExptHitCounts : otherExptHitCounts; 
		if(modelType==-1){//Genome-wide
			if(countThreshold==0){
				countThreshold = calcCountThreshold();
			}
		}else if(modelType==0){//Current region
			double sum=0;
			for(int i=0; i<stackedHitCounts.length; i++)
				sum+=stackedHitCounts[i];
			totalReads = scaling*(sum/(((readLength)/binStep)+1));
			regionLength = currReg.getWidth();
			mappableRegion=1.0;//Big assumption
			countThreshold = calcCountThreshold();
		}else{//Window around current position
			int win = modelType;
			int istart = currOffset-(win/2)<0 ? 0 :currOffset-(win/2);
			int istop = currOffset+(win/2)>=currReg.getWidth() ? currReg.getWidth()-1 :currOffset+(win/2);
			int winTrueSize=0;
			double sum=0;
			for(int i=istart; i<=istop; i++){
				int binid = (int)Math.max(0, ((double)(i)/binStep));
				sum+=stackedHitCounts[binid]; 
				winTrueSize+=binStep;
			}
			totalReads = scaling*(sum/(((readLength)/binStep)+1));
			regionLength=winTrueSize;
			//mappableRegion=1.0; //any need for this assumption?
			countThreshold = calcCountThreshold();
		}
	}	
}
