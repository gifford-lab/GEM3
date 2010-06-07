package edu.mit.csail.cgs.deepseq;

import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * <tt>ReadHit</tt> is a type of <tt>StrandedRegion</tt> with a <tt>weight</tt> <br>
 * It represents the stranded region of the referenced genome where a read maps to <br>
 * <u>Note</u>: <br>
 *              1) For a single read, we can have multiple (read) hits <br>
 *              2) The <tt>start < end</tt> of the <tt>ReadHit</tt> ALWAYS   <br>
 *                 strand <tt>str</tt> is used to indicate whether the read is
 *                 forward or reverse.
 * @author shaun
 *
 */
public class ReadHit extends StrandedRegion{

	protected static int hitcount=1;
	protected int id=1;
	protected int misMatch; //For now, misMatch is only used to filter sub-optimal results out of certain file formats
	protected double weight;
		
	public ReadHit(Genome g, int id,String c, int s, int e, char str){this(g,id,c,s,e,str,1.0, 0);}
	public ReadHit(Genome g, int id,String c, int s, int e, char str, double w){this(g,id,c,s,e,str,w, 0);}
	public ReadHit(Genome g, int id,String c, int s, int e, char str, double w, int m){
		super(g,c,s,e>s ? e:s+1,str);
		this.id= id==-1 ? ++this.hitcount : id;
		weight = w;
		misMatch=m;
	}
	
	//Accessors
	public int getID(){return(id);}
	public double getWeight(){return(weight);}
	public int getMisMatch(){return(misMatch);}
	public void setWeight(double w){weight=w;}
	public double getReadLength(){return(this.getEnd()-this.getStart()+1);}
}
