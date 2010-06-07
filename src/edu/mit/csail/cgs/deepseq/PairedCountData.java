package edu.mit.csail.cgs.deepseq;

import edu.mit.csail.cgs.utils.models.Model;

public class PairedCountData extends Model{
	public Double x,y;
	public PairedCountData(double a, double b){
		x=a;
		y=b;
	}
}
