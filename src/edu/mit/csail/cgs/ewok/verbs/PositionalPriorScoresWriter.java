package edu.mit.csail.cgs.ewok.verbs;

import java.util.Formatter;
import java.util.Iterator;
import java.util.regex.*;
import java.io.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.iterators.SingleIterator;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.BayesLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;


public class PositionalPriorScoresWriter<X extends Region> implements Sink<X>, Closeable {

    private PrintStream ps;
    private FileOutputStream os;
    private int lineLength;
    private ChipChipBayes ccb;
    private double minProbeThreshold;
	private double minPrior = Double.MAX_VALUE;
	private double minPosterior = Double.MAX_VALUE;

    
    public PositionalPriorScoresWriter (String fname, Genome genome, String expt, String version, double minProbeThreshold) throws IOException, FileNotFoundException {
        os = new FileOutputStream(fname);
        ps = new PrintStream(os);
        lineLength = 100;

        //get the jbd output
        this.minProbeThreshold = minProbeThreshold;
		BayesLocator bloc = new BayesLocator(genome, expt, version);
		ccb = bloc.createObject();

    }

    public PositionalPriorScoresWriter(PrintStream ps, Genome genome, String expt, String version, double minProbeThreshold) {
        this.ps = ps;
        os = null;
        lineLength = 100;
        
		//get the jbd output
        this.minProbeThreshold = minProbeThreshold;
		BayesLocator bloc = new BayesLocator(genome, expt, version);
		ccb = bloc.createObject();
    }
    
    public void setLineLength(int ll) { lineLength = ll; }

	public void init() {}

    public void consume(Iterator<X> iter) {
        while (iter.hasNext()) {
            consume(iter.next());
        }
    	System.out.println("min posterior: " + minPosterior);
        System.out.println("min prior: " + minPrior);
    }

	public void finish() { close(); }

    public void consume(X r) {
        ps.println(">" + r.toString());
        try {
        	ccb.window(r.getChrom(), r.getStart(), r.getEnd());
        	if (ccb.getCount() < 2) {
        		throw new IllegalArgumentException();
        	}
        	
        	double[] posteriors = new double[r.getWidth()];
        	
        	int prevProbeIndex = 0;
        	int prevProbePos = ccb.getPos(prevProbeIndex);
        	double prevProbePosterior = ccb.getPosterior(prevProbeIndex);
        	if (prevProbePosterior < minProbeThreshold) {
        		prevProbePosterior = minProbeThreshold;
        	}
        	
        	int nextProbeIndex = 0;
        	int nextProbePos = prevProbePos;
        	double nextProbePosterior = prevProbePosterior;

        	if (prevProbePosterior < minPosterior) {
        		minPosterior = prevProbePosterior;
        	}
        	double maximumPosterior = prevProbePosterior;
        	double posteriorSum = 0;
        	
        	for(int i = 0; i < r.getWidth(); i++) {
        		int curpos = i + r.getStart();
        		if (curpos < ccb.getPos(0)) {
        			posteriors[i] = ccb.getPosterior(0);
        		}
        		else if (curpos > ccb.getPos(ccb.getCount() - 1)) {
        			posteriors[i] = ccb.getPosterior(ccb.getCount() - 1);
        		}
        		else if (curpos == nextProbePos) {
        			posteriors[i] = nextProbePosterior;
        			
        			prevProbeIndex = nextProbeIndex;
        			prevProbePos = nextProbePos;
        			prevProbePosterior = nextProbePosterior;
        			

        			if ((nextProbeIndex + 1) < ccb.getCount()) {
        				nextProbeIndex = nextProbeIndex + 1;
        				nextProbePos = ccb.getPos(nextProbeIndex);
        				nextProbePosterior = ccb.getPosterior(nextProbeIndex);
        	        	if (nextProbePosterior < minProbeThreshold) {
        	        		nextProbePosterior = minProbeThreshold;
        	        	}

        				if (nextProbePosterior > maximumPosterior) {
        					maximumPosterior = nextProbePosterior;
        				}
        				if (nextProbePosterior < minPosterior) {
        					minPosterior = nextProbePosterior;
        				}
        			}
        		}
        		else {
        			//curpos is in between two probes so interpolate...
        			double mix = ((double)(curpos - prevProbePos)) / ((double)(nextProbePos - prevProbePos));
        			posteriors[i] = prevProbePosterior + (mix * (nextProbePosterior - prevProbePosterior));
        		}
    			posteriorSum = posteriorSum + posteriors[i];
        	}
        	
        	double scaleFactor = maximumPosterior / posteriorSum;
        	
        	StringBuffer buf = new StringBuffer();
        	for (int i = 0; i < posteriors.length; i++) {
        		posteriors[i] = posteriors[i] * scaleFactor;
        		if (posteriors[i] < minPrior) {
        			minPrior = posteriors[i];
        		}
        		
        		String strVal = String.format("%1.7f", posteriors[i]);
        		if ((buf.length() == 0) || ((buf.length() + strVal.length() + 1) <= lineLength)) {
        			buf.append(strVal + " ");
        		}
        		else {
        			ps.println(buf.toString());
        			buf = new StringBuffer();
        			buf.append(strVal + " ");
        		}
        	}
        	//print the last buffer
        	ps.println(buf.toString());
        }
        catch (NotFoundException ex) {
			// TODO Auto-generated catch block
			ex.printStackTrace();
		}        
    }

    public void close() {
        try {
            ps.close();
            if (os != null) {
                os.close();                
            }
            ps.close();
            os = null;
            ps = null;
        } catch (IOException ex) {
            throw new RuntimeException(ex.toString(),ex);
        }
    }

    public boolean isClosed() {return ps == null;}
}
