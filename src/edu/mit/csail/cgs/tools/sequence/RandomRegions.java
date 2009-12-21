package edu.mit.csail.cgs.tools.sequence;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.RepeatMaskedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.RepeatMaskedGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Generate, on stdout, a list of random regions in a genome.
 *
 */

public class RandomRegions {
	
	public static void main(String[] args) throws NotFoundException{
        int numSamples = Args.parseInteger(args,"numSamples",5000);
        int validSamples=0;
        int regionSize = Args.parseInteger(args,"regionSize",200);
        Genome gen = Args.parseGenome(args).cdr();
        RepeatMaskedGenerator repMask = null; 
        if (Args.parseFlags(args).contains("repeatMask")) {
			repMask = new RepeatMaskedGenerator(gen);
        }        
    
		double genomeSize=0;
		long [] chromoSize;
		String [] chromoNames;
		int numChroms=0;
		ArrayList<Region> regList = new ArrayList<Region>();
		Random rand = new Random();
		
        //First see how big the genome is:
        chromoSize = new long[gen.getChromList().size()];
        chromoNames = new String[gen.getChromList().size()];
        Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
        while (chroms.hasNext()) {
            NamedRegion currentChrom = chroms.next();
            genomeSize += (double)currentChrom.getWidth();
            chromoSize[numChroms]=currentChrom.getWidth();
            chromoNames[numChroms]=currentChrom.getChrom();
            //System.out.println(chromoNames[numChroms]+"\t"+chromoSize[numChroms]);
            numChroms++;				
        }//System.out.println(genomeSize);
			
        //Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
        while(validSamples<numSamples){
            Region potential;				
            long randPos = (long)(1+(rand.nextDouble()*genomeSize));
            //find the chr
            boolean found=false;
            long total=0;
            for(int c=0; c<numChroms && !found; c++){
                if(randPos<total+chromoSize[c]){
                    found=true;
                    if(randPos+regionSize<total+chromoSize[c]){
                        potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+regionSize-total));
							
                        //is this overlapping a repeat?
                        boolean repOver=false;
                        if (repMask != null) {
                            Iterator<RepeatMaskedRegion> repItr = repMask.execute(potential);
                            while(repItr.hasNext()){
                                RepeatMaskedRegion currRep = repItr.next();
                                if(currRep.overlaps(potential)){
                                    repOver=false;
                                }
                            }
                        }
							
                        if(!repOver){
                            validSamples++;
                            regList.add(potential);
                            System.out.println(potential.getChrom()+":"+potential.getStart()+"-"+potential.getEnd());
                        }
                    }
                }
                total+=chromoSize[c];
            }				
        }						
	}

}
