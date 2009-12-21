package edu.mit.csail.cgs.tools.sequence;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Generates two fasta files that contain homologous sequence from two different species.
 *
 * PullPairedFasta --speciesone "$SC;SGDv1" --speciestwo "$SC;Sigmav6" --fastaone s288c.fasta --fastatwo sigma.fasta 
 *   reads list of regions on stdin, expands each by 1000bp (override with, eg, --expand 2000) to get the
 *   source region.  Expands source by --alignexpand (default 5000) to get the query region.  Finds best alignment
 *   for the query region and uses that to get the sequence corresponding to the source region.  Writes
 *   the query and the corresponding region to the FASTA file.
 * The reason for the two expansions is to have alignexpand be fairly large such that you actually find the
 *   right spot in the other genome.  If the query region is fairly small (few hundred by), it may be a good
 *   match to several spots
 * 
 * You can also give --fasta and all sequences will be written to the same fasta.  They'll be interleaved, so it may
 * be harder to process a file with more than one source region (but easier to handle if there's only one source region)
 *
 * You can also omit a --fasta* on the command line; in this case, one file will be created per input region and the
 * file name will be based on the region name
 */


public class PullPairedFASTA {

    public static void main(String args[]) throws Exception {

        int expand = Args.parseInteger(args,"expand",1000);
        int alignExpand = Args.parseInteger(args,"alignexpand",5000);
        String alignType = Args.parseString(args,"aligntype","blast");

        String sone = Args.parseString(args,"speciesone",null);
        String stwo = Args.parseString(args,"speciestwo",null);
        String fone = Args.parseString(args,"fastaone",null);
        String ftwo = Args.parseString(args,"fastatwo",null);
        String fcombined = Args.parseString(args,"fasta",null);

        Genome sourceGenome, targetGenome;

        String pieces[] = sone.split(";");
        sourceGenome = new Genome(pieces[0], pieces[1]);
        pieces = stwo.split(";");
        targetGenome = new Genome(pieces[0], pieces[1]);
        MultiZAlignGenerator generator = new MultiZAlignGenerator(sourceGenome, targetGenome);
        if (alignType != null) {
            generator.setAlignPrefix(alignType);
        }

        FastaWriter<Region> sourceFasta = fone != null ? new FastaWriter<Region>(fone) : null;
        FastaWriter<Region> targetFasta = ftwo != null ? new FastaWriter<Region>(ftwo) : null;
        FastaWriter<Region> combinedFasta = fcombined != null ? new FastaWriter<Region>(fcombined) : null;
        
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            try {
                Region source = Region.fromString(sourceGenome,line).expand(expand,expand);
                Region query = source.expand(alignExpand,alignExpand);
                Iterator<MultiZAlignRegion> iter = generator.execute(query);
                ArrayList<MultiZAlignRegion> blocks = new ArrayList<MultiZAlignRegion>();
                while (iter.hasNext()) {
                    blocks.add(iter.next());
                }
                AlignmentStitcher<MultiZAlignRegion> stitcher = new AlignmentStitcher<MultiZAlignRegion>(blocks);
                Region target = stitcher.getBestAlignedRegion(source);
            
                if (target == null || source == null) {
                    continue;
                }
                if (target.getWidth() > 10 * source.getWidth()) {
                    continue;
                }



                if (sourceFasta != null && targetFasta != null) {
                    sourceFasta.consume(source);
                    targetFasta.consume(target);
                }
                if (combinedFasta != null) {
                    combinedFasta.consume(source);
                    combinedFasta.consume(target);
                }
                if (combinedFasta == null && (sourceFasta == null || targetFasta == null)) {
                    FastaWriter<Region> thisFasta = new FastaWriter<Region>(source.toString() + ".fasta");
                    thisFasta.consume(source);
                    thisFasta.consume(target);
                    thisFasta.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }




        }
        if (sourceFasta != null) {
            sourceFasta.close();
        }
        if (targetFasta != null) {
            targetFasta.close();       
        }
        if (combinedFasta != null) {
            combinedFasta.close();
        }


    }

}