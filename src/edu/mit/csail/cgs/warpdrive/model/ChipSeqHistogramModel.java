package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * Data model for chipseq histogram.  Separate methods for retrieving
 * plus and minus strand results
 */
public class ChipSeqHistogramModel extends WarpModel implements RegionModel, Runnable {
    
    private Client client;
    private TreeMap<Integer,Float> resultsPlus, resultsMinus;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    private ChipSeqHistogramProperties props;

    private Region region;
    private boolean newinput;

    public ChipSeqHistogramModel (ChipSeqAlignment a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.add(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        ids = new HashSet<String>();
        ids.add(Integer.toString(a.getDBID()));
    }
    public ChipSeqHistogramModel (Collection<ChipSeqAlignment> a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.addAll(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        ids = new HashSet<String>();
        for (ChipSeqAlignment align : alignments) {
            ids.add(Integer.toString(align.getDBID()));
        }
    }    
    public ChipSeqHistogramProperties getProperties() {return props;}
    
    
    public void clearValues() {
        resultsPlus = null;
        resultsMinus = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public Map<Integer,Float> getPlus() {return resultsPlus;}
    public Map<Integer,Float> getMinus() {return resultsMinus;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    int width = props.BinWidth;
                    boolean extension = props.ReadExtension;
                    // for GaussianKernel, get 1bp resolution data
                    if (props.GaussianKernelWidth!=0 && region.getWidth()<=1000){ 
                    	width = 1;
                    }
                    resultsPlus = null;
                    resultsMinus = null;
                    if (props.UseWeights) {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                            resultsPlus = Aggregator.mergeHistogramsFF(resultsPlus,
                                                                       client.getWeightHistogram(ids,
                                                                                                 region.getGenome().getChromID(region.getChrom()),
                                                                                                 false,
                                                                                                 extension,
                                                                                                 width,
                                                                                                 (int)props.DeDuplicate,
                                                                                                 region.getStart(),
                                                                                                 region.getEnd(),
                                                                                                 null,
                                                                                                 true));
                            
                            resultsMinus = Aggregator.mergeHistogramsFF(resultsMinus,
                                                                        client.getWeightHistogram(ids,
                                                                                                  region.getGenome().getChromID(region.getChrom()),
                                                                                                  false,
                                                                                                  extension,
                                                                                                  width,
                                                                                                  (int)props.DeDuplicate,
                                                                                                  region.getStart(),
                                                                                                  region.getEnd(),
                                                                                                  null,
                                                                                                  false));
                        }
                        if (props.ShowPairedReads) {
                            resultsPlus = Aggregator.mergeHistogramsFF(resultsPlus,
                                                                    client.getWeightHistogram(ids,
                                                                                              region.getGenome().getChromID(region.getChrom()),
                                                                                              true,
                                                                                              extension,
                                                                                              width,
                                                                                              (int)props.DeDuplicate,
                                                                                              region.getStart(),
                                                                                              region.getEnd(),
                                                                                              null,
                                                                                              true));
                            
                            resultsMinus = Aggregator.mergeHistogramsFF(resultsMinus,
                                                                     client.getWeightHistogram(ids,
                                                                                               region.getGenome().getChromID(region.getChrom()),
                                                                                               true,
                                                                                               extension,
                                                                                               width,
                                                                                               (int)props.DeDuplicate,
                                                                                               region.getStart(),
                                                                                               region.getEnd(),
                                                                                               null,
                                                                                               false));

                        }

                    } else {
                        if (!props.ShowPairedReads || props.ShowSingleReads) {
                            resultsPlus = Aggregator.mergeHistogramsIF(client.getHistogram(ids,
                                                                                           region.getGenome().getChromID(region.getChrom()),
                                                                                           false,
                                                                                           extension,
                                                                                           width,
                                                                                           (int)props.DeDuplicate,
                                                                                           region.getStart(),
                                                                                           region.getEnd(),
                                                                                           null,
                                                                                           true),
                                                                    resultsPlus);
                            
                            resultsMinus = Aggregator.mergeHistogramsIF(client.getHistogram(ids,
                                                                                            region.getGenome().getChromID(region.getChrom()),
                                                                                            false,
                                                                                            extension,
                                                                                            width,
                                                                                            (int)props.DeDuplicate,
                                                                                            region.getStart(),
                                                                                            region.getEnd(),
                                                                                            null,
                                                                                            false),
                                                                     resultsMinus);
                        }
                        if (props.ShowPairedReads) {
                            resultsPlus = Aggregator.mergeHistogramsIF(
                                                                    client.getHistogram(ids,
                                                                                        region.getGenome().getChromID(region.getChrom()),
                                                                                        true,
                                                                                        extension,
                                                                                        width,
                                                                                        (int)props.DeDuplicate,
                                                                                        region.getStart(),
                                                                                        region.getEnd(),
                                                                                        null,
                                                                                        true),
                                                                       resultsPlus);
                            
                            resultsMinus = Aggregator.mergeHistogramsIF(
                                                                     client.getHistogram(ids,
                                                                                         region.getGenome().getChromID(region.getChrom()),
                                                                                         true,
                                                                                         extension,
                                                                                         width,
                                                                                         (int)props.DeDuplicate,
                                                                                         region.getStart(),
                                                                                         region.getEnd(),
                                                                                         null,
                                                                                         false),
                                                                       resultsMinus);

                        }

                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    // assign empty output.  This is useful because Client
                    // throws an exception for non-existant chromosomes, such
                    // as those for which there were no alignment results
                    resultsPlus = new TreeMap<Integer,Float>();
                    resultsMinus = resultsPlus;
                }
                newinput = false;
                notifyListeners();
            }
        }
        System.err.println("ChipSeqHistogram Model is closing");
        client.close();
    }                     
 }