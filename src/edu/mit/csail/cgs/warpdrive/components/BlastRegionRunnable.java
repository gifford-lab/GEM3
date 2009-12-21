/*
 * Created on Mar 2, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;


import javax.swing.JFrame;

import org.xml.sax.*;
import org.xml.sax.helpers.XMLReaderFactory;

import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.warpdrive.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class BlastRegionRunnable implements Runnable {
    private Genome base, target;
    private BLASTInterface blastInterface;    
    private Region region;
    private RegionList regionList;
        
    public BlastRegionRunnable(Genome baseGenome, Genome targetGenome, BLASTInterface inter) {
        base = baseGenome;
        target = targetGenome;
        blastInterface = inter;
        region = null;
    }

    public void run() { 
        SequenceGenerator reg2seq = new SequenceGenerator();
        String seq = reg2seq.execute(region);

        try {
            String queryResponse = blastInterface.doQuery(seq);

            BLASTHitHandler handler = new BLASTHitHandler();

            XMLReader oParser = XMLReaderFactory.createXMLReader(); 
            oParser.setContentHandler(handler);
            EntityResolver resolver = new EntityResolver() {
                    public InputSource resolveEntity(String publicId, String systemId) {
                        
                        System.out.println("Resolve: \"" + publicId + "\", \"" + systemId + "\"");                       

                        String dtdPath = "edu.mit.csail.cgs/dtd/";
                        String[] recs = { "NCBI_BlastOutput.dtd", "NCBI_BlastOutput.mod", "NCBI_Entity.mod" };
                        String[] paths = { "NCBI_BlastOutput.dtd", 
                                           "NCBI_BlastOutput.mod", "NCBI_Entity.mod" };
                        

                        for(int i = 0; i < recs.length; i++) { 
                            if(systemId.endsWith(recs[i])) { 
                                URL url = ClassLoader.getSystemClassLoader().getResource(dtdPath + paths[i]);
                                try {
                                    return new InputSource(url.openStream());
                                } catch (IOException ex) {
                                    System.err.println("Couldn't open " + url);
                                    ex.printStackTrace();
                                }
                            }
                        }

                        
                        System.out.println("No match.");
                        return null;
                    }
                };
            oParser.setEntityResolver(resolver);

            StringReader reader = new StringReader(queryResponse);
            System.err.println("Parsing");
            try { 
                oParser.parse(new InputSource(reader));
            } catch(MalformedURLException e) { 
                System.err.println("Message: \"" + e.getMessage() + "\"");
                e.printStackTrace(System.err);
            }
            System.err.println("extracting hit regions");
            LinkedList<ScoredStrandedRegion> hits = handler.extractHitRegions(target);
            for(ScoredStrandedRegion hit : hits) {       
                System.err.println("Adding " + hit);
                regionList.addRegion(hit);
            }
            System.err.println("Done");
            
            //             RegionFrame targetFrame = null;
            //             if(RegionFrame.genomeFrameMap.containsKey(target.getVersion())) { 
            //                 targetFrame = RegionFrame.genomeFrameMap.get(target.getVersion());
            //             } else { 
            //                 WarpOptions opts = new WarpOptions(target.getVersion());

            //                 Region startRegion = null;
            //                 if(!hits.isEmpty()) {
            //                     startRegion= hits.getFirst();
            //                 } else { 
            //                     List<String> chroms = target.getChromList();
            //                     String firstChrom = chroms.get(0);
            //                     int start = 1;
            //                     int end = Math.min(10000, target.getChromLength(firstChrom));
            //                     startRegion = new Region(target, firstChrom, start, end);
            //                 }

            //                 opts.chrom = startRegion.getChrom();
            //                 opts.start = startRegion.getStart();
            //                 opts.stop = startRegion.getEnd();
                
            //                 targetFrame = new RegionFrame(opts);
            //             }
            
            //             RegionPanel targetPanel = targetFrame.getRegionPanel();
            //             RegionListPanel listPanel = new RegionListPanel(targetPanel, hits);
            //             JFrame listFrame = RegionListPanel.makeFrame(listPanel, "BLAST Hits");

        } catch (IOException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            e.printStackTrace();
        }
    }
    
    public void startInThread(Region r, RegionList rl) { 
        region = r;
        regionList = rl;
        Thread t = new Thread(this);
        t.start();
    }    
}
