package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

//for testing
import java.util.List;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.Pair;

public class BEDPEFileReader extends PairedAlignmentFileReader {

    public BEDPEFileReader(File f) {
        super(f);
    }

    public BEDPEFileReader(File f, Genome g, int mis, boolean nonUnique,
            int idSeed, HashMap<String, Integer> chrom2id,
            HashMap<Integer, String> id2Chrom) {
        super(f, g, mis, nonUnique, idSeed, chrom2id, id2Chrom);
    }

    @Override
    protected void estimateGenome(File f) {
        HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(f));
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if(line.charAt(0)!='#'){
                    String[] words = line.split("\\s+");
                    String readOneChr = words[0];
                    String[] tmp = readOneChr.split("\\.");
                    readOneChr=tmp[0].replaceFirst("chr", "");
                    readOneChr=readOneChr.replaceFirst("^>", "");
                    String readTwoChr = words[3];
                    tmp = readTwoChr.split("\\.");
                    readTwoChr=tmp[0].replaceFirst("chr", "");
                    readTwoChr=readOneChr.replaceFirst("^>", "");
                    if (readOneChr.equals(readTwoChr)) {
                        int readOneMax = Math.max(new Integer(words[1]).intValue(), new Integer(words[2]).intValue());
                        int readTwoMax = Math.max(new Integer(words[4]).intValue(), new Integer(words[5]).intValue());
                        int max = Math.max(readOneMax, readTwoMax);
                        
                        if(!chrLenMap.containsKey(readOneChr) || chrLenMap.get(readOneChr)<max)
                            chrLenMap.put(readOneChr, max+1);
                    }
                }
            }
            gen=new Genome("Genome", chrLenMap);
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (NumberFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
    @Override
    //Return the total reads and weight
    protected void countReads() {
        try {
            readLengthOne=-1;
            readLengthTwo=-1;
            insertLength=-1;
            totalHits=0;
            totalWeight=0;
            BufferedReader reader = new BufferedReader(new FileReader(inFile));
            String line, lastID="";
            double currReadHitCount=0;
            Read currReadOne=null;
            Read currReadTwo=null;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if(line.charAt(0)!='#'){
                    String[] words = line.split("\\s+");
                    if (words.length<10){
                        System.err.println("Line "+(currID+1));
                        System.err.println(line+"\nBEDPE format should have at least 10 fields!");
                        return;
                    }
                        
                    String chrOne="."; char strandOne = '.';
                    String chrTwo="."; char strandTwo = '.';
                    int startOne=0, endOne=0;
                    int startTwo=0, endTwo=0;
                    
                    //String ID = words[3]; //No reliable ID for BED format, so treat EVERY hit as a new/unique read

                    if(currReadOne!=null){
                        currReadOne.setNumHits(currReadHitCount);
                        //Add the hits to the data structure
                        addHits(currReadOne, currReadTwo);
                        currReadOne=null;
                    }
                    currReadHitCount=1;                     
                    try{
                        chrOne = words[0];
                        String[] tmp = chrOne.split("\\.");
                        chrOne=tmp[0].replaceFirst("chr", "");
                        chrOne=chrOne.replaceFirst("^>", "");
                        
                        chrTwo = words[3];
                        String[] tmp2 = chrTwo.split("\\.");
                        chrTwo=tmp2[0].replaceFirst("chr", "");
                        chrTwo=chrTwo.replaceFirst("^>", "");
                        // http://genome.ucsc.edu/FAQ/FAQformat.html#format1
                        //BED format is half open - The chromEnd base is not included  
                        // For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
                        startOne = new Integer(words[1]).intValue();
                        endOne = new Integer(words[2]).intValue();
                        
                        startTwo = new Integer(words[4]).intValue();
                        endTwo = new Integer(words[5]).intValue();
                        
                        if(readLengthOne==-1)
                            readLengthOne = endOne-startOne;
                        if(readLengthTwo==-1)
                            readLengthTwo = endTwo-startTwo;
                        if(insertLength==-1)
                            insertLength = startTwo-endOne;
                        strandOne = words[8].charAt(0);
                        strandTwo = words[9].charAt(0);
                        if (!chrOne.matches(chrTwo)||strandOne==strandTwo) {//wrong chr or improper strand pairing
                            //System.out.println("invalid pair");
                            continue;
                        }
                        ReadHit currHitOne = new ReadHit(gen,currID,chrOne, startOne, endOne-1, strandOne);
                        currID++;
                        ReadHit currHitTwo = new ReadHit(gen, currID, chrTwo, startTwo, startTwo-1, strandTwo);
                        currReadOne = new Read((int)totalWeight);
                        currReadTwo = new Read((int)totalWeight);
                        totalWeight++;
                        currReadOne.addHit(currHitOne);
                        currReadTwo.addHit(currHitTwo);
                    } catch (NumberFormatException e){
                        // skip reading this line for header or comment lines
                    }
                }
            }
            if(currReadOne!=null){
                currReadOne.setNumHits(currReadHitCount);
                //Add the hits to the data structure
                addHits(currReadOne, currReadTwo);
            }
            reader.close();
            populateArrays();
            
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }//end of countReads method

    /**
     * @param args
     */
    public static void main(String[] args) {
        File notPEfile = new File("/Users/jennylin/Documents/Jenny/UROP/chr1.bed");
        File testFile = new File("/Users/jennylin/Documents/Jenny/UROP/chr1.bedpe");
        HashMap<String, Integer> c2i = new HashMap<String, Integer>();
        HashMap<Integer, String> i2c = new HashMap<Integer, String>();
        System.out.println("Hi world once");
        BEDPEFileReader testReading = new BEDPEFileReader(testFile);
        System.out.println("Hi world");
        Genome g = testReading.getGenome();
        List<String> chromList = g.getChromList();
        int i=0; 
        for(String c:chromList){
            c2i.put(c, i);
            i2c.put(i, c);
            i++;
        }
        Region testRegion = new Region(g, chromList.get(0), 51338000, 51340000);
        List<File> fileList = new ArrayList<File>();
        fileList.add(notPEfile);
        
        //FileReadLoader frl = new FileReadLoader(fileList, "BED");
        //List<ReadHit> notPEreads = frl.loadHits(testRegion);
        
        //BEDFileReader frl = new BEDFileReader(notPEfile, g, false, 1, c2i, i2c);
        //List<ReadHit> notPEreads = frl.loadHits(testRegion);
        
        testReading = new BEDPEFileReader(testFile, g, 5, false, 1, c2i, i2c);
        List<ReadHit> testReads = testReading.loadHits(testRegion);
        Pair<Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>,ArrayList<ArrayList<Float>>> matrix = 
                testReading.loadStrandedFivePrimeCounts(testRegion, '+');
        ArrayList<Integer> fiveEnd = matrix.getFirst().getFirst();
        ArrayList<ArrayList<Integer>> threeEndMatrix = matrix.getFirst().getLast();
        ArrayList<ArrayList<Float>> weightMatrix = matrix.getLast();
        
        StringBuilder sb = new StringBuilder();
        //StringBuilder hitsb = new StringBuilder();
        for (int a = 0; a<fiveEnd.size(); a++) {
            
            //hitsb.append(fiveEnd.get(a));
            //hitsb.append('\t');
            for (int b = 0; b<threeEndMatrix.get(a).size(); b++) {
                sb.append(fiveEnd.get(a));
                sb.append('\t');
                sb.append(threeEndMatrix.get(a).get(b));
                sb.append('\t');
                sb.append(weightMatrix.get(a).get(b));
                sb.append('\n');
                //hitsb.append('\n');
//                if (b<threeEndMatrix.get(a).size()-1) {
//                    sb.append('\t');
//                    hitsb.append('\t');
//                }
//                else {
//                    sb.append('\n');
//                    hitsb.append('\n');
//                }
            }
        }
        PrintWriter reads;
        //PrintWriter weights;
        try {
            reads = new PrintWriter("/Users/jennylin/Documents/Jenny/UROP/bindResults2.csv");
            reads.write(sb.toString());
            reads.close();
            //weights = new PrintWriter("/Users/jennylin/Documents/Jenny/UROP/weightVect.csv");
            //weights.write(hitsb.toString());
            //weights.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        System.out.println("Done");
    }

}
