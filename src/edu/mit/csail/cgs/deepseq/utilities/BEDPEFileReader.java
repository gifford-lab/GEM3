package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

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
                    if (readOneChr == readTwoChr) {
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
            Read currRead=null;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if(line.charAt(0)!='#'){
                    String[] words = line.split("\\s+");
                    if (words.length<10){
                        System.err.println("Line "+(currID+1));
                        System.err.println(line+"\nBEDPE format should have at least 10 fields!");
                        return;
                    }
                        
                    String chr="."; char strand = '.';
                    int startOne=0, endOne=0;
                    int startTwo=0, endTwo=0;
                    
                    //String ID = words[3]; //No reliable ID for BED format, so treat EVERY hit as a new/unique read

                    if(currRead!=null){
                        currRead.setNumHits(currReadHitCount);
                        //Add the hits to the data structure
                        addHits(currRead);
                        currRead=null;
                    }
                    currReadHitCount=1;                     
                    try{
                        chr = words[0];
                        String[] tmp = chr.split("\\.");
                        chr=tmp[0].replaceFirst("chr", "");
                        chr=chr.replaceFirst("^>", "");
                        // http://genome.ucsc.edu/FAQ/FAQformat.html#format1
                        //BED format is half open - The chromEnd base is not included  
                        // For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
                        startOne = new Integer(words[1]).intValue();
                        endOne = new Integer(words[2]).intValue();
                        if(readLengthOne==-1)
                            readLengthOne = endOne-startOne;
                        strand = words[5].charAt(0);
                        ReadHit currHit = new ReadHit(gen,currID,chr, startOne, endOne-1, strand);
                        currID++;
                        currRead = new Read((int)totalWeight);
                        totalWeight++;
                        currRead.addHit(currHit);
                    } catch (NumberFormatException e){
                        // skip reading this line for header or comment lines
                    }
                }
            }
            if(currRead!=null){
                currRead.setNumHits(currReadHitCount);
                //Add the hits to the data structure
                addHits(currRead);
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
        File testFile = new File("/Users/jennylin/Documents/Jenny/UROP/chr1.bedpe");
        BEDPEFileReader testReading = new BEDPEFileReader(testFile);
    }

}
