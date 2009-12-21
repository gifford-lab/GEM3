package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.datasets.general.Region;

/** Given a set of aligned regions between two genomes, stitch them together to
 *  get the optimal overall alignment between the two genomes.  The individual
 *  aligned regions come from Blast or Blat or similar; any part of the
 *  query may occur in any number of the aligned regions.  Any part
 *  of the query will occur in only one of the output blocks.
 *
 *  AlignmentStitcher uses a two-phase DP algorithm to find the subset
 *  of the input blocks that represent the best overall alignment of
 *  the query to the target.
 */

public class AlignmentStitcher<MZAR extends MultiZAlignRegion> {
    
    public static final int ENDBLOCK = Integer.MAX_VALUE;

    /* in figuring out the best set of alignment blocks to use,
       AlignmentStitcher tries to maximize the sum of the
       scores of the blocks it keeps (subject obviously to conditions like
       the blocks can't overlap in the query genome). The constants
       below define score penalties in that maximization.

       CHROMOSOMESWITCHPENALTY is assessed when the query region
       aligns to two target chromosomes.  The penalty is high
       as this is assumed to be an unlikely event.

       SPLITPENALTY is for reordering blocks or otherwise
       deviating from a single-chromosome in order
       best path.
    */
    public static final int DEFAULTCHROMOSOMESWITCHPENALTY = 80000;    
    public static final int DEFAULTSPLITPENALTY = 500;
    private int chromosomeSwitchPenalty = DEFAULTCHROMOSOMESWITCHPENALTY;
    private int splitPenalty = DEFAULTSPLITPENALTY;
    Map<String,MZAR[]> byChrom;
    String bestChrom;
    Map<String,List<MZAR>> bestAlign;

    /* these are reset every time you call stitch */
    private int n, size;
    private int[][] inorderBlocks, inorderstartswith, inorderendswith;
    private double[][] inorderScores;
    private int[] blocks, startswith, endswith;
    private double[] score;
    private int[][] originalIndex;
    private Class mzarclass, mzararrayclass;

    @SuppressWarnings({"unchecked"})
    public AlignmentStitcher(Collection<MZAR> input) {
        byChrom = new HashMap<String,MZAR[]>();
        Map<String,List<MZAR>> tbc = new HashMap<String,List<MZAR>>();
        mzarclass = null;
        mzararrayclass = null;
        for (MZAR r : input) {
            if (!tbc.containsKey(r.getChrom())) {
                tbc.put(r.getChrom(),new ArrayList<MZAR>());
            }
            tbc.get(r.getChrom()).add(r);
            if (mzarclass == null) {
                mzarclass = r.getClass();                
            }

        }        
        for (String c : tbc.keySet()) {
            List<MZAR> l = tbc.get(c);
            Collections.sort(l);
            MZAR[] a = (MZAR[]) java.lang.reflect.Array.newInstance(mzarclass, l.size());
            for (int i = 0; i < l.size(); i++) {
                a[i] = l.get(i);
            }
            byChrom.put(c, a);
            if (mzararrayclass == null) {
                mzararrayclass = a.getClass();
            }

        }
        bestAlign = new HashMap<String,List<MZAR>>();
        size = -1;
        stitch();
    }
    public void splitPenalty(int s) {splitPenalty = s;}
    public void chromSwitchPenalty(int c) {chromosomeSwitchPenalty = c;}
    public void stitch() {
        double bestChromScore = 0;
        bestChrom = null;

        long totalcomb = 0;
        long scoretest = 0;
        long overlap = 0; 
        
        for (String sourceChrom: byChrom.keySet()) {
            MZAR[] regions = byChrom.get(sourceChrom);            

            int chromcount = 0;
            HashMap<String,Integer> targetChroms = new HashMap<String,Integer>();
            Map<Integer,Integer> chromCounts = new HashMap<Integer,Integer>();
            for (MZAR mzar : regions) {
                String targetChrom = mzar.getOtherChrom();
                if (!targetChroms.containsKey(targetChrom)) {
                    chromCounts.put(chromcount,0);
                    targetChroms.put(targetChrom,
                                     chromcount++);
                }
                int c = targetChroms.get(targetChrom);
                chromCounts.put(c, chromCounts.get(c) + 1);
            }            

            MZAR[][] byTargetChrom = (MZAR[][])java.lang.reflect.Array.newInstance(mzararrayclass, chromcount);
            originalIndex = new int[chromcount][];
            int seen[] = new int[chromcount];
            for (int i = 0; i < chromcount; i++) {
                byTargetChrom[i] = (MZAR[])java.lang.reflect.Array.newInstance(mzarclass, chromCounts.get(i));
                originalIndex[i] = new int[chromCounts.get(i)];
                seen[i] = 0;
            }
            int count = 0;
            for (MZAR mzar : regions) {
                String targetChrom = mzar.getOtherChrom();
                int chromindex = targetChroms.get(targetChrom);
                byTargetChrom[chromindex][seen[chromindex]] = mzar;
                originalIndex[chromindex][seen[chromindex]] = count;
                seen[chromindex]++;
                count++;
            }
            int numchroms = targetChroms.size();
            /* first step is to compute the best in-order path from i to j. 
               This allows skipping aligned blocks between i and j, but no 
               reordering in the source genome or target genome.  Also, each
               path uses exactly one chromosome in each genome.

               The score of a path is just the sum of the scores of each included block; there
               are no gap penalties or anything like that

               To make this algorithm work, we run it on each target chromosome individually and then
               merge the results back into the full set of results
            */
            inorderBlocks = new int[numchroms][];
            inorderstartswith = new int[numchroms][];
            inorderendswith = new int[numchroms][];
            inorderScores = new double[numchroms][];

            count = 0;
            //            System.err.println("BLOCKS ARE ");
            for (MZAR mzar : regions) {
                //                System.err.println("\t" + count + ": " + mzar + "    " + (int)mzar.getScore());
                count++;
            }
            System.err.println("Source chrom " + sourceChrom + " at " + System.currentTimeMillis());
            for (String targetChrom : targetChroms.keySet()) {                
                // at this point we're working with one source and one target chromosome
                int thischromindex = targetChroms.get(targetChrom);
                MZAR[] localregions = byTargetChrom[thischromindex];
                int[] originalIndices = originalIndex[thischromindex];
                int localsize = localregions.length;
                int localn = (localsize * localsize + localsize) / 2;
                int[] localInOrderBlocks = new int[localn];
                int[] localStartsWith = new int[localn];
                int[] localEndsWith = new int[localn];
                double[] localInOrderScores = new double[localn];
                inorderBlocks[thischromindex] = localInOrderBlocks;
                inorderstartswith[thischromindex] = localStartsWith;
                inorderendswith[thischromindex] = localEndsWith;
                inorderScores[thischromindex] = localInOrderScores;
                for (int diff = 0; diff <= localsize; diff++) {
                    for (int start = 0; start <= localsize - 1 - diff; start++) {
                        int stop = start + diff;
                        int thisindex = getIndex(start,stop);
                        if (start == stop) {
                            localInOrderBlocks[thisindex] = ENDBLOCK;
                            localInOrderScores[thisindex] = localregions[start].getScore();
                            localStartsWith[thisindex] = start;
                            localEndsWith[thisindex] = start;
                        } else {
                            /* Find best straight path: either includes the current block (regions[start]) or not */
                            double bestscorewith = 0, bestscorewithout = 0;
                            int beststartwith = -1, beststartwithout = -1;
                            int nextblockwith = -1, nextblockwithout = -1;
                        
                            int se = localregions[start].getEnd();
                            String c = localregions[start].getOtherChrom();

                            int withoutindex = getIndex(start+1, stop);
                            bestscorewithout = localInOrderScores[withoutindex];
                            beststartwithout = localStartsWith[withoutindex];
                            nextblockwithout = localEndsWith[withoutindex];

                            bestscorewith = localregions[start].getScore();
                            beststartwith = start;
                            nextblockwith = ENDBLOCK;  

                            for (int otherstart = start+1; otherstart <= stop; otherstart++) {
                                int otherindex = getIndex(otherstart,stop);
                                MZAR otherstartswith = localregions[localStartsWith[otherindex]];                                                                
                                if (otherstartswith.getStart() > se) {
                                    double scoreWith = localregions[start].getScore() + localInOrderScores[otherindex];
                                    bestscorewith = scoreWith;
                                    beststartwith = otherstart;
                                    nextblockwith = otherstart;
                                    break;
                                } 
                            }
                            if (bestscorewith <= 0 && bestscorewithout <= 0) {
                                throw new RuntimeException("start=" + start + " stop=" + stop + "  and couldn't find anything");
                            }                       
                            if (bestscorewith > bestscorewithout) {
                                if (nextblockwith < 0) {
                                    throw new RuntimeException("Didn't initialize nextblockwith for " + start + "," + stop);
                                }                            
                                localInOrderScores[thisindex] = bestscorewith;
                                localInOrderBlocks[thisindex] = nextblockwith;
                                localStartsWith[thisindex] = start;
                                localEndsWith[thisindex] = nextblockwith == ENDBLOCK ? start : localEndsWith[getIndex(beststartwith,stop)];
                            } else {
                                if (nextblockwithout < 0) {
                                    throw new RuntimeException("Didn't initialize nextblockwithout for " + start + "," + stop);
                                }
                                localInOrderScores[thisindex] = bestscorewithout;
                                localInOrderBlocks[thisindex] = -1 * nextblockwithout;
                                localStartsWith[thisindex] = localStartsWith[getIndex(beststartwithout,stop)];
                                localEndsWith[thisindex] = localEndsWith[getIndex(beststartwithout,stop)];
                            }
                        }
                    }
                }
            }
            System.err.println("Finished first pass at " + System.currentTimeMillis());
            int size = regions.length;
            int n = (size * size + size) / 2;
            /* block contains, for each i,j, the list of blocks
               (indices into regions) that is the optimal path through
               this source chromosome.  It can reference blocks on
               different target chromosomes.
            */
            blocks = new int[n];
            startswith = new int[n];
            endswith = new int[n];
            score = new double[n];
            for (int i = 0; i < n; i++) {
                score[i] = 0;
            }
            int localStarts[][] = new int[numchroms][size];
            int localStops[][] = new int[numchroms][size];
            
            for (int c = 0; c < numchroms; c++) {
                int[] orig = originalIndex[c];
                int s = orig.length;
                for (int i = 0; i < size; i++) {
                    int startindex = Arrays.binarySearch(orig,i);
                    if (startindex < 0) {
                        startindex = -1 - startindex;
                    }
                    int stopindex = Arrays.binarySearch(orig,i);
                    if (stopindex < 0) {
                        stopindex = -2 - stopindex;
                    }
                    if (startindex >= s) {
                        startindex = s - 1;
                    }
                    if (stopindex >= s) {
                        stopindex = s -1;
                    }
                    localStarts[c][i] = startindex;
                    localStops[c][i] = stopindex;
                }
            }


            System.err.println("Starting second pass with size " + size);
            for (int diff = 0; diff <= size; diff++) {

                //                 if (diff % 200 == 0) {
                //                     System.err.println("allorder " + diff);
                //                 }

                for (int start = 0; start <= size - 1 - diff; start++) {
                    int stop = start + diff;
                    int index = getIndex(start,stop);
                    if (start == stop) {
                        score[index] = regions[start].getScore();
                        blocks[index] = ENDBLOCK;
                        startswith[index] = start;
                        endswith[index] = start;
                        continue;
                    }
                    double straightscore = 0; 
                    int straightchrom = -1, straightstartswith = -1, straightendswith = -1;
                    for (int c = 0; c < numchroms; c++) {
                        int startindex = localStarts[c][start];
                        int stopindex = localStops[c][stop];
                        if (startindex > stopindex) {
                            continue;
                        }
                        int localindex = getIndex(startindex, stopindex);
                        if (inorderScores[c][localindex] > straightscore) {
                            straightscore = inorderScores[c][localindex];
                            straightchrom = c;
                            straightstartswith = originalIndex[c][inorderstartswith[c][localindex]];
                            straightendswith = originalIndex[c][inorderendswith[c][localindex]];
                        }
                    }

                    double bestscore = 0;
                    int bestk = -1;
                    for (int k = start; k < stop; k++) {
                        totalcomb++;
                        int indexleft = getIndex(start,k);
                        int indexright = getIndex(k+1,stop);

                        MZAR l = regions[endswith[indexleft]];
                        MZAR r = regions[startswith[indexright]];

                        if (l.overlaps(r)) {
                            overlap++;
                            continue;
                            //                            splitscore -= l.getScore() + r.getScore();
                        }

                        double splitscore = score[indexleft] + score[indexright] - splitPenalty;
                        if (splitscore < bestscore) {
                            scoretest++;
                            continue;
                        }
                        if (!l.getOtherChrom().equals(r.getOtherChrom())) {
                            splitscore -= chromosomeSwitchPenalty;
                        }
                        if (splitscore > bestscore) {
                            bestscore = splitscore;
                            bestk = k;
                        }
                    }

                    //                     if (start == 0 && stop >= (size - 2)) {
                    //                         System.err.println("start=" + start + "  stop=" + stop + "  straightscore= " + straightscore + 
                    //                                            "bestscore = "+bestscore + " split at " + bestk);
                    //                     }

                    if (bestscore > straightscore && bestk >= 0) {
                        score[index] = bestscore;
                        blocks[index] = bestk;
                        int indexleft = getIndex(start,bestk);
                        int indexright = getIndex(bestk+1,stop);
                        startswith[index] = startswith[indexleft];
                        endswith[index] = endswith[indexright];
                    } else {
                        score[index] = straightscore;
                        blocks[index] = -1 * (straightchrom + 1);
                        startswith[index] = straightstartswith;
                        endswith[index] = straightendswith;
                    }
                }                                
            }
            System.err.println("Finished second pass at " + System.currentTimeMillis());
            
            ArrayList<MZAR> bestPath = new ArrayList<MZAR>();
            List<Integer> bestblockindices = getBestAlignment(0,size-1);
            //            System.err.println("INDICES ARE " + bestblockindices);
            for (int i : bestblockindices) {
                bestPath.add(regions[i]);
            }
            count = 0;
            //            System.err.println("PATH IS ");
            for (MZAR mzar : bestPath) {
                //                System.err.println("\t" + count + ": " + mzar);
                count++;
            }
            int wholeindex = getIndex(0,size-1);
            bestAlign.put(sourceChrom,bestPath);
            if (score[wholeindex] > bestChromScore) {
                bestChromScore = score[wholeindex];
                bestChrom = sourceChrom;
            }
            //            System.err.println("STITCHER : " + sourceChrom + "," + score[wholeindex] + "," + bestPath);
        }
        //        System.err.println(" Aligned to " + bestChrom);
        System.err.println(String.format("Stats total %d, score %d, overlap %d", totalcomb, scoretest, overlap));
    }

    private void getInOrderBlocks(int start,
                                  int end,
                                  int c,
                                  List<Integer> out) {
        int index = getIndex(start,end);
        if (inorderstartswith[c][index] != start) {
            getInOrderBlocks(inorderstartswith[c][index],
                             end,
                             c,
                             out);
        } else {
            int next = inorderBlocks[c][index];
            if (next > 0) {
                out.add(originalIndex[c][start]);
            } else {
                next *= -1;
            }
            if (next != ENDBLOCK) {
                getInOrderBlocks(next,
                                 end,
                                 c,
                                 out);
            }
        }
    }
    
    private List<Integer> getBestAlignment(int start,
                                           int end) {
        int index = getIndex(start,end);
        if (blocks[index] < 0) {
            int c = -1 * blocks[index] - 1;
            List<Integer> out = new ArrayList<Integer>();
            Pair<Integer,Integer> pair = getLocalIndices(start,
                                                         end,
                                                         originalIndex[c]);
            getInOrderBlocks(pair.car(),
                             pair.cdr(),
                             c,
                             out);
            //            System.err.println("STRAIGHT from " + start + " to " + end + ": " + out +  "   score=" + inorderScores[c][getIndex(pair.car(), pair.cdr())]);
            return out;
        } else if (blocks[index] == ENDBLOCK) {
            //            System.err.println("END at " + start + " to " + end);
            List<Integer> out = new ArrayList<Integer>();
            out.add(start);
            return out;
        } else {
            //            System.err.println("SPLIT " + start + " to " + end + " at " + blocks[index]);
            List<Integer> left = getBestAlignment(start,
                                                  blocks[index]);
            left.addAll(getBestAlignment(blocks[index] + 1,
                                         end));
            return left;
        }
    }

    private Pair<Integer,Integer> getLocalIndices(int start,
                                                  int stop,
                                                  int[] originalIndex) {
        int startindex = Arrays.binarySearch(originalIndex,start);
        if (startindex < 0) {
            startindex = -1 - startindex;
        }
        int stopindex = Arrays.binarySearch(originalIndex,stop);
        if (stopindex < 0) {
            stopindex = -2 - stopindex;
        }
        int s = originalIndex.length;
        if (startindex >= s) {
            startindex = s - 1;
        }
        if (stopindex >= s) {
            stopindex = s -1;
        }
        return new Pair<Integer,Integer>(startindex,stopindex);
    }

    public String getBestAlignedChrom() {
        return bestChrom;
    }
    /* gets the best alignment with the specified source chromosome */
    public List<MZAR> getBestAlignment(String c) {
        if (bestAlign.get(c) == null) {
            //            System.err.println("RETURNING EMPTY");
            return new ArrayList<MZAR>();
        } else {
            return bestAlign.get(c);
        }
    }
    /* maps a numeric coordinate to the coresponding coordinate in the other genome.
       This assumes that hte list of MZARs are a coherent alignment along a single 
       chromosome.
    */
    public Pair<Integer, String> mapCoord(List<MZAR> blocks, int sourceCoord) {
        if (blocks.size() == 0) {
            return null;
        }

        if (sourceCoord < blocks.get(0).getStart()) {
            return new Pair<Integer,String>(blocks.get(0).getStart(),
                                            blocks.get(0).getOtherChrom());
        }
        if (sourceCoord > blocks.get(blocks.size()-1).getEnd()) {
            return new Pair<Integer,String>(blocks.get(blocks.size()-1).getEnd(),
                                            blocks.get(blocks.size()-1).getOtherChrom());
        }
        for (int i = 0; i < blocks.size(); i++) {
            MZAR block = blocks.get(i);            
            if (sourceCoord >= block.getStart() && sourceCoord <= block.getEnd()) {                
                double factor = (sourceCoord - block.getStart()) / ((double)block.getWidth());
                if (block.getStrand() == '+') {
                    int coord = block.getOtherStart() + (int)(block.getOtherWidth() * factor);                   
                    return new Pair<Integer,String>(coord, block.getOtherChrom());
                } else {
                    return new Pair<Integer,String>((int)(block.getOtherEnd() - block.getOtherWidth() * factor),
                                                    block.getOtherChrom());
                }
            }
            if (i < blocks.size() - 1) {
                MZAR next = blocks.get(i+1);
                if (!next.getOtherChrom().equals(block.getOtherChrom())) {
                    /* alignment shifts chromosomes; no good way to do this,
                       so throw an error */
                    throw new IllegalArgumentException("Trying to map a coordinate between chromosomes");
                }
                if (sourceCoord >= block.getEnd() && sourceCoord <= next.getStart()) {
                    float w = next.getStart() - block.getEnd();
                    float otherw = next.getOtherStart() - block.getOtherEnd();
                    if (block.getOtherEnd() < next.getOtherStart()) {
                        return new Pair<Integer,String>((int)(block.getOtherEnd() + otherw * (sourceCoord - block.getEnd()) / w),
                                                        next.getOtherChrom());
                    } else {
                        return new Pair<Integer,String>((int)(block.getOtherStart() - otherw * (sourceCoord - block.getEnd()) / w),
                                                        next.getOtherChrom());   
                    }
                }
            }
        }        
        throw new RuntimeException("Couldn't figure out where " + sourceCoord + " goes");
    }
    /* tries to map the specified input region into a corresponding region in the other
       genome (as returned by getBestAlignment) by interpolationg the input region's
       boundaries within alignment blocks.  This will throw an IllegalArgumentException
       if the input region spans a chromosome boundary in the alignment
    */
    public Region getBestAlignedRegion(Region input) {
        List<MZAR> blocks = getBestAlignment(getBestAlignedChrom());
        Pair<Integer,String> start = mapCoord(blocks,input.getStart());
        Pair<Integer,String> stop = mapCoord(blocks,input.getEnd());
        if (start == null || stop == null) {
            return null;
        }

        if (!stop.cdr().equals(start.cdr())) {
            throw new IllegalArgumentException("Trying to map a region that spans target chromosomes");
        }
        int startpos = start.car();
        int stoppos = stop.car();
        if (stoppos < startpos) {
            startpos = stop.car();
            stoppos = start.car();
        }

        Region output = new Region(blocks.get(0).getOtherGenome(),
                                   blocks.get(0).getOtherChrom(),
                                   startpos,
                                   stoppos);
        return output;
        
    }

    private int getIndex(int start, int stop) {
        if (start > stop) {
            throw new IllegalArgumentException("start > stop");
        }
        int rowlen = (stop * stop + stop) / 2;
        int index = rowlen + start;
        return index;

    }
}