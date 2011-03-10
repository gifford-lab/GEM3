package edu.mit.csail.cgs.ewok.verbs;

import java.io.IOException;
import java.io.File;
import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;

/** 
 * <code>SequenceGenerator</code> maps a Region to the genomic
 * sequence included in that Region.
 */
public class SequenceGenerator<X extends Region> implements Mapper<X,String>, SelfDescribingVerb {

    private static Map<Integer,String> cache;
    private boolean useCache = false;
    private boolean useLocalFiles = true;

    private static Map<String, String[]> regionCache;
    private static Map<String, int[]> regionStarts;
    private static boolean regionIsCached = false;
    
    // no longer used, but kept for compatibility 
    public SequenceGenerator (Genome g) {        
    }
    public SequenceGenerator() {}
    public void useCache(boolean b) {
        if (b && cache == null) {
            cache = new HashMap<Integer,String>();
        } 
        useCache = b;
    }
    public void useLocalFiles(boolean b) {
        useLocalFiles = b;
    }
    private void cache(X region) throws SQLException, IOException {
        int chromid = region.getGenome().getChromID(region.getChrom());
        synchronized(cache) {
            if (cache.containsKey(chromid)) {
                return;
            }
        }
        String chromseq = null;
        if (useLocalFiles) {
            File f = new File("/scratch/" + region.getGenome().getVersion() + "/chr" + region.getChrom() + ".fa");
            if (!f.exists()) {
                f = new File("/scratch/" + region.getGenome().getVersion() + "/chr" + region.getChrom() + ".fasta");
            }
            if (f.exists()) {
                FASTAStream stream = new FASTAStream(f);
                while (stream.hasNext()) {
                    Pair<String,String> pair = stream.next();
                    String pairchrom = pair.car().replaceFirst("^chr","");
                    if (pairchrom.equals(region.getChrom())) {
                        chromseq = pair.cdr();
                        break;
                    }                    
                }
                stream.close();                
            }
        }
        if (chromseq == null) {
            java.sql.Connection cxn = DatabaseFactory.getConnection("core");
            PreparedStatement ps = cxn.prepareStatement("select sequence from chromsequence where id = ?");
            ps.setInt(1,chromid);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                chromseq = rs.getString(1);
            }   
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
        }
        if (chromseq == null) {
            return;
        }
        synchronized(cache) {
            if (!cache.containsKey(chromid)) {
                cache.put(chromid, chromseq);
            }
        }
    }
    public String execute(X region) {
        String result = null;
        String chromname = region.getChrom();
        
        try {
            Genome genome = region.getGenome();
            int chromid = genome.getChromID(chromname);
            if (useCache) {
            	if (regionIsCached)
            		return getRegionCacheSequence(region);
            	
                cache(region);
                String chromString = null;
                synchronized(cache) {
                    if (!cache.containsKey(chromid)) {
                        return null;
                    }
                    chromString = cache.get(chromid);
                }
                result = chromString.substring(region.getStart(), region.getEnd() + 1);
            }
            if (result == null) {
                java.sql.Connection cxn =
                DatabaseFactory.getConnection("core");
                PreparedStatement ps;
                int start = Math.max(region.getStart() + 1,0);
                ps = cxn.prepareStatement("select substr(sequence,?,?) from chromsequence where id = ?");
                ps.setInt(1,start);
                ps.setInt(2,region.getEnd() - region.getStart() + 1);
                ps.setInt(3,chromid);                   
                ResultSet rs = ps.executeQuery();
                if (rs.next()) {
                    result = rs.getString(1);
                } 
                rs.close();
                ps.close();
                cxn.commit();
                DatabaseFactory.freeConnection(cxn);
            }
        } catch (SQLException ex) {
            ex.printStackTrace();           
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect to core",ex);
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException("Couldn't load file to cache " + ex.toString(), ex);
        }

        if (result == null) {
            throw new DatabaseException("Couldn't get any sequence for " + region);
        }

        if (result.length() != region.getWidth()) {
            System.err.println("Wanted " + region + "(" + 
            		region.getWidth() + ") but only got " + result.length());
        }

        return result;
    }
    /**
     * Compact the cache of genome sequences to cover only the specified regions
     * @param regions sorted, non-overlapping regions
     */
    public void compactRegionCache(ArrayList<Region> regions){
    	if (regions==null||regions.isEmpty())
    		return;
    	
    	useCache(true);
    	regionCache = new HashMap<String, String[]>();
    	regionStarts = new HashMap<String, int[]>();
    	Genome g = regions.get(0).getGenome();
    	Region lastRegion = regions.get(regions.size()-1);
    	// setup the space
    	String chrom = regions.get(0).getChrom();
    	int count = 0;
    	for(Region r: regions){
    		count ++;
    		if (!r.getChrom().equals(chrom)){	
    			regionCache.put(chrom, new String[count]);
    			regionStarts.put(chrom, new int[count]);
    			chrom = r.getChrom();
    			count = 1;
    		}
    	}
		regionCache.put(chrom, new String[count]);
		regionStarts.put(chrom, new int[count]);
    	
		chrom = regions.get(0).getChrom();
    	count = 0;
    	for (Region r:regions){
    		count ++;
    		synchronized(regionCache) {
    			regionStarts.get(chrom)[count]=r.getStart(); 
        		regionCache.get(chrom)[count]=execute((X)r);    			
    		}
    		if (!r.getChrom().equals(chrom)){	// new Chrom
    			System.out.println("Compact sequence cache: finish Chrom " + chrom);
    			synchronized(cache) {
    				cache.put(g.getChromID(chrom), null);
    				cache.remove(g.getChromID(chrom));	// clean cach for last chrom
    			}
    	    	System.gc();
    			chrom = r.getChrom();
    			count = 1;
    		}
    	}
    	synchronized(cache) {
    		cache.put(g.getChromID(lastRegion.getChrom()), null);
    		cache.remove(g.getChromID(lastRegion.getChrom()));
    	}
    	cache=null;
    	System.gc();
    	
    	regionIsCached = true;
    }
    
    private String getRegionCacheSequence(Region r){
    	int[] starts = regionStarts.get(r.getChrom());
    	int idx = Arrays.binarySearch(starts, r.getStart());
    	if( idx < 0 ) { idx = -idx - 1; }
    	if (!regionCache.containsKey(r.getChrom()))
    		return null;
	    synchronized(regionCache) {
	    	return regionCache.get(r.getChrom())[idx].substring(r.getStart()-starts[idx], r.getEnd()-starts[idx]+1);
    	}
    }
    
    public static void clearCache() {
        synchronized(cache) {
            cache.clear();
        }
    }
    private static final String[] inputNames = { "Regions" };
    private static final EchoType[] inputTypes = { new ClassType(Region.class) };
    private static final EchoType outputType = new ClassType(String.class);

    public EchoType[] getInputClasses() { return inputTypes; }

    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }    
}
