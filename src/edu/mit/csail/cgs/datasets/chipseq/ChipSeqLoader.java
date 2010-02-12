/*
 * Created on May 16, 2007
 */
package edu.mit.csail.cgs.datasets.chipseq;

import java.util.*;
import java.sql.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.stats.StatUtil;

import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;

/**
 * @author tdanford
 */
public class ChipSeqLoader implements edu.mit.csail.cgs.utils.Closeable {

	public static String role = "chipseq";


	public static void main(String[] args) throws Exception{
		try {
			ChipSeqLoader loader = new ChipSeqLoader();
			Collection<ChipSeqExpt> expts = loader.loadAllExperiments();
			for (ChipSeqExpt expt : expts) {
				Collection<ChipSeqAlignment> aligns = loader.loadAllAlignments(expt);
				for (ChipSeqAlignment align : aligns) {
					System.out.println(expt.getDBID() + "\t" + expt.getName() + ";"+ expt.getReplicate()+"\t"+align.getName()+"\t"+align.getDBID()+"\t"+align.getGenome());
				}				
			}
		}
		catch (SQLException e) {
			e.printStackTrace();
		}
	}

	private MetadataLoader metaLoader;
	private boolean closeMetaLoader;
	private java.sql.Connection cxn;
    private Client client=null;
    
    public ChipSeqLoader() throws SQLException, IOException{this(true);}
	public ChipSeqLoader(boolean openClient) throws SQLException, IOException {
		metaLoader = new MetadataLoader();
		closeMetaLoader = true;
		if(openClient){
	        try {
	            client = new Client();
	        } catch (ClientException e) {
	            throw new IllegalArgumentException(e);
	        }
		}
		cxn = DatabaseFactory.getConnection(role);
	}


	public MetadataLoader getMetadataLoader() {
		return metaLoader;
	}


	public Collection<Genome> loadExperimentGenomes(ChipSeqExpt expt) throws SQLException {
		LinkedList<Genome> genomes = new LinkedList<Genome>();
		String query = String.format("select genome from chipseqalignments where expt=%d", expt.getDBID());
		Statement s = cxn.createStatement();
		ResultSet rs = s.executeQuery(query);
		while (rs.next()) {
			int gid = rs.getInt(1);
			try {
				Genome g = Organism.findGenome(gid);
				genomes.add(g);
			}
			catch (NotFoundException e) {
				e.printStackTrace();
			}
		}
		rs.close();
		s.close();
		return genomes;
	}


	public Collection<ChipSeqExpt> loadAllExperiments() throws SQLException {
		PreparedStatement ps = ChipSeqExpt.createLoadAll(cxn);
		LinkedList<ChipSeqExpt> expts = new LinkedList<ChipSeqExpt>();
		ResultSet rs = ps.executeQuery();
		while (rs.next()) {
			expts.addLast(new ChipSeqExpt(rs, this));
		}
		rs.close();
		ps.close();

		return expts;
	}


	public ChipSeqExpt loadExperiment(String name, String rep) throws NotFoundException, SQLException {
		PreparedStatement ps = ChipSeqExpt.createLoadByNameReplicate(cxn);
		ps.setString(1, name);
		ps.setString(2, rep);
		ResultSet rs = ps.executeQuery();
		ChipSeqExpt expt = null;
		if (rs.next()) {
			expt = new ChipSeqExpt(rs, this);
		}
		rs.close();
		ps.close();

		if (expt == null) { throw new NotFoundException(name); }
		return expt;
	}


	public Collection<ChipSeqExpt> loadExperiments(String name) throws SQLException {
		PreparedStatement ps = ChipSeqExpt.createLoadByName(cxn);
		LinkedList<ChipSeqExpt> expts = new LinkedList<ChipSeqExpt>();
		ps.setString(1, name);
		ResultSet rs = ps.executeQuery();
		ChipSeqExpt expt = null;
		while (rs.next()) {
			expt = new ChipSeqExpt(rs, this);
			expts.add(expt);
		}
		rs.close();
		ps.close();

		return expts;
	}


	public ChipSeqExpt loadExperiment(int dbid) throws NotFoundException, SQLException {
		PreparedStatement ps = ChipSeqExpt.createLoadByDBID(cxn);
		ps.setInt(1, dbid);
		ResultSet rs = ps.executeQuery();
		ChipSeqExpt expt = null;
		if (rs.next()) {
			expt = new ChipSeqExpt(rs, this);
		}
		rs.close();
		ps.close();

		if (expt == null) {
			String err = String.format("No such ChipPet Experiment %d", dbid);
			throw new NotFoundException(err);
		}
		return expt;
	}

    public Collection<ChipSeqAlignment> loadAlignments (Genome g) throws SQLException {
		Collection<ChipSeqAlignment> aligns = new LinkedList<ChipSeqAlignment>();
		PreparedStatement ps = ChipSeqAlignment.createLoadAllByGenomeStatement(cxn);
		ps.setInt(1, g.getDBID());
        ResultSet rs = ps.executeQuery();
		while (rs.next()) {
            try {
                ChipSeqAlignment align = new ChipSeqAlignment(rs, loadExperiment(rs.getInt(2)));
                aligns.add(align);
            } catch (NotFoundException e) {
                // this only happens if we get back an invalid expt ID, which shouldn't happen.
                e.printStackTrace();
                throw new DatabaseException(e.toString());
            }
		}
		rs.close();
		ps.close();
		return aligns;
    }

	public Collection<ChipSeqAlignment> loadAllAlignments(ChipSeqExpt expt) throws SQLException {
		Collection<ChipSeqAlignment> aligns = new LinkedList<ChipSeqAlignment>();
		PreparedStatement ps = ChipSeqAlignment.createLoadAllByExptStatement(cxn);
		ps.setInt(1, expt.getDBID());

		ResultSet rs = ps.executeQuery();
		while (rs.next()) {
			ChipSeqAlignment align = new ChipSeqAlignment(rs, expt);
			aligns.add(align);
		}
		rs.close();

		ps.close();
		return aligns;
	}


	public ChipSeqAlignment loadAlignment(ChipSeqExpt expt, String n) throws NotFoundException, SQLException {
		ChipSeqAlignment align = null;
		PreparedStatement ps = ChipSeqAlignment.createLoadByNameAndExptStatement(cxn);
		ps.setString(1, n);
		ps.setInt(2, expt.getDBID());

		ResultSet rs = ps.executeQuery();
		if (rs.next()) {
			align = new ChipSeqAlignment(rs, expt);
		}
		else {
			throw new NotFoundException("Couldn't find alignment " + n + " for " + expt);
		}

		rs.close();

		ps.close();
		return align;
	}


	public Collection<ChipSeqAlignment> loadAlignments(ChipSeqLocator locator) throws SQLException, NotFoundException {
		List<ChipSeqAlignment> output = new ArrayList<ChipSeqAlignment>();
		for (String rep : locator.getReplicates()) {
			try {
				ChipSeqExpt expt = loadExperiment(locator.getExptName(), rep);
				ChipSeqAlignment align = loadAlignment(expt, locator.getAlignName());
				if (align != null) {
					output.add(align);
				}
			}
			catch (IllegalArgumentException e) {
				throw new NotFoundException("Couldn't find experiment for " + locator);
			}
		}
		return output;
	}

    private void instantiateHits(Collection<ChipSeqHit> output,
                                 int[] positions,
                                 float[] weights,
                                 Genome g,
                                 String chrom,
                                 char strand,
                                 ChipSeqAlignment align) {
        int readlen = align.getExpt().getReadLength();
        if (strand == '+') {
            for (int i = 0; i < positions.length; i++) {
                output.add(new ChipSeqHit(align.getGenome(),
                                          chrom,
                                          positions[i],
                                          positions[i] + readlen,
                                          strand,
                                          align,
                                          weights[i]));        
            }
        } else {
            for (int i = 0; i < positions.length; i++) {
                output.add(new ChipSeqHit(align.getGenome(),
                                          chrom,
                                          positions[i] - readlen,
                                          positions[i],
                                          strand,
                                          align,
                                          weights[i]));        
            }
        }
    }
                                 
	public Vector<ChipSeqHit> loadAllHits(ChipSeqAlignment a) throws IOException {
		Vector<ChipSeqHit> data = new Vector<ChipSeqHit>();
        try {
            for (String chrom : client.getChroms(Integer.toString(a.getDBID()))) {
                char strand = chrom.charAt(chrom.length() - 1);
                String justchrom = chrom.substring(0,chrom.length() - 1);
                int[] positions = client.getHits(Integer.toString(a.getDBID()),
                                                 chrom);
                float[] weights = client.getWeights(Integer.toString(a.getDBID()),
                                                    chrom);            
                instantiateHits(data, positions, weights, a.getGenome(), justchrom, strand, a);
            }
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
    
	public Vector<ChipSeqHit> loadByChrom(ChipSeqAlignment a, String chrom) throws IOException {
		Vector<ChipSeqHit> data = new Vector<ChipSeqHit>();
        try {
            char strand = chrom.charAt(chrom.length() - 1);
            String justchrom = chrom.substring(0,chrom.length() - 1);
            int[] positions = client.getHits(Integer.toString(a.getDBID()),
                                             chrom);
            float[] weights = client.getWeights(Integer.toString(a.getDBID()),
                                                chrom);            
            instantiateHits(data, positions, weights, a.getGenome(), justchrom, strand, a);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
	
	public int[] getStartsByChrom(ChipSeqAlignment a, String chrom){
		int[] data = null;
        try {
            int[] positions = client.getHits(Integer.toString(a.getDBID()),
                                             chrom);
            float[] weights = client.getWeights(Integer.toString(a.getDBID()),
                                                chrom);    
            // assume weights are integers
            int totalCount = 0;
            for (float w:weights){
            	totalCount += (int)w;
            }
            data = new int[totalCount];
            for (int i=0; i<positions.length;i++){
            	for (int j=0; j<(int)weights[i];j++){
                	data[i+j]= positions[i];            		
            	}
            }            
        } catch (Exception e) {
            System.err.println(e.toString());
        }
        
		return data;
	}
		
	public Collection<ChipSeqHit> loadByRegion(ChipSeqAlignment align, Region r) throws IOException {
		Genome g = r.getGenome();
		ArrayList<ChipSeqHit> data = new ArrayList<ChipSeqHit>();
        try {
            int[] positions = client.getHitsRange(Integer.toString(align.getDBID()),
                                                  r.getChrom() + '+',
                                                  r.getStart(),
                                                  r.getEnd());
            float[] weights = client.getWeightsRange(Integer.toString(align.getDBID()),
                                                     r.getChrom() + '+',
                                                     r.getStart(),
                                                     r.getEnd());
            instantiateHits(data, positions,weights, g, r.getChrom(), '+', align);
            positions = client.getHitsRange(Integer.toString(align.getDBID()),
                                            r.getChrom() + '-',
                                            r.getStart(),
                                            r.getEnd());
            weights = client.getWeightsRange(Integer.toString(align.getDBID()),
                                             r.getChrom() + '-',
                                             r.getStart(),
                                             r.getEnd());
            instantiateHits(data, positions,weights, g, r.getChrom(), '-', align);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
        return data;
	}


	public Collection<ChipSeqHit> loadByRegion(List<ChipSeqAlignment> alignments, Region r) throws IOException {
		if (alignments.size() < 1) {
			throw new IllegalArgumentException("Alignment List must not be empty.");
		}
        Collection<ChipSeqHit> output = null;
        for (ChipSeqAlignment a : alignments) {
            if (output == null) {
                output = loadByRegion(a,r);
            } else {
                output.addAll(loadByRegion(a,r));
            }
        }
		return output;
	}

	public int countByRegion(ChipSeqAlignment align, Region r) throws IOException {
        try {
            return client.getCountRange(Integer.toString(align.getDBID()),
                                        r.getChrom()+'+',
                                        r.getStart(),
                                        r.getEnd()) + 
                client.getCountRange(Integer.toString(align.getDBID()),
                                     r.getChrom()+'-',
                                     r.getStart(),
                                     r.getEnd());
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    

	public int countByRegion(List<ChipSeqAlignment> alignments, Region r) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        int total = 0;
        for (ChipSeqAlignment a : alignments) {
            total += countByRegion(a,r);
        }
        return total;
	}
	public int countByRegion(ChipSeqAlignment align, StrandedRegion r) throws IOException {
        try {
            return client.getCountRange(Integer.toString(align.getDBID()),
                                        r.getChrom()+r.getStrand(),
                                        r.getStart(),
                                        r.getEnd());
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}
	public int countByRegion(List<ChipSeqAlignment> alignments, StrandedRegion r) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        int total = 0;
        for (ChipSeqAlignment a : alignments) {
            total += countByRegion(a,r);
        }
        return total;
	}

	
	public double weightByRegion(List<ChipSeqAlignment> alignments, Region r) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        double total = 0;
        for (ChipSeqAlignment a : alignments) {
            try {
                total += client.getWeightRange(Integer.toString(a.getDBID()),
                                               r.getChrom() + "+",
                                               r.getStart(),
                                               r.getEnd());
                total += client.getWeightRange(Integer.toString(a.getDBID()),
                                               r.getChrom() + "-",
                                               r.getStart(),
                                               r.getEnd());
            } catch (ClientException e) {
                throw new IllegalArgumentException(e);
            }            
        }
        return total;
	}
	public double weightByRegion(List<ChipSeqAlignment> alignments, StrandedRegion r) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        double total = 0;
        for (ChipSeqAlignment a : alignments) {
            try {
                total += client.getWeightRange(Integer.toString(a.getDBID()),
                                               r.getChrom() + r.getStrand(),
                                               r.getStart(),
                                               r.getEnd());
            } catch (ClientException e) {
                throw new IllegalArgumentException(e);
            }            
        }
        return total;
	}


	/**
	 * @param align
	 * @return
	 * @throws SQLException
	 */
	public int countAllHits(ChipSeqAlignment align) throws IOException {
        try {
            return client.getCount(Integer.toString(align.getDBID()));
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}


	public double weighAllHits(ChipSeqAlignment align) throws IOException {
        try {
            return client.getWeight(Integer.toString(align.getDBID()));
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}


	public void addAlignmentParameters(ChipSeqAlignment align, File paramsfile) throws SQLException, IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(paramsfile)));
		Map<String, String> params = new HashMap<String, String>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			int p = line.indexOf('=');
			String key = line.substring(0, p);
			String value = line.substring(p + 1);
			params.put(key, value);
		}
		reader.close();
		addAlignmentParameters(align, params);
	}


	public void addAlignmentParameters(ChipSeqAlignment align, Map<String, ? extends Object> params) throws SQLException {
		PreparedStatement insert = cxn.prepareStatement("insert into alignmentparameters(alignment,name,value) values(?,?,?)");
		insert.setInt(1, align.getDBID());
		for (String k : params.keySet()) {
			insert.setString(2, k);
			Object val = params.get(k);
			if (val == null) {
				val = "";
			} else {
				val = val.toString();
				if (val == null) {
					val = "";
				}
			}


			insert.setString(3, (String)val);
			insert.execute();
		}
		insert.close();
	}


	public Map<String, String> getAlignmentParameters(ChipSeqAlignment align) throws SQLException {
		Statement get = cxn.createStatement();
		ResultSet rs = get.executeQuery("select name, value from alignmentparameters where alignment = " + align.getDBID());
		Map<String, String> output = new HashMap<String, String>();
		while (rs.next()) {
			output.put(rs.getString(1), rs.getString(2));
		}
		rs.close();
		get.close();
		return output;

	}

    /**
     * Get the total # of hits and weight for an alignment but only include reads
     * on the specified strand.  
     */
    public Pair<Long,Double> getAlignmentStrandedCountWeight(ChipSeqAlignment align, char strand) throws IOException {
        long count = 0;
        double sum = 0;
        try {
            for (String chrom : client.getChroms(Integer.toString(align.getDBID()))) {
                if (chrom.charAt(chrom.length() - 1) != strand) {
                    continue;
                }
                count += client.getCount(Integer.toString(align.getDBID()), chrom);
                sum += client.getWeight(Integer.toString(align.getDBID()), chrom);
            }
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
        Pair<Long,Double> output = new Pair<Long,Double>(count,sum);
        return output;
    }

    /** Generates a histogram of the total weight of reads mapped to each bin.
     * Output maps bin center to weight centered around that bin.  Each read
     * is summarized by its start position.
     */
    public Map<Integer,Float> histogramWeight(ChipSeqAlignment align, char strand, Region r, int binWidth) throws IOException {
        try {
            return client.getWeightHistogram(Integer.toString(align.getDBID()),
                                             r.getChrom() + strand,
                                             r.getStart(),
                                             r.getEnd(),
                                             binWidth);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    /** Generates a histogram of the total number of reads mapped to each bin.
     * Output maps bin center to weight centered around that bin.  Each read
     * is summarized by its start position.
     */
    public Map<Integer,Integer> histogramCount(ChipSeqAlignment align, char strand, Region r, int binWidth) throws IOException {        
        try {
            return client.getHistogram(Integer.toString(align.getDBID()),
                                       r.getChrom() + strand,
                                       r.getStart(),
                                       r.getEnd(),
                                       binWidth);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    /** Generates a probability density of all reads in given region.
     *  Each read is summarized by its start position.
     */
    public double[] kernalDensityInRegion(ChipSeqAlignment align, char strand, Region r, int std) throws SQLException {
        
        PreparedStatement stmt = cxn.prepareStatement("select startpos as pos, sum(weight) from chipseqhits " +
        	"where alignment = ? and chromosome = ? and startpos >= ? and startpos < ? and strand = ? group by startpos");
        stmt.setInt(1,align.getDBID());
        stmt.setInt(2, r.getGenome().getChromID(r.getChrom()));
        stmt.setInt(3, r.getStart());
        stmt.setInt(4, r.getEnd());
        stmt.setString(5, Character.toString(strand));
        double[] readCount = new double[r.getWidth()];
        ResultSet rs = stmt.executeQuery();
        int start=r.getStart();
        while (rs.next()) {
        	readCount[rs.getInt(1)-start] = rs.getDouble(2) ;
        }
        rs.close();
        stmt.close();
        
        return StatUtil.gaussianSmoother(readCount, std);
    }

	public void close() {
		if (closeMetaLoader && !metaLoader.isClosed()) {
			metaLoader.close();
            metaLoader = null;
		}
        if (client != null) {
            client.close();
            client = null;
        }
		DatabaseFactory.freeConnection(cxn);
		cxn = null;
	}


	public boolean isClosed() {
		return cxn == null;
	}

}