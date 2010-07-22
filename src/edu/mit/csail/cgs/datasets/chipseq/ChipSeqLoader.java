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
import edu.mit.csail.cgs.projects.readdb.SingleHit;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
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

    public List<ChipSeqHit> convert(Collection<SingleHit> input, ChipSeqAlignment align) {
        Genome g = align.getGenome();
        ArrayList<ChipSeqHit> output = new ArrayList<ChipSeqHit>();
        for (SingleHit s : input) {
            int start = s.pos;
            int end = s.strand ? s.pos + s.length : s.pos - s.length;
            output.add(new ChipSeqHit(g, g.getChromName(s.chrom), Math.min(start,end), Math.max(start,end),
                                      s.strand ? '+' : '-', align, s.weight));
        }
        return output;
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


	public ChipSeqAlignment loadAlignment(ChipSeqExpt expt, String n, Genome g) throws NotFoundException, SQLException {
		ChipSeqAlignment align = null;
		PreparedStatement ps = ChipSeqAlignment.createLoadByNameAndExptStatement(cxn);
		ps.setString(1, n);
		ps.setInt(2, expt.getDBID());

		ResultSet rs = ps.executeQuery();        
		while (align == null && rs.next()) {
			align = new ChipSeqAlignment(rs, expt);
            if (!align.getGenome().equals(g)) {
                align = null;
            }
		}
		rs.close();
		ps.close();
        if (align == null) {
            throw new NotFoundException("Couldn't find alignment " + n + " for " + expt + " in genome " + g);
        }
        return align;

	}
	public ChipSeqAlignment loadAlignment(int dbid) throws NotFoundException, SQLException {
		ChipSeqAlignment align = null;
		PreparedStatement ps = ChipSeqAlignment.createLoadByIDStatement(cxn);
		ps.setInt(1, dbid);

		ResultSet rs = ps.executeQuery();
		if (rs.next()) {
			align = new ChipSeqAlignment(rs, this);
		}
		else {
			throw new NotFoundException("Couldn't find alignment by id = " + dbid);
		}
		rs.close();
		ps.close();
		return align;
	}


	public Collection<ChipSeqAlignment> loadAlignments(ChipSeqLocator locator, Genome genome) throws SQLException, NotFoundException {
		List<ChipSeqAlignment> output = new ArrayList<ChipSeqAlignment>();
        for (String rep : locator.getReplicates()) {
            try {
                ChipSeqExpt expt = loadExperiment(locator.getExptName(), rep);
                ChipSeqAlignment align = loadAlignment(expt, locator.getAlignName(), genome);
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
                                 
	public List<ChipSeqHit> loadAllHits(ChipSeqAlignment a) throws IOException {
		List<ChipSeqHit> data = new ArrayList<ChipSeqHit>();
        String alignid = Integer.toString(a.getDBID());
        try {
            for (int chromid : client.getChroms(alignid, false, false)) {
                data.addAll(convert(client.getSingleHits(alignid, chromid,null,null,null,null),a));
            }
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
                            
	public Collection<ChipSeqAlignment> loadAlignments(String name, String replicate, String align,
                                                       Integer factor, Integer cells, Integer condition,
                                                       Genome genome) throws SQLException {
        String query = "select id, expt, name, genome from chipseqalignments";
        if (name != null || replicate != null || align != null || factor != null || cells != null || condition != null || genome != null) {
            query += " where ";
        }
        boolean and = false;
        if (name != null || replicate != null || factor != null || cells != null || condition != null) {
            query += " expt in ( select id from chipseqexpts where ";
            if (name != null) { query += " name = ? "; and = true;}
            if (replicate != null) { query += (and ? " and " : " ") + " replicate = ? "; and = true;}
            if (factor != null) { query += (and ? " and " : " ") + " factor = " + factor; and = true;}
            if (cells != null) { query += (and ? " and " : " ") + " cells = " + cells; and = true;}
            if (condition != null) { query += (and ? " and " : " ") + " condition = " + condition; and = true;}
            query += ")";
            and = true;
        }
        if (genome != null) {query += (and ? " and " : " ") + " genome = " + genome.getDBID(); and = true; }
        if (align != null) {query += (and ? " and " : " ") + " name = ? "; and = true; }

        PreparedStatement ps = cxn.prepareStatement(query);
        int index = 1;
        if (name != null || replicate != null) {
            if (name != null) { ps.setString(index++,name);}
            if (replicate != null) { ps.setString(index++,replicate);}
        }
        if (align != null) {ps.setString(index++,align);}
        
        ResultSet rs = ps.executeQuery();
        Collection<ChipSeqAlignment> output = new ArrayList<ChipSeqAlignment>();
        while (rs.next()) {
            try {
                output.add(new ChipSeqAlignment(rs,this));
            } catch (NotFoundException e) {
                throw new DatabaseException(e.toString(),e);
            }
        }
        rs.close();
        ps.close();
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
                                 
	public List<ChipSeqHit> loadByChrom(ChipSeqAlignment a, int chromid) throws IOException {
		List<ChipSeqHit> data = new ArrayList<ChipSeqHit>();
        String alignid = Integer.toString(a.getDBID());
        try {
            data.addAll(convert(client.getSingleHits(alignid, chromid,null,null,null,null),a));
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
			
	public List<ChipSeqHit> loadByRegion(ChipSeqAlignment align, Region r) throws IOException {
        try {
            return convert(client.getSingleHits(Integer.toString(align.getDBID()),
                                                r.getGenome().getChromID(r.getChrom()),
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                null), align);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
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
    
    /* if Region is a StrandedRegion, then the positions returned are only for that strand */
    public List<Integer> positionsByRegion(List<ChipSeqAlignment> alignments, Region r) throws IOException, ClientException {
		if (alignments.size() < 1) {
			throw new IllegalArgumentException("Alignment List must not be empty.");
		}
        List<Integer> output = new ArrayList<Integer>();
        for (ChipSeqAlignment a : alignments) {
            int[] pos = client.getPositions(Integer.toString(a.getDBID()),
                                            r.getGenome().getChromID(r.getChrom()),
                                            false,
                                            r.getStart(),
                                            r.getEnd(),
                                            null,
                                            null,
                                            r instanceof StrandedRegion ? null : (((StrandedRegion)r).getStrand() == '+'));
            for (int i = 0; i < pos.length; i++) {
                output.add(pos[i]);
            }                                            
        }
        return output;
    }
    public List<Integer> positionsByRegion(ChipSeqAlignment alignment, Region r) throws IOException {
        List<Integer> output = new ArrayList<Integer>();
        try {
            int[] pos = client.getPositions(Integer.toString(alignment.getDBID()),
                                            r.getGenome().getChromID(r.getChrom()),
                                            false,
                                            r.getStart(),
                                            r.getEnd(),
                                            null,
                                            null,
                                            r instanceof StrandedRegion ? null : (((StrandedRegion)r).getStrand() == '+'));
            for (int i = 0; i < pos.length; i++) {
                output.add(pos[i]);
            }                                            
            return output;
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }

	public int countByRegion(ChipSeqAlignment align, Region r) throws IOException {
        try {
            return client.getCount(Integer.toString(align.getDBID()),
                                   r.getGenome().getChromID(r.getChrom()),
                                   false,
                                   r.getStart(),
                                   r.getEnd(),
                                   null,
                                   null,
                                   null);
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
            return client.getCount(Integer.toString(align.getDBID()),
                                   r.getGenome().getChromID(r.getChrom()),
                                   false,
                                   r.getStart(),
                                   r.getEnd(),
                                   null,
                                   null,
                                   r.getStrand() == '+');
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
                total += client.getWeight(Integer.toString(a.getDBID()),
                                          r.getGenome().getChromID(r.getChrom()),
                                          false,
                                          r.getStart(),
                                          r.getEnd(),
                                          null,
                                          null,
                                          null);
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
                total += client.getWeight(Integer.toString(a.getDBID()),
                                          r.getGenome().getChromID(r.getChrom()),
                                          false,
                                          r.getStart(),
                                          r.getEnd(),
                                          null,
                                          null,
                                          r.getStrand() == '+');
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
            return client.getCount(Integer.toString(align.getDBID()),
                                   false,false,false);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}


	public double weighAllHits(ChipSeqAlignment align) throws IOException {
        try {
            return client.getWeight(Integer.toString(align.getDBID()),
                                    false,false,false);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}

    public static Map<String,String> readParameters(BufferedReader reader) throws IOException {
		Map<String, String> params = new HashMap<String, String>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			int p = line.indexOf('=');
			String key = line.substring(0, p);
			String value = line.substring(p + 1);
			params.put(key, value);
		}
		reader.close();
        return params;
    }

	public void addAlignmentParameters(ChipSeqAlignment align, File paramsfile) throws SQLException, IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(paramsfile)));
		addAlignmentParameters(align, readParameters(reader));
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
        try {
            long count = client.getCount(Integer.toString(align.getDBID()), false, false, strand=='+');
            double weight = client.getWeight(Integer.toString(align.getDBID()), false, false, strand=='+');
            Pair<Long,Double> output = new Pair<Long,Double>(count,weight);
            return output;
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /** Generates a histogram of the total weight of reads mapped to each bin.
     * Output maps bin center to weight centered around that bin.  Each read
     * is summarized by its start position.
     */
    public Map<Integer,Float> histogramWeight(ChipSeqAlignment align, char strand, Region r, int binWidth) throws IOException {
        try {
            return client.getWeightHistogram(Integer.toString(align.getDBID()),
                                             r.getGenome().getChromID(r.getChrom()),
                                             false,
                                             false,
                                             binWidth,
                                             r.getStart(),
                                             r.getEnd(),
                                             null,
                                             strand == '+');
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
                                       r.getGenome().getChromID(r.getChrom()),
                                       false,
                                       false,
                                       binWidth,
                                       r.getStart(),
                                       r.getEnd(),
                                       null,
                                       strand == '+');
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
        if (closeMetaLoader && client != null) {
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