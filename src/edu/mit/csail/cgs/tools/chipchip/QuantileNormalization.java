package edu.mit.csail.cgs.tools.chipchip;

import java.sql.*;
import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Quantile normalization takes a set of input experiments and produces a set of output experiments.
 * It works by ranking the probes in the input experiments by their ratio and then computing
 * the average ratio for each rank.  The ratio for a probe in the output experiment is
 * the average ratio for a probe of that rank.  This class does not alter the MOR or channelratio because MOR comes only from
 * the feature-extracted file and channelratio is defined to be the ratio of the channel after normalization
 * that effects the intensities but *before* any normalization that works only on the ratio measurements.
 *
 * QuantileNormalization doesn't extend NormalizeExperiment because that class assumes a single
 * input and single output experiment.
 *
 * The normalization uses two temporary tables:
 * tempA collects all of the data that you're going to normalize.  It's the same as the data table except
 * it adds a rank column that indicates the rank of the ratio within the experiment.
 * tempB maps ranks to the average ratio for that rank.
 *
 * To run this normalization, provide "--version 'new data version'" on the command line and a list of experiments
 * (either as "name;version;replicate" or "name;version") on stdin.  New experiments will be created with the normalization 
 * data and the specified version string.
 *
 * 
 * This class is Oracle-specific.
 *
 * 
 */

public class QuantileNormalization {

    private List<Experiment> inputExpts, outputExpts;
    private Map<Integer,Integer> exptIDMap;
    private String newVersion;
    private ChipChipMetadataLoader loader;
    private MetadataLoader general;
    private String tempA, tempB;
    private java.sql.Connection cxn;

    /* newVersion is the version that the new experimetns
       will have
    */
    public QuantileNormalization(String newVersion) throws SQLException {
        inputExpts = new ArrayList<Experiment>();
        outputExpts = new ArrayList<Experiment>();
        exptIDMap = new HashMap<Integer,Integer>();
        this.newVersion = newVersion;
        loader = new ChipChipMetadataLoader();
        general = new MetadataLoader();
        cxn = DatabaseFactory.getConnection("chipchip");
        System.err.println("New Version is " + newVersion);
    }

    /* fills in inputExpts, outputExpts, and exptIDMap
       lines should be semicolon or tab separated and contain
       either 'name;version' or 'name;version;replicate'
     */
    public void readExptList(BufferedReader reader) throws IOException, NotFoundException, SQLException {
        String line;
        PreparedStatement getToGenome = cxn.prepareStatement("select genome from exptToGenome where experiment = ?");
        PreparedStatement addToGenome = cxn.prepareStatement("insert into exptToGenome(genome,experiment) values(?,?)");
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("[;\\t]");
            ExptNameVersion env;
            if (pieces.length == 2) {
                env = new ExptNameVersion(pieces[0], pieces[1]);
            } else if (pieces.length == 3) {                
                env = new ExptNameVersion(pieces[0], pieces[1], pieces[2]);
            } else {
                throw new IllegalArgumentException("lines must contain two or three tab or semicolon separated fields to specify name, version, and optionally replicate");
            }
            Collection<Experiment> expts = loader.loadExperiment(env);
            for (Experiment oldexpt : expts) {
                FragDist fd = loader.loadFragDist(oldexpt.getFragDist());
                Organism org = new Organism(oldexpt.getSpecies());

                Experiment newexpt = loader.getExperiment(oldexpt.getName(),
                                                          newVersion,
                                                          oldexpt.getReplicate(),
                                                          org.getName(),
                                                          fd.getName(),
                                                          fd.getVersion(),
                                                          general.loadFactor(oldexpt.getFactorOne()).getName(),
                                                          general.loadFactor(oldexpt.getFactorTwo()).getName(),
                                                          general.loadCells(oldexpt.getCellsOne()).getName(),
                                                          general.loadCells(oldexpt.getCellsTwo()).getName(),
                                                          general.loadCondition(oldexpt.getConditionOne()).getName(),
                                                          general.loadCondition(oldexpt.getConditionTwo()).getName(),
                                                          true);
                inputExpts.add(oldexpt);
                outputExpts.add(newexpt);
                exptIDMap.put(oldexpt.getDBID(),
                              newexpt.getDBID());
                System.err.println(String.format("Mapping %d to %d",
                                                 oldexpt.getDBID(),
                                                 newexpt.getDBID()));
                getToGenome.setInt(1,oldexpt.getDBID());
                ResultSet rs = getToGenome.executeQuery();
                while (rs.next()) {
                    addToGenome.setInt(1,rs.getInt(1));
                    addToGenome.setInt(2,newexpt.getDBID());
                    addToGenome.execute();
                }
                rs.close();
            }
        }
        getToGenome.close();
        addToGenome.close();
    }
    
    /* create the names for the two temporary tables.
       The names use the username so they're created in your personal schema
       rather than the chipchip schema.  Thsi way they don't clutter up the
       public areas if they're left around for some reason.
    */
    public void createTempTables() throws SQLException {
        /* don't use the full time because it might exceed the
           number of permissible characters in a table name */
        long time = System.currentTimeMillis() % 100000;
        //        String user = System.getProperty("user.name");
        String user = cxn.getMetaData().getUserName();
        tempA = user + ".tempA" + time;
        tempB = user + ".tempB" + time;
    }
    public void dropTempTables() throws SQLException {
        Statement stmt = cxn.createStatement();
        stmt.execute("drop table " + tempA);
        stmt.execute("drop table " + tempB);
        stmt.close();
    } 

    /* runs the normalization 
     */
    public void runNorm() throws SQLException {
        Statement stmt = cxn.createStatement();
        StringBuffer exptSelection = new StringBuffer();
        for (int i = 0; i < inputExpts.size(); i++) {
            if (i > 0) {
                exptSelection.append(",");
            }
            exptSelection.append(inputExpts.get(i).getDBID());
        }
        System.err.println("Expt string is " + exptSelection.toString());
        stmt.execute("create table " + tempA + 
                     "(experiment, probe, channelone, channeltwo, mor, channelratio, ratio, rank) nologging " + 
                     " as select experiment, probe, channelone, channeltwo, mor, channelratio, ratio, rank() " +
                     " over (partition by experiment order by ratio) from data where experiment in (" +
                     exptSelection.toString() + ")");
        stmt.execute("create index " + tempA + "_rank on " + tempA + "(rank)");
        stmt.execute("create table " + tempB + " (rank, mean)  nologging as select rank, avg(ratio) from " + tempA + " group by rank");
        stmt.execute("create index " + tempB + "_rank on " + tempB + "(rank)");
        stmt.execute("update " + tempA + " set ratio = (select mean from " + tempB +" where " + tempB +".rank = " + tempA + ".rank)");
        for (int oldid : exptIDMap.keySet()) {
            stmt.execute("update " + tempA + " set experiment = " + exptIDMap.get(oldid) + " where experiment = " + oldid);
        }
        stmt.execute("insert into data (experiment, probe, channelone, channeltwo, mor, channelratio, ratio) " + 
                     " (select experiment, probe, channelone, channeltwo, mor, channelratio, ratio from " + tempA + ")");
        stmt.close();
    }

    public static void main(String args[]) throws Exception {
        QuantileNormalization qn = new QuantileNormalization(Args.parseString(args,"version","quantile norm"));
        qn.readExptList(new BufferedReader(new InputStreamReader(System.in)));
        qn.createTempTables();
        qn.runNorm();
        qn.dropTempTables();
    }
}
