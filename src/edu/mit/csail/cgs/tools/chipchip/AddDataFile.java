package edu.mit.csail.cgs.tools.chipchip;

import java.io.*;
import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.parsing.textfiles.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * AddDataFile loads one or more files of microarray data.  The design for the arrays must
 * have been previously loaded.
 *
 * Mandatory arguments
 * --expt 'exptname;version;replicate'
 * --designname array_design_name
 * --species "species;genome"
 * --fragdistname 'name' --fragdistversion 'version'
 * --conditionone YPD --conditiontwo YPD --cellsone BY4741 --cellstwo BY4741
 * --factorone H3K4me3 --factortwo WCE
 *
 * Optional Arguments
 * --median : do median normalization.  sets ip = ip * (median(wce) / median(ip))
 * --mean : 
 * --linefit : do linefitting normalization
 * --fournorms : do crosstalk, median, and linefitting normalization and store the results after each step.  Uses the provided
 *      version and appends "", "crosstalk norm", "median norm", and "median linefit"
 * --controlexpt 'cname;cver;crep' : does control experiment normalization of the ratios
 * --iptowce : sets the coefficient of IP to WCE crosstalk for the crosstalk normalization.  Crosstalk norm is done if 
 *   --iptowce and --wcetoip are set.
 * --wcetoip : sets the WCE to IP coefficient.
 * --scatchTablespace : Oracle tablespace to use for temporary tables.
 * 
 * AddDataFile uses some oracle-specific SQL.  Not only do you need write permissions in the chipchip schema,
 * you need to be able to create tables in your own schema as well.
 */

public class AddDataFile {

    private String factorone, factortwo, cellsone, cellstwo, conditionone, conditiontwo;
    private String name, version, replicate, designname, species, fragdistname, fragdistversion;
    private boolean mean, median, linefit, nofar, noafe, verbose, fournorms, channelswap, bgsub;
    private int blockoffset, designid;
    private String normalization;
    private List<String> ipfiles, wcefiles, bothfiles;
    private java.sql.Connection core, chipchip;
    private ChipChipMetadataLoader loader;
    private Experiment expt;
    private Experiment noNormExpt, crosstalkNormExpt, medianNormExpt, medianLinefitNormExpt;
    private List<ExptNameVersion> controlexpt;
    private List<Region> medianNormRegions, medianLFNormRegions;
    private String dataTable, galTable, galIndex;
    private PreparedStatement insertStmt;
    private String galFileOverride, scratchTablespace;
    private double crossIPtoWCE = Double.NaN, crossWCEtoIP = Double.NaN;
    
    public static void main(String args[]) {
        try {
            AddDataFile add = new AddDataFile();
            if (args.length == 0) {
                add.argsGUI();
            } else {
                add.parseArgs(args);
                add.go();
            }
        } catch (Exception e) {
            System.err.println("Error during setup " + e.toString());
            System.err.println("Came from ");
            e.printStackTrace();
        }
    }

    public void argsGUI() throws SQLException {
        //        JFrame frame = new JFrame("Add Datafiles");
        
    }

    /* main method after command line args have been parsed */
    public void go() throws SQLException, NotFoundException {
        if (fournorms) {
            noNormExpt = getExperiment(("no norm " + version).trim());
            crosstalkNormExpt = getExperiment(("crosstalk norm " + version).trim());
            medianNormExpt = getExperiment(("median " + version).trim());
            medianLinefitNormExpt = getExperiment(("median linefit " + version).trim());
            mapToGenomes(noNormExpt.getDBID());
            mapToGenomes(medianNormExpt.getDBID());
            mapToGenomes(medianLinefitNormExpt.getDBID());
        } else if (version != null && !version.equals("")) {
            expt = getExperiment(version);            
            mapToGenomes(expt.getDBID());                   
        } else {
            throw new IllegalArgumentException("No Version provided and can't used defaults");
        }


        /* bothfiles is a list of filenames that contain data for both channels.
           Compare to ipfiles and wcefiles, which are paired lists that contain filenames where each
           file only contains data for one channel of the experiment.
        */
        try {
            createTempTables();
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        }

        for (String fname : bothfiles) {
            try {
                int galfileid = addFile(fname);    
                assignProbeIDs(galfileid,true);
            } catch (Exception e) {
                System.err.println("Error Adding file " + fname);
                e.printStackTrace();
            }
        }
        for (int i = 0; i < ipfiles.size(); i++) {
            try {
                int galfileid = addFiles(ipfiles.get(i),
                                         wcefiles.get(i));
                assignProbeIDs(galfileid,false);
            } catch (Exception e) {
                System.err.println("Error adding " + ipfiles.get(i));
                e.printStackTrace();
            } 
        }
        try {
            Statement stmt = chipchip.createStatement();
            if (channelswap) {
                System.err.println("Swapping Channels");
                System.err.println("NOTE: this happens *before* the crosstalk normalization, so I'm swapping the constants you gave");
                double t = crossIPtoWCE;
                crossIPtoWCE = crossWCEtoIP;
                crossWCEtoIP = t;
                stmt.execute("update " + dataTable + " set channelone = channeltwo, channeltwo = channelone, mor = 1.0 / mor");
            }
            stmt.execute("update " + dataTable + " set channelratio = channelone / channeltwo, ratio = channelone / channeltwo");
            stmt.close();
            /* run through the normalization steps, saving the intermediate results if 
               requested */
            if (noNormExpt != null) {
                System.err.println("Saving no-norm");
                moveFromTemp(noNormExpt.getDBID());
            }
            crossTalkNormalization();
            if (crosstalkNormExpt != null) {
                System.err.println("Saving crosstalk");
                moveFromTemp(crosstalkNormExpt.getDBID());
            }
            medianOrMeanNorm();
            if (medianNormExpt != null) {
                System.err.println("saving median");
                moveFromTemp(medianNormExpt.getDBID());
            }
            LinefitNorm();
            if (medianLinefitNormExpt != null) {
                System.err.println("saving median linefit");
                moveFromTemp(medianLinefitNormExpt.getDBID());
            }
            if (expt != null) {
                System.err.println("Saving final");
                controlExptNormalize();
                moveFromTemp(expt.getDBID());
            }
        } catch (Exception e) {
            System.err.println("Error doing normalizations");
            e.printStackTrace();
        }
        try {
            chipchip.commit();
        } catch (SQLException e) {
            System.err.println("Couldn't commit at end!");
            e.printStackTrace();
        } 
        try {
            dropTempTables();
        } catch (SQLException e) {
            System.err.println("Couldn't drop temp tables " + e.toString());
            e.printStackTrace();
        }
    }

    public AddDataFile () throws UnknownRoleException, SQLException {
        core = DatabaseFactory.getConnection("core");
        chipchip = DatabaseFactory.getConnection("chipchip");
        chipchip.setAutoCommit(false);
        loader = new ChipChipMetadataLoader();
        /* don't use the full time because it might exceed the
           number of permissible characters in a table name */
        long time = System.currentTimeMillis() % 100000;
        //        String user = System.getProperty("user.name");
        String user = chipchip.getMetaData().getUserName();
        dataTable = user + ".dtmp" + time;
        galTable = user + ".gtmp" + time;
        galIndex = user + ".ix_gtmp" + time;
    }

    /* 
       TODO this SQL is oracle specific
    */
    public void createTempTables() throws SQLException {
        String createDataSQL = "create table " + dataTable + "(" + 
            "experiment number(10) not null," +
            "id number(10)," +
            "probeid varchar(200),"+
            "blockno number(10),"+
            "rowno number(10),"+
            "colno number(10),"+
            "channelone binary_float,"+
            "channeltwo binary_float,"+
            "mor binary_float,"+
            "channelratio binary_float,"+
            "ratio binary_float,"+
            "controlratio binary_float) tablespace " + scratchTablespace + " nologging";
        String createGalSQL = "create table " + galTable + "(" + 
            "id number(10),"+
			"probeid varchar2(200),"+
			"arraydesign number(10),"+
			"galfile number(10),"+
			"blockno number(10),"+
			"colno number(10),"+
			"rowno number(10)"+
            ") tablespace " + scratchTablespace + " nologging";
        Statement stmt = chipchip.createStatement();
        System.err.println("Trying to create " + dataTable + ", " + galTable);
        stmt.execute(createDataSQL);
        stmt.execute(createGalSQL);

        insertStmt = chipchip.prepareStatement("insert into " + dataTable + 
                                               "(experiment,id,probeid,blockno,rowno,colno,channelone,channeltwo,mor) values " +
                                               "(?,NULL,?,?,?,?,?,?,?)");
        if (fournorms) {
            insertStmt.setInt(1,noNormExpt.getDBID());
        } else {
            insertStmt.setInt(1,expt.getDBID());
        }
    }
    public void dropTempTables() throws SQLException {
        Statement stmt = chipchip.createStatement();
        stmt.execute("drop table " + dataTable);
        stmt.execute("drop table " + galTable);
    }

    public void parseArgs(String[] args) throws NotFoundException {
        controlexpt = new ArrayList<ExptNameVersion>();
        Set<String> flags = Args.parseFlags(args);       
        System.err.println("FLAGS are " +flags);
        fournorms = flags.contains("fournorms") || flags.contains("stdnorms");
        median = flags.contains("median") || fournorms;
        mean = flags.contains("mean");
        linefit = flags.contains("linefit") || fournorms;
        nofar = flags.contains("nofar");
        noafe = flags.contains("noafe");
        verbose = flags.contains("verbose");
        bgsub = !flags.contains("nobgsub");
        channelswap = flags.contains("channelswap") || flags.contains("dyeswap");
        galFileOverride = Args.parseString(args,"galfile",null);
        crossIPtoWCE = Args.parseDouble(args,"iptowce",Double.NaN);
        crossWCEtoIP = Args.parseDouble(args,"wcetoip",Double.NaN);
        controlexpt = Args.parseENV(args,"controlexpt");
        designname = Args.parseString(args,"designname",null);
        species = Args.parseGenome(args).getFirst().getName();
        cellsone = Args.parseString(args,"cellsone",null);
        cellstwo = Args.parseString(args,"cellstwo",null);
        conditionone = Args.parseString(args,"conditionone",null);
        conditiontwo = Args.parseString(args,"conditiontwo",null);
        factorone = Args.parseString(args,"factorone",null);
        factortwo = Args.parseString(args,"factortwo",null);
        blockoffset = Args.parseInteger(args,"blockoffset",0);
        List<ExptNameVersion> expts = Args.parseENV(args,"expt");
        name = expts.get(0).getName();
        version = expts.get(0).getVersion();
        replicate = expts.get(0).getReplicate();
        fragdistname = Args.parseString(args,"fragdist",null);
        scratchTablespace = Args.parseString(args,"scratchTablespace","scratch");
        try {
            medianNormRegions = Args.readLocations(args,"medianregions");
            medianLFNormRegions = Args.readLocations(args,"medianlfregions");
        } catch (IOException e) {
            System.err.println("Couldn't open file specified as --medianregions : " + e.toString());
        }
        if (fragdistname != null) {
            String[] pieces = fragdistname.split(";");
            if (pieces.length < 2) {
                throw new IllegalArgumentException("Can't parse fragdist name and version from " + fragdistname);
            }
            fragdistname = pieces[0];
            fragdistversion  = pieces[1];            
        } else {
            fragdistname = Args.parseString(args,"fragdistname",null);
            fragdistversion = Args.parseString(args,"fragdistversion",null);
        }
        ipfiles = Args.parseList(args,"ip");
        wcefiles = Args.parseList(args,"wce");
        bothfiles = Args.parseFile(args);

        if (ipfiles.size() != wcefiles.size()) {
            System.err.println("There are " + ipfiles.size() + " IP files but there are " + 
                               wcefiles.size() + " WCE files");
            throw new IllegalArgumentException("Number of IP files doesn't match number of WCE files");
        }
        try {
            designid = loader.loadArrayDesign(designname,null).getDBID();            
        } catch (SQLException e) {
            throw new NotFoundException("Can't get array design " + designname);
        }
    }

    public Experiment getExperiment (String exptversion) throws SQLException, NotFoundException {
        return loader.getExperiment(name,exptversion,replicate,species,fragdistname,fragdistversion,
                                    factorone,factortwo,cellsone,cellstwo,conditionone,conditiontwo,true);
    }

    /* fills in exptToGenome automatically by looking at which genomes the probes for 
       this experiment's array design have been mapped to */
    public void mapToGenomes(int dbid) throws SQLException {
        PreparedStatement ps = chipchip.prepareStatement("select unique(chromosome) from probelocation where id in (select id from probedesign where arraydesign = ? and rownum < 1000)");
        ps.setInt(1,designid);
        ResultSet rs = ps.executeQuery();
        ArrayList<Integer> chromids = new ArrayList<Integer>();
        while (rs.next()) {
            chromids.add(rs.getInt(1));
        }
        System.err.println("Got chromids for design " + designid + " as " + chromids);
        rs.close();
        ps.close();
        Statement stmt = core.createStatement();
        StringBuffer buf = new StringBuffer();
        buf.append(chromids.get(0));
        for (int i = 1; i < chromids.size(); i++) {
            buf.append("," + chromids.get(i));
        }
        System.err.println("SQL is select unique(genome) from chromosome where id in (" + buf.toString() + ")");
        rs = stmt.executeQuery("select unique(genome) from chromosome where id in (" + buf.toString() + ")");
        PreparedStatement exists = chipchip.prepareStatement("select count(*) from exptToGenome where experiment = ? and genome = ?");
        PreparedStatement insert = chipchip.prepareStatement("insert into exptToGenome(experiment,genome) values(?,?)");
        exists.setInt(1,dbid);
        insert.setInt(1,dbid);
        while (rs.next()) {
            exists.setInt(2,rs.getInt(1));            
            ResultSet e2grs = exists.executeQuery();
            e2grs.next();
            if (e2grs.getInt(1) == 0) {
                insert.setInt(2,rs.getInt(1));
                insert.execute();
            }
            e2grs.close();
        }
        rs.close();
        exists.close();
        insert.close();
    }

    /* do crosstalk normalization if crossIPtoWCE and crossWCEtoIP have been set.
       This tries to account for the Cy5 material that's picked up in the Cy3 channel and
       vice versa.
     */
    public void crossTalkNormalization() throws SQLException {
        if (!Double.isNaN(crossIPtoWCE) && !Double.isNaN(crossWCEtoIP)) {
            System.err.println("Doing crosstalk norm");
            Statement stmt = chipchip.createStatement();
            double z = (1 - crossIPtoWCE * crossWCEtoIP);
            stmt.execute(String.format("update " + dataTable + " set channelone = (channelone - %f * channeltwo)/%f, " +
                                       " channeltwo = (channeltwo - %f * channelone)/%f",
                                       crossWCEtoIP, z,
                                       crossIPtoWCE, z));
            stmt.execute("update " + dataTable + " set ratio = channelone / channeltwo");
            stmt.close();
        }        
    }
    /*  perform a multaplicative normalization such that the median or mean intensity in
        both channels is the same.  This helps to balance the channels in the case where
        one contained more dye than the other
    */
    public void medianOrMeanNorm() throws SQLException {
        if (median || mean) {
            System.err.println("doing median norm");
            Statement stmt = chipchip.createStatement();
            double ip = 1, wce = 1;
            ResultSet rs = null;
            System.err.println("Doing medianOrMeanNorm with median=" + median + " and mean =" + mean);
            String where = "";
            if (medianNormRegions != null && medianNormRegions.size() > 0) {
                where = " where " + getProbesClause(medianNormRegions);
            }
            if (median) {
                rs = stmt.executeQuery("select median(channelone), median(channeltwo) from " + dataTable + where);
            } else if (mean) {
                rs = stmt.executeQuery("select mean(channelone), mean(channeltwo) from " + dataTable + where);
            }
            rs.next();
            ip = rs.getDouble(1);
            wce = rs.getDouble(2);
            rs.close();
            double diff = wce / ip;
            if (Double.isNaN(ip) || Double.isNaN(wce) || Double.isNaN(diff)) {
                throw new RuntimeException("Got a NaN : ip=" + ip + " wce=" + wce + " diff=" + diff);
            }
            System.err.println("Doing median normalization with diff=" + diff);
            stmt.execute("update " + dataTable + " set channelone = channelone * " + diff);            
            stmt.execute("update " + dataTable + " set ratio = channelone / channeltwo");
            stmt.close();
        }        
    }

    /* The linefitting normalization assumes that the best fit line to the 2D data (cy3 on one axis and
       cy5 on the other) should be cy5 = cy3.  This method fits a line to the data and then
       rotates/recenters the data based on that fitted line such that the best fit would
       now be cy5 = cy3.
       Note that this method can turn low intensities into negative intensities, which will make a mess
       later on if you deal with log-intensity.

       TODO this is an oracle specific method
    */
    public void LinefitNorm() throws SQLException, NotFoundException {
        if (linefit) {
            Statement stmt = chipchip.createStatement();
            double factor, intercept;
            stmt.execute("update " + dataTable + " set channelone = log(2,channelone), channeltwo = log(2,channeltwo)");
            String sql = "select to_number(atan(REGR_SLOPE(channelone, channeltwo)) - atan(1)), to_number(REGR_INTERCEPT(channelone,channeltwo)) from "
                + dataTable + " where channelone - channeltwo < 3 and channelone - channeltwo > -3";
            if (medianLFNormRegions != null && medianLFNormRegions.size() > 0) {
                sql += " and " + getProbesClause(medianLFNormRegions);
            }
            System.err.println("Going to normalize with " + sql);
            ResultSet rs = stmt.executeQuery(sql);
            rs.next();
            factor = rs.getDouble(1);
            intercept = rs.getDouble(2);
            System.err.println(String.format("Doing linefit with intercept=%.2f and slope=%.2f", intercept, factor));

            stmt.execute("update " + dataTable + " set channelone = channelone - " + intercept);                
            stmt.execute("update " + dataTable + " set channelone = power(2,sqrt(channelone*channelone + channeltwo*channeltwo) * sin(atan(channelone/channeltwo) - " + factor + ")), " +
                         "channeltwo = power(2,sqrt(channelone*channelone + channeltwo*channeltwo) * cos(atan(channelone/channeltwo) - " + factor + " ))");
            stmt.execute("update " + dataTable + " set ratio = channelone / channeltwo");
            stmt.close();
        }        
    }

    /* control experiment normalization assumes that one or more control experiments represent a signal
       that is present (at some level) in this experiment.  This method performs linear regression, trying to predict
       this experiment based on the control experiment and then subtracts out the portion of the ratio
       that it believes came from the control experiment.  
       This method is useful, eg, for chipchip to account for a lack of specificity in the antibody.  The control
       experiment might be a mock IP or proteinG or similar to represent the parts of the genome that will be pulled
       down by any ip

       TODO this is oracle specific
    */
    private void controlExptNormalize() throws SQLException, NotFoundException {
        if (controlexpt.size() == 0) {
            return;
        }
        System.err.println("Doing control expt normalization");
        long time = System.currentTimeMillis() % 100000;
        String user = chipchip.getMetaData().getUserName();
        String normTable = user + ".normtmp" + time;
        String createNormTableSQL = "create table " + normTable + " (probe number(10) primary key, expt binary_float, control binary_float)  organization index tablespace "
                                     + scratchTablespace + "  nologging";
        Statement stmt = chipchip.createStatement();
        stmt.execute(createNormTableSQL);

        String ixname = dataTable + "_ix";
        stmt.execute("create index " + ixname + " on " + dataTable + "(id)");

        String exptidstring;
        List<Integer> controlexptids = new ArrayList<Integer>();
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        for (ExptNameVersion env : controlexpt) {
            try {
                Collection<Experiment> expts = loader.loadExperiment(env);
                for (Experiment expt : expts) {
                    controlexptids.add(expt.getDBID());
                }
            } catch (NotFoundException e) {
                stmt.execute("drop table " + normTable);
                stmt.close();
                throw e;
            }
        }
        if (controlexptids.size() == 1) {
            exptidstring = " data.experiment = " + controlexptids.get(0);
        } else {
            StringBuilder sb = new StringBuilder("experiment in (");
            for (int i = 0; i < controlexptids.size() - 1; i++) {
                sb.append(controlexptids.get(i));
                sb.append(",");
            }
            sb.append(controlexptids.get(controlexptids.size() -1));
            sb.append(")");
            exptidstring = sb.toString();
        }
        System.err.println("Data table is " + dataTable + " and norm table is " + normTable);
        String getControlDataSQL = "insert into " + normTable + "(probe, control) (select probe, ln(avg(ratio)) from data where " +
                                   exptidstring + " group by probe)";
        System.err.println("Executing " + getControlDataSQL);
        stmt.execute(getControlDataSQL);
        String updateTmpSQL = " update " + normTable + " set expt = (select ln(" + dataTable + ".ratio) from " + dataTable + "  where " + dataTable +
                              ".id = " + normTable+".probe)";
        System.err.println("Executing " + updateTmpSQL);
        stmt.execute(updateTmpSQL);
        String regressSQL = "select to_number(REGR_SLOPE(expt, control)), to_number(REGR_INTERCEPT(expt,control)), REGR_R2(expt,control) from " + normTable + " where expt is not NaN and control is not NaN";
        System.err.println("Executing " + regressSQL);
        ResultSet rs = stmt.executeQuery(regressSQL);
        rs.next();
        Double slope = rs.getDouble(1);
        Double intercept = rs.getDouble(2);
        Double r2 = rs.getDouble(3);
        System.err.println("r2 from control to expt regression is " + r2);
        String subtractPredictedSQL = " update " + normTable + " set expt = exp(expt - (" + intercept + " + control * " + slope + "))";
        System.err.println("Executing " + subtractPredictedSQL);
        stmt.execute(subtractPredictedSQL);
        String saveBackSQL = " update " + dataTable + " set ratio = (select expt from " + normTable + " where " + dataTable + ".id = " + normTable + ".probe)";
        System.err.println("Executing " + saveBackSQL);
        stmt.execute(saveBackSQL);
        stmt.execute("drop table " + normTable);
        rs.close();
        stmt.close();
    }

    /* Copies all rows from the temporary data table to the permanent data table using the specified experiment ID
    */
    public void moveFromTemp(int dbid) throws SQLException {
        Statement stmt = chipchip.createStatement();
        stmt.execute("insert into data(experiment,probe,channelone,channeltwo,mor,channelratio,ratio) " +
                     "(select " + dbid + ", id, channelone, channeltwo, mor, channelratio, ratio from " + dataTable + ")");

    }
    /* assigns probe DBIDs to the temporary data table based on the temporary gal table.
       If fixGeometry is true, then try to figure ouf the correct geometry transformation.  This
       has the side effect of checking that the block/row/col in the galfile matches the block/row/col
       in the data file for a given probeid.

       If fixGeometry is false, then probeids are assigned *soley* based on the block/row/col.

       TODO this may be oracle specific
    */
    public void assignProbeIDs(int galfile, boolean fixGeometry) throws SQLException, UnknownGeometryException {
        Statement stmt = chipchip.createStatement();
        try {
            stmt.execute("drop index " + galIndex);
        } catch (SQLException e) {
            // ignore it, since the index may not exist.
        }
        stmt.execute("delete from " + galTable);
        stmt.execute("insert into " + galTable + "(id,probeid,blockno,rowno,colno) " +
                     "(select id, probeid, blockno, rowno, colno from probedesign where " + 
                     " arraydesign = " + designid + " and galfile = " + galfile + ")");
        stmt.execute("create index " + galIndex + " on " + galTable + "(blockno,rowno,colno)");
        if (fixGeometry) {
            fixGeometry();
        }
        stmt.execute("update "+dataTable +" set id = (select id from "+galTable+" where " +
                     ""+galTable+".blockno = "+dataTable+".blockno and "+galTable+".colno = "+dataTable+".colno " +
                     "and "+galTable+".rowno = "+dataTable+".rowno) where id is null");
        ResultSet rs = stmt.executeQuery("select count(*) from " + dataTable + " where id is null");
        rs.next();
        int noID = rs.getInt(1);
        System.err.println("Going to delete " + noID + " probe values because they don't have an ID.  galfileid was " + galfile);
        rs.close();
        stmt.execute("delete from " + dataTable + " where id is null");
        stmt.close();
    }
    /* examines a subset of the entries in the temporary data table to see which geometry transformation
       makes the most data entries match galfile entries.  This transformation is then applied to 
       all of the temporary data table.

       Examples of geometry transformations are
       - flipped and rotated
       - interleaved blocks
       - side by side blocks
    */
    public void fixGeometry () throws SQLException, UnknownGeometryException{
        boolean isplain = false, isafe = false, isfar = false, ishalf = false, isside = false;
        int maxrow, maxcol;
        Statement stmt = chipchip.createStatement();       
        HashSet<ProbeCacheEntry> probecache = new HashSet<ProbeCacheEntry>();
        ResultSet rs = stmt.executeQuery("select blockno, rowno, colno, probeid from " + galTable);
        while (rs.next()) {
            ProbeCacheEntry entry = new ProbeCacheEntry(rs.getInt(1),rs.getInt(2),rs.getInt(3),rs.getString(4));
            probecache.add(entry);
        }
        rs.close();
        rs = stmt.executeQuery("select max(rowno), max(colno) from " + galTable);
        rs.next();
        maxrow = rs.getInt(1);
        maxcol = rs.getInt(2);
        rs.close();
        if (verbose) {
            System.err.println("maxrow= " + maxrow + ",  maxcol=" + maxcol);
        }


        PreparedStatement actualPos = chipchip.prepareStatement("select blockno, colno, rowno from probedesign where arraydesign = " +
                                                                designid + " and probeid = ?");
        
        int plain = 0, far = 0, afe = 0, afeandfar = 0, half = 0, sidebyside = 0, sideandfar = 0;
        final int max = 3000;
        ArrayList<Integer> blockOffsetCandidates = new ArrayList<Integer>();
        rs = stmt.executeQuery("select probeid, blockno, rowno, colno from (select * from " + dataTable + " where id is null) where rownum < " + max);
        while (rs.next()) {
            String thisid = rs.getString(1);
            int thisblock = rs.getInt(2);
            int thisrow = rs.getInt(3);
            int thiscol = rs.getInt(4);
            if (thisid.matches(".*DarkCorner.*") ||
                thisid.matches(".*3xSLv1.*")) {continue;}
            boolean any = false;
            ProbeCacheEntry trial = new ProbeCacheEntry();
            trial.probeid = thisid;

            trial.block = thisblock; trial.row = thisrow; trial.col = thiscol;
            if (probecache.contains(trial)) {
                any = true;
                plain++;
            }

            trial.row = maxrow - thiscol + 1;
            trial.col = maxcol - thisrow + 1;
            if (probecache.contains(trial)) {
                any = true;
                far++;
            }

            trial.row = (int)((thisrow + 1) / 2);
            trial.col = thiscol * 2 - thisrow % 2;
            if (probecache.contains(trial)) {
                any = true;
                afe++;
            }

            trial.row = maxrow - (thiscol * 2 - (thisrow % 2)) + 1;
            trial.col = maxcol - ((int)((thisrow + 1) / 2)) + 1;
            //            System.err.println("  AFE and FAR gives col=" + trial.col + " and row= " + trial.row + " from col= " + thiscol +", row=" + thisrow + "    maxcol=" + maxcol + ",  maxrow=" + maxrow);
            if (probecache.contains(trial)) {
                any = true;
                afeandfar++;
            }

            trial.row = thisrow;
            trial.col = thiscol * 2 - thisrow % 2;
            if (probecache.contains(trial)) {
                any = true;
                half++;
            }

            if (blockoffset != 0) {
                trial.row = thisrow + (thisblock - 1) * blockoffset;
                trial.col = thiscol;
                trial.block = 1;
                if (probecache.contains(trial)) {
                    any = true;
                    sidebyside++;                    
                }
            } else {
                int minoff, maxoff;
                trial.block = 1;
                if (thisblock == 1) {
                    minoff = 0; maxoff = 1;
                } else {
                    minoff = 0; maxoff = 1000;
                }
                for (int off = minoff; off < maxoff; off++) {
                    trial.row = thisrow + (thisblock -1) * off;
                    trial.col = thiscol;
                    if (probecache.contains(trial)) {
                        any = true;
                        sidebyside++;
                        blockOffsetCandidates.add(off);
                    }
                }
            }

            trial.row = thisrow + (thisblock - 1) * (blockoffset > 0 ? blockoffset : (maxcol / 2));
            trial.col = thiscol;
            trial.row = maxrow - trial.col + 1;
            trial.col = maxcol - trial.row + 1;
            trial.block = 1;
            if (probecache.contains(trial)) {
                any = true;
                sideandfar++;
            }

            if (!any && verbose) {
                actualPos.setString(1,thisid);
                ResultSet aprs = actualPos.executeQuery();
                System.err.print("No transformation can put " + thisid +" at " + thisblock + "," + thiscol + "," + thisrow);
                if (aprs.next()) {
                    System.err.println(" instead of " + aprs.getInt(1) + "," + aprs.getInt(2) + "," + aprs.getInt(3))                    ;
                } else {
                    System.err.println(" because probeid doesn't exist in design");
                }
                aprs.close();
            }
        }

        boolean anygeometry = false;
        if (plain > far &&
            plain > afe && 
            plain > half &&
            plain >= sidebyside &&
            plain > sideandfar &&
            plain > afeandfar) {
            String message = "I think this file is plain: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (plain < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isplain = true;
            anygeometry = true;
        }
        
        if (far > plain &&
            far > afe && 
            far > half &&
            far > sidebyside &&
            far >= sideandfar &&
            far > afeandfar) {
            String message = "I think this file is far: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (far < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isfar = true;
            anygeometry = true;
        }
        
        if (afe > plain &&
            afe > far && 
            afe > half &&
            afe > sidebyside &&
            afe > sideandfar &&
            afe > afeandfar) {
            String message = "I think this file is afe: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (afe < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isafe = true;
            anygeometry = true;
        }
        
        if (afeandfar > plain &&
            afeandfar > far && 
            afeandfar > half &&
            afeandfar > sidebyside &&
            afeandfar > sideandfar &&
            afeandfar > afe) {
            String message = "I think this file is afeandfar: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (afeandfar < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isafe = true;
            isfar = true;
            anygeometry = true;
        }
        if (half > plain &&
            half > far && 
            half > afeandfar &&
            half > sidebyside &&
            half > sideandfar &&
            half > afe) {
            String message = "I think this file is half: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (half < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            ishalf = true;
            anygeometry = true;
        }
        if (sidebyside > plain &&
            sidebyside > far && 
            sidebyside > afeandfar &&
            sidebyside > half &&
            sidebyside > sideandfar &&
            sidebyside > afe) {
            String message = "I think this file is sidebyside: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (sidebyside < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isside = true;
            anygeometry = true;
            if (blockoffset == 0) {
                System.err.println("Guessing block offset value as mode of values that worked");
                int bestval = -1, most = 0;
                for (int i : blockOffsetCandidates) {
                    int count = Collections.frequency(blockOffsetCandidates,i);
                        if (count > most) {
                            bestval = i;
                            most = count;
                        }
                }
                System.err.println("Guesssing " + bestval);
                if (most < (int)(.8 * max)) {
                    throw new UnknownGeometryException("Guessing didn't work because " + bestval+ " didn't work on enough probes");
                }
                blockoffset = bestval;
            }            
        }
        if (sideandfar > plain &&
            sideandfar > far && 
            sideandfar > afeandfar &&
            sideandfar > half &&
            sideandfar > sideandfar &&
            sideandfar > afe) {
            String message = "I think this file is sideandfar: " + plain + "," + far + "," + afe + "," + 
                    afeandfar + "," + half + "," + sidebyside + "," + sideandfar;
            System.err.println(message);
            if (sideandfar < (int)(.8 * max)) {
                throw new UnknownGeometryException("Too few matched: " + message);
            }
            isside = true;
            isfar = true;
            anygeometry = true;
        }                       

        if (!anygeometry) {
            throw new UnknownGeometryException("Couldn't figure out any geometry from plain=" + plain + ", far=" + far + ", afeandfar=" +
                                               afeandfar + ", sideandfar=" + sideandfar + ", half=" + half + ", afe=" + afe);
        }

        actualPos.close();
        int offset;
        if (isside && isfar) {
            offset = blockoffset > 0 ? blockoffset : (maxcol / 2);
        } else if (isside) {
            offset = blockoffset > 0 ? blockoffset : maxrow;
        }
        if (isafe) {
            stmt.execute("update " + dataTable + " set rowno = trunc((rowno + 1)/2), colno = (colno * 2 - mod(rowno,2))");
        }
        if (isside) {
            stmt.execute("update "+dataTable+" set blockno = 1, rowno = (rowno + (blockno - 1) * "+blockoffset+") ");
        }
        if (isfar) {
            stmt.execute("update "+dataTable+" set rowno = "+maxrow+" - colno + 1, colno = "+maxcol+" - rowno + 1 ");
        }
        if (ishalf) {
            stmt.execute("update "+dataTable+" set colno = trunc(colno * 2) - mod(rowno,2) ");
        }
        stmt.close();
    }

    /** 
     * returns an SQL clause of the form
     *    id in (...) 
     * that selects on the probes in the input regions.  If the regions are StrandedRegions, then
     * only probes on the appropriate strand are used.
     */
    private String getProbesClause (List<Region> regions) throws SQLException {
        Set<Integer> pids = new HashSet<Integer>();
        PreparedStatement getProbes = chipchip.prepareStatement("select pd.id from probedesign pd, probelocation pl where pl.id = pd.id and pl.chromosome = ? and pl.startpos >= ? and pl.startpos <= ?");
        PreparedStatement getProbesStranded = chipchip.prepareStatement("select pd.id from probedesign pd, probelocation pl where pl.id = pd.id and " +
                                                                        " pl.chromosome = ? and pl.startpos >= ? and pl.startpos <= ? and pl.strand = ?");
        ResultSet rs;
        for (Region r : regions) {
            if (r instanceof StrandedRegion) {
                StrandedRegion sr = (StrandedRegion) r;
                getProbesStranded.setInt(1, r.getGenome().getChromID(r.getChrom()));
                getProbesStranded.setInt(2, r.getStart());
                getProbesStranded.setInt(3, r.getEnd());
                getProbesStranded.setString(4, (new Character(sr.getStrand())).toString());
                rs = getProbesStranded.executeQuery();
                System.err.println("Got Stranded Region " + sr);
            } else {
                getProbes.setInt(1, r.getGenome().getChromID(r.getChrom()));
                getProbes.setInt(2, r.getStart());
                getProbes.setInt(3, r.getEnd());
                rs = getProbes.executeQuery();
            }
            while (rs.next()) {
                pids.add(rs.getInt(1));
            }
            rs.close();
        }
        getProbes.close();
        getProbesStranded.close();
        StringBuffer clause = new StringBuffer();
        clause.append(" id in (");
        boolean first = true;
        for (int p : pids) {
            if (first) {
                clause.append(p);
            } else {
                clause.append(", " + p);
            }
            first = false;
        }
        clause.append(")");
        clause.toString();
        System.err.println("Using the following probes for median/mean normalization " + clause.toString());
        return clause.toString();
    }

    // stores the probe observations into the temporary data table and returns the galfileid
    public int addFile(String fname) throws SQLException, IOException, UnknownFileTypeException {
        int galid;
        if (fname.matches(".*gpr$")) {            
            AddGPRHandler handler = new AddGPRHandler(chipchip,
                                                      insertStmt);
            GPRFile file = new GPRFile(fname,handler);            
            file.parse();
            String galfname;
            if (galFileOverride != null) {
                galfname = galFileOverride;
            } else {
                galfname = handler.getHeader("GalFile");
            }

            System.err.println("Galfile name is " + galfname);
            galid = loader.loadGALFile(galfname).getDBID();
        } else if (fname.matches(".*afe") ||
                   fname.matches(".*txt")) {
            AddAFEHandler handler = new AddAFEHandler(chipchip,
                                                      insertStmt);
            handler.setBackgroundSubtraction(bgsub);
            AFEFile file = new AFEFile(fname,
                                       handler);
            file.parse();
            String galfname = handler.getHeader("Grid_Name");
            System.err.println("Galfile name is " + galfname);
            galid = loader.loadGALFile(galfname).getDBID();
        } else {
            throw new UnknownFileTypeException("Can't figure out file type of " + fname);
        }
        
        
        return galid;
    }
    // stores the probe observations into the temporary data table and returns the galfileid
    public int addFiles(String ipfname, String wcefname) throws SQLException, IOException, UnknownFileTypeException {
        int galid;
        if (ipfname.matches(".*cel$") || ipfname.matches(".*cel.txt$")) {
            throw new UnknownFileTypeException("Can't handle cel files yet");
        } else if (ipfname.matches(".*pair") ||ipfname.matches(".*pair.txt")) {
            AddPairHandler handler = new AddPairHandler(chipchip,
                                                        insertStmt);
            PairFile file = new PairFile(ipfname,
                                         wcefname,
                                         handler);
            file.parse();
            String galfname = handler.getHeaderOne("designfile");
            galid = loader.loadGALFile(galfname).getDBID();
        } else {
            throw new UnknownFileTypeException("Can't figure out file type of " + ipfname);
        }

        return galid;
    }
}

/* a ProbeCacheEntry is combination of array coordinates (block, row column) and a probe id.
   When a design is loaded, we load all four values for each spot in the database.  When we load
   a datafile and look up a spot based on the coordinates, the probeids should match.  Similarly, if
   we look up based on probeid, the coordinates should match.  The ProbeCacheEntry class is data used
   in this matching and checking process */   
class ProbeCacheEntry {
    public int block, row, col;
    public String probeid;
    public ProbeCacheEntry() {}
    public ProbeCacheEntry (int b, int r, int c, String i) {
        block = b;
        row = r;
        col = c;
        probeid = i;
    }
    public boolean equals(Object o) {
        if (o instanceof ProbeCacheEntry) {
            ProbeCacheEntry other = (ProbeCacheEntry)o;
            return (other.block == block &&
                    other.row == row &&
                    other.col == col &&
                    other.probeid.equals(probeid));
        } else {
            return false;
        }
    }
    public int hashCode() {
        return block * 17 + row * 5 + col + probeid.hashCode();
    }
}
class UnknownGeometryException extends Exception {
    public UnknownGeometryException(String e) {super(e);}
}
