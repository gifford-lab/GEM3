/*
 * Created on Apr 5, 2007
 */
package edu.mit.csail.cgs.datasets.expression.parsing;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.TimePoint;
import edu.mit.csail.cgs.utils.ArgParser;

public class TextTableExpressionParser {
    
    public static void main(String[] args) { 
        ArgParser ap = new ArgParser(args);
        try {
            MetadataLoader metaLoader = new MetadataLoader();
            ExpressionInserter inserter = new ExpressionInserter();
            ExpressionLoader loader = inserter.getLoader();
                        
            String exptName = ap.getKeyValue("expt");
            int index = ap.hasKey("index") ? Integer.parseInt(ap.getKeyValue("index")) : 0;
            boolean logScale = ap.hasKey("scale") ? ap.getKeyValue("scale").equals("log") : true;
            int type = ap.hasKey("type") ? Integer.parseInt(ap.getKeyValue("type")) : -1;
            Cells cells = metaLoader.getCells(ap.getKeyValue("cells"));
            Condition condition = metaLoader.getCondition(ap.getKeyValue("condition"));

            File inputFile = new File(ap.getKeyValue("input"));
            ProbePlatform platform = loader.loadPlatform(ap.getKeyValue("platform"));

            TextTableExpressionParser parser = new TextTableExpressionParser(inputFile, platform);
            parser.setLogScale(logScale);
            parser.setType(type);
            parser.setCells(cells);
            parser.setCondition(condition);
            
            parser.insertIntoDB(inserter, exptName, index);
            
            loader.close();
            inserter.close();
            metaLoader.close();

        } catch (SQLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private Map<String,Vector<Double>> values;
    
    private Map<String,String> params;
    private Cells cells;
    private Condition condition;
    private TimePoint timepoint;
    private ProbePlatform platform;
    private boolean logScale;
    private int type;

    public TextTableExpressionParser(File tableFile, ProbePlatform plat) throws IOException {
        values = new HashMap<String,Vector<Double>>();
        params = new HashMap<String,String>();
        cells = null;
        condition = null;
        timepoint = null;
        platform = plat;
        type = -1;
        logScale = false;
        
        BufferedReader br = new BufferedReader(new FileReader(tableFile));
        String line = null;
        while((line = br.readLine()) != null) { 
            String[] array = line.split("\\s+");
            String id = array[0];
            Vector<Double> vals = new Vector<Double>();
            for(int i = 1; i < array.length; i++) { vals.add(Double.parseDouble(array[i])); }
            values.put(id, vals);
        }
        br.close();
    }

    public void addParameter(String k, String v) { params.put(k, v); }
    public void setLogScale(boolean ls) { logScale = ls; }
    public void setType(int t) { type = t; }
    public void setCells(Cells c) { cells = c; }
    public void setCondition(Condition c) { condition = c; }
    
    public void insertIntoDB(ExpressionInserter inserter, String name, int index) throws SQLException {
        ExpressionLoader loader = inserter.getLoader();
        inserter.beginTransaction();

        int exptID = inserter.insertExperiment(name, type, logScale, cells, condition, timepoint, platform);
        System.out.println("Inserted Experiment \"" + name + "\"");
        
        Experiment expt = loader.loadExperiment(exptID);
        
        System.out.println("Building Probe->Value Map...");
        Map<String,Probe> probeMap = loader.loadProbes(platform, values.keySet());
        Map<Probe,Double> probeValues = new HashMap<Probe,Double>();
        for(String k : values.keySet()) {
            if(probeMap.containsKey(k) && values.get(k).size() > index) { 
                Probe p = probeMap.get(k);
                probeValues.put(p, values.get(k).get(index));
            }
        }
        
        System.out.println("Inserting Measurements...");
        inserter.insertMeasurements(expt, probeValues);
        
        inserter.commitTransaction();
        loader.close();
        System.out.println("Done.");
    }
}
