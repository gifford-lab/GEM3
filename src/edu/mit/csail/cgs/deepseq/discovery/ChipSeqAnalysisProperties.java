package edu.mit.csail.cgs.deepseq.discovery;

import java.io.*;
import edu.mit.csail.cgs.utils.json.*;
import edu.mit.csail.cgs.utils.models.*;

/** Class for ChIP-Seq binding mixture properties.
 * Uses Tim's JSON model stuff to handle serializing and deserializing.
 */
public class ChipSeqAnalysisProperties extends Model {

	// default values, will be overwritten when loading the property file.
	// 
    public Boolean TF_binding = true;
    public Integer window_size_factor = 8;
    public Integer max_hit_per_bp = -1;
    
	// parameters for BindingMixture
    public Double q_value_threshold = 2.0;
    public String model="CTCF.empiricalDistribution.txt";
    public Integer mappable_genome_length = (int) 2.08E9; // mouse genome 

//    public Integer update_model_round = 0;		// extra rounds to update empirical dist.
	public Double max_shape = 2.0;
	public Double min_strength = 30.0;
    public Boolean use_single_event = true;
    
    // text file holding the mean and stddev of normal distribution of logKL value for some number of reads
    public String peak_shape_distribution="PeakShapeDistribution.txt";
    
    public Double sparseness=6.0;
    public Boolean use_dynamic_sparseness = true;
//    public Boolean make_hard_assignment = true;  
    
    public Integer bmverbose=1;		// BindingMixture verbose mode
	public Boolean print_mixing_probability = false;

    
    // controls whether some exception stack traces are printed
    private boolean debugging = true;
    /* true when a window is currently open to configure this WarpProperties.
       Used to prevent multiple windows from working on the same properties at once.
    */
    private Boolean configuring;

    public ChipSeqAnalysisProperties() {
        super();
        configuring = Boolean.FALSE;
    }

    /** Saves these properties to the specified file
     */
    public void saveToFile(File f) {
        JSONObject json = asJSON();
        try {
            PrintWriter pw = new PrintWriter(f);
            pw.println(json.toString());
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /** Load these properties from the specified file
     */
    public void loadFromFile(File f) throws FileNotFoundException {
        loadFromStream(new FileReader(f));
    }
    
    public void loadFromStream(InputStreamReader reader) {
        try {
            BufferedReader breader = new BufferedReader(reader);
            JSONObject json = new JSONObject(breader.readLine());
            breader.close();
            reader.close();
            setFromJSON(json);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args)
    {
    	ChipSeqAnalysisProperties obj = new ChipSeqAnalysisProperties();
    	File f = new File("test_properties.txt");
    	obj.saveToFile(f);
    }


}