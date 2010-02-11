/*
 * Created on May 24, 2005
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.io.parsing.FunctionalCategories;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class MIPSTextFunctionalCategories implements FunctionalCategories {

	public static void main(String[] args) { 
		ArgParser ap = new ArgParser(args);
		if(!ap.hasKey("scheme") || !ap.hasKey("data")) { 
			System.err.println("USAGE: --scheme and --data options are needed.");
			System.exit(1);
		}

		File schemeFile = new File(ap.getKeyValue("scheme"));
		File dataFile = new File(ap.getKeyValue("data"));
		MIPSTextFunctionalCategories fc = 
			new MIPSTextFunctionalCategories(schemeFile, dataFile);
		interactive(fc);
	}

	public static void interactive(FunctionalCategories fc) { 
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			System.out.print(">"); System.out.flush();
			String line = null;
			while((line = br.readLine()) != null) { 
				line = line.trim();	
				if(!fc.hasObject(line)) { 
					System.out.println("No Such Object: " + line);
				} else { 
					Collection<String> cats = fc.classify(line);
					for(String cat : cats) { 
						System.out.println(cat + " --> " + fc.getDescription(cat));
					}
				}
				System.out.print(">"); System.out.flush();
			}
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
			throw new RuntimeException(ie);
		}
	}
    
    private File scheme, data;
    private Map<String,MIPSCategORF> orf2Categ;
    private Map<String,MIPSCategory> name2Categ;
    
    public MIPSTextFunctionalCategories(File scheme, File data) { 
        parseScheme(scheme);
        parseData(data);
    }
    
    private void parseScheme(File input) { 
        try {
            int count = 0;
            name2Categ = new HashMap<String,MIPSCategory>();
            BufferedReader br = new BufferedReader(new FileReader(input));
            String line = null;
            while((line = br.readLine()) != null) {
                line = line.trim();
                if(!line.startsWith("#")) { 
                    MIPSCategory categ = new MIPSCategory(line);
                    name2Categ.put(categ.getName(), categ);
                    count += 1;
                }
            }
            br.close();
            System.err.println("# Categories: " + count);
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
            throw new RuntimeException(ie);
        }
    }
    
    private void parseData(File input) { 
        try {
            int count = 0;
            int ccount = 0;
            orf2Categ = new HashMap<String,MIPSCategORF>();
            BufferedReader br = new BufferedReader(new FileReader(input));
            String line = null;
            while((line = br.readLine()) != null) {
                line = line.trim();
				//System.out.println("Line {" + line + "}");
                if(!line.startsWith("#")) {
                    String[] entries = line.split("\\|");
					String key = entries[0].toUpperCase();
					//System.out.print("[" + key + "]");
                    if(orf2Categ.containsKey(key)) { 
						//System.out.println(" --> EXISTS");
                        MIPSCategORF orf = orf2Categ.get(key);
                        orf.addLine(line);
                    } else { 
						//System.out.println(" --> NEW");
                        MIPSCategORF orf = new MIPSCategORF(line);
                        orf2Categ.put(orf.getORF(), orf);
						//System.out.println("\tADDED: " + orf.getORF());
                        count += 1;
                    }
                    ccount += 1;
                }
            }            
            br.close();
            System.err.println("# Categorized ORFs: " + count);
            System.err.println("# ORF-Categories: " + ccount);
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
            throw new RuntimeException(ie);
        }        
    }

	public boolean hasObject(String object) { 
		return orf2Categ.containsKey(object);
	}

	public boolean hasCategory(String name) { 
		return name2Categ.containsKey(name);
	}
    
    public Collection<String> getCategories() { 
        return name2Categ.keySet();
    }

	public Collection<String> getCategoryObjects(String catName) { 
		LinkedList<String> lst = new LinkedList<String>();
		for(String orfName : orf2Categ.keySet()) { 
			MIPSCategORF orf = orf2Categ.get(orfName);
			if(orf.hasCategory(catName)) { 
				lst.addLast(orfName);
			}
		}
		return lst;
	}
    
    public String getDescription(String name) { 
        MIPSCategory categ = name2Categ.get(name);
        return categ.getDesc();
    }
    
    public Collection<String> getParents(String name) { 
        MIPSCategory categ = name2Categ.get(name);
        Set<String> parents = new HashSet<String>();
        for(String key : name2Categ.keySet()) { 
            MIPSCategory c = name2Categ.get(key);
            if(c.isParent(categ)) { parents.add(c.getName()); }
        }
        return parents;
    }
    
    public Collection<String> classify(String object) { 
        return orf2Categ.get(object).getCategories();
    }    
}
