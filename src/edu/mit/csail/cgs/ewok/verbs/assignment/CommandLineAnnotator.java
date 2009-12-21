package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;
import java.io.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingScanGenerator;
import edu.mit.csail.cgs.ewok.nouns.*;

public class CommandLineAnnotator {

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);

		if(!ap.hasKey("genome")) { 
			System.err.println("USAGE: CommandLineAnnotator " + 
					"--genome <genome name> " + 
					"[--bindVersion <binding-scan version> --bindType <binding-scan type>] " + 
					"[--upstream <# bp>] " + 
                    "[--genes genetype] " +
					"[--downstream <# bp>]\n\n" + 
					"\tThe only required argument is the name of the genome (i.e., hg17, or mm6, or some other " + 
					"genome version).  If just the genome name is given, " +
                    "then the tool will print out a list of " + 
					"all the binding scan version/type pairs that are available for that genome, and then quit.\n " + 
                    "\tThe normal operation of the tool occurs when a version/type pair are specified with the " + 
					"--bindVersion and --bindType options.  When this happens, the tool will find every binding " + 
					"event from the scan, as well as all (refGene) genes that are assigned to that event.  " + 
					"Each event is output as a line, with its location, size, and confidence.  Every gene " + 
					"assigned to that event is then output below the event; one line per gene, each line indented " + 
					"by one tab.  Each gene line consists of the location of the gene, its refGene identifier, and " + 
					"the corresponding Entrez Gene ID.\n" + 
					"\tThe assignment of the binding event to the gene depends on the relative distance between the " + 
					"two.  By default, we assign an event to a gene if the event is within 8k bp upstream of the " + 
					"gene or 2k bp downstream of the gene.  These parameters can be changed, however, using the " + 
					"--upstream and --downstream command-line options.\"");
			return;
		}

		String gname = ap.getKeyValue("genome");
		int up = ap.hasKey("upstream") ? Integer.parseInt(ap.getKeyValue("upstream")) : 8000;
		int down = ap.hasKey("downstream") ? Integer.parseInt(ap.getKeyValue("downstream")) : 2000;

		String version = ap.getKeyValue("bindVersion");
		String type = ap.getKeyValue("bindType");
		
		try {
			Genome g = Organism.findGenome(gname);
			java.sql.Connection cxn = g.getUcscConnection();
			TableXRef xref = new TableXRef(cxn, "refLink", "mrnaAcc", "locusLinkId");
			
	        BindingScanLoader loader = new BindingScanLoader();
            RefGeneGenerator refgenes;
            if (ap.hasKey("genes")) {
                refgenes = new RefGeneGenerator(g,ap.getKeyValue("genes"));	        
            } else {
                refgenes = new RefGeneGenerator(g);	        
            }


            refgenes.setUpstreamDownstream(up,down);
	        Collection<BindingScan> scans = loader.loadScans(g);
	        BindingScan selectedScan = null;
	        
	        for(BindingScan bs : scans) { 
	        	if(bs.getVersion().equals(version) && bs.getType().equals(type)) { 
	        		selectedScan = bs;
	        	}
	        }
            
	        if(selectedScan != null) { 
	        	System.err.println("Selected scan: " + selectedScan.toString());
	        
	        	Expander<Region,BindingEvent> caller = new BindingScanGenerator(loader, selectedScan);
	        	BindingEventAnnotations annots = new BindingEventAnnotations(g, caller);
                annots.addAnnotations("genes",
                                      refgenes);

	        	for(int i = 0; i < annots.getNumItems(); i++) { 
	        		BindingEvent evt = annots.getItem(i);

	        		//System.out.println(evt.toString());
	        		System.out.println(evt.getLocationString() + "\t" + 
	        				evt.getSize() + "\t" + 
	        				evt.getConf());
	        		
	        		for(Region r : annots.getAnnotations(evt, "genes")) { 
	        			Gene gene = (Gene)r;
	        			Iterator<String> entrezItr = xref.execute(gene.getID());
	        			
	        			System.out.print("\t" + gene.getLocationString() + "\t" + 
	        					gene.getID()); 
	        			if(entrezItr.hasNext()) { System.out.print("\t"); }
	        			boolean first = true;
	        			while(entrezItr.hasNext()) { 
	        				String entrez = entrezItr.next();
	        				if(!first) { System.out.print(","); }
	        				System.out.print(entrez);
	        				first = false;
	        			}
	        			System.out.println();
	        		}
	        	}
	        } else { 
	        	System.err.println("No scan matches: \"" + type + "\", \"" + version + "\"");
	        	for(BindingScan bs : scans) { 
	        		System.err.println(bs.getVersion() + "\t" + bs.getType());
	        	}
	        }
	        
	        xref.close();
	        DatabaseFactory.freeConnection(cxn);

		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (UnknownRoleException e) {
			e.printStackTrace();
		}
	}
}
