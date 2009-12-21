package edu.mit.csail.cgs.tools.hypotheses;

import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.HarbisonRegCodeProbe;
import edu.mit.csail.cgs.ewok.verbs.CastingMapper;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.HarbisonRegCodeProbeGenerator;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingExpander;
import edu.mit.csail.cgs.tools.hypotheses.utils.HarbisonCachedBindingExpander;
import edu.mit.csail.cgs.tools.hypotheses.utils.HarbisonFactorIterator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.SetTools;

public class BindingExplorer {
    
    private static String[][] signalArray = { 
        { "Sc Tpk1:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "Tpk1"},
        { "Sc Tpk2:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "Tpk2" },
        { "Sc Kss1:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "Kss1" },
        { "Sc Fus3:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "Fus3" },
        { "Sc Ste5:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "Ste5" },
        { "Sc Hog1:Sc:YPDFA vs WCE:Sc:YPDFA", "11/8/06, default params", "Hog1" },
        { "Sc GCN4:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "GCN4" },
        { "Sc PolII:Sc:YPD vs WCE:Sc:YPD", "11/8/06, default params", "PolII" }
    };
    
    public static void omain(String[] args) { 
        try { 
            MetadataLoader metaloader = new MetadataLoader();
            Genome genome = Organism.findGenome("sacCer1");
            Condition cond = metaloader.getCondition("YPD");
            
            BindingScanLoader bsLoader = new BindingScanLoader();
            LinkedList<BindingScan> scans = new LinkedList<BindingScan>();
            BindingResources recs = new BindingResources(genome);
            
            for(int i = 0; i < signalArray.length; i++) { 
                String version = signalArray[i][0] + "," + signalArray[i][1];
                String type = "BayesBindingGenerator";
                LinkedList<BindingScan> scs = new LinkedList<BindingScan>(bsLoader.loadScans(genome, version, type));
                if(!scs.isEmpty()) { 
                    BindingScan bs = scs.getFirst();
                    scans.addLast(bs);
                    Factor f = metaloader.getFactor(signalArray[i][2]);
                    recs.addExpander(f, new BindingExpander(bsLoader, bs));
                } else { 
                    System.out.println("Couldn't find a BindingScan for " + version + " ** " + type);
                }
            }

            System.out.println("Creating base region & probe lists...");
            Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
            RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(genome, "sgdGene");
            Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
            
            Vector<Region> regions = new Vector<Region>();
            
            while(genes.hasNext()) {
                Gene g = genes.next();
                int tss = g.getTSS();
                Region r = new Region(genome, g.getChrom(), tss-5, tss+5);
                regions.add(r);
            }
            
            System.out.println("Building Binding Summaries...");
            BindingExplorer explorer = new BindingExplorer(recs, regions);
            
            Vector<Factor> fs = new Vector<Factor>(recs.getFactors());
            Map<String,Factor> factorMap = new HashMap<String,Factor>();
            for(Factor f : fs) { factorMap.put(f.getName(), f); }
            
            FactorSurvey survey = new FactorSurvey(factorMap.get("GCN4"), explorer);
            HypothesisTree tree = survey.createTree(0);
            
            LinkedList<BindingHypothesis> optPath = tree.findOptimalPath();
            int optScore = tree.findOptimalScore();
            tree.printTree();
            
            System.out.println("Optimal Path: " + optPath);
            System.out.println("Optimal Score: " + optScore);
            
            metaloader.close();
            
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
	
	public static void main(String[] args) { 
		try { 
			MetadataLoader metaloader = new MetadataLoader();
			Genome genome = Organism.findGenome("sacCer1");
			Condition cond = metaloader.getCondition("YPD");

			HarbisonRegCodeProbeGenerator gen = new HarbisonRegCodeProbeGenerator(genome);

			System.out.println("Creating base region & probe lists...");
			Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
			Iterator<Region> chromRegions = 
				new MapperIterator<NamedRegion,Region>(new CastingMapper<NamedRegion,Region>(), chroms);
			Iterator<HarbisonRegCodeProbe> probeitr = 
				new ExpanderIterator<Region,HarbisonRegCodeProbe>(gen, chromRegions);

			Vector<Region> regions = new Vector<Region>();
			Vector<HarbisonRegCodeProbe> probes = new Vector<HarbisonRegCodeProbe>();
			
			while(probeitr.hasNext()) {
				HarbisonRegCodeProbe p = probeitr.next();
				regions.add(p);
				probes.add(p);
			}
			
			System.out.println("Retrieving Factor list...");
			HarbisonFactorIterator itr = new HarbisonFactorIterator(genome);
			BindingResources recs = new BindingResources(genome);

			int i = 0;
			while(itr.hasNext()) { 
				Factor f = itr.next();
				Expander<Region,BindingEvent> exp = 
					new HarbisonCachedBindingExpander(probes, f, cond);
				recs.addExpander(f, exp);
				i += 1;
			}
			System.out.println("\tAdded " + i + " factors.");
			
			System.out.println("Building Binding Summaries...");
			BindingExplorer explorer = new BindingExplorer(recs, regions);
			
			Vector<Factor> fs = new Vector<Factor>(recs.getFactors());
			Map<String,Factor> factorMap = new HashMap<String,Factor>();
			for(Factor f : fs) { factorMap.put(f.getName(), f); }
			
			FactorSurvey survey = new FactorSurvey(factorMap.get("STE12"), explorer);
            HypothesisTree tree = survey.createTree(0);
            
            LinkedList<BindingHypothesis> optPath = tree.findOptimalPath();
            int optScore = tree.findOptimalScore();
            tree.printTree();
            
            System.out.println("Optimal Path: " + optPath);
            System.out.println("Optimal Score: " + optScore);
			
			metaloader.close();
			
		} catch(SQLException se) { 
			se.printStackTrace(System.err);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
    
    private static SetTools<Region> regionTools;
    
    static { 
        regionTools = new SetTools<Region>();
    }
	
	private BindingResources resources;
	private Vector<Factor> factors;
	private Vector<Region> regions;
    private HashSet<Region> regionSet;
	private Map<Region,LocalBindingSummary> summaries;
	private GlobalBindingSummary globalSummary;
    private Map<Factor,Set<Region>> boundRegions;
	
	public BindingExplorer(BindingResources recs, Vector<Region> regs) { 
		resources = recs;
		factors = new Vector<Factor>(resources.getFactors());
		regions = new Vector<Region>(regs);
        regionSet = new HashSet<Region>(regs);
		summaries = new HashMap<Region,LocalBindingSummary>();
		globalSummary = resources.getSummary();
        boundRegions = new HashMap<Factor,Set<Region>>();
        
        for(Factor f : factors) { boundRegions.put(f, new HashSet<Region>()); }
		
		int i = 0;
		for(Region r : regions) { 
            LocalBindingSummary summ = globalSummary.createSummary(r);
            for(Factor f : summ.getBound()) { 
                boundRegions.get(f).add(r);
            }
			summaries.put(r, summ);
			i += 1;
		}
        
	}
	
	public BindingResources getResources() { return resources; }
	public int getNumFactors() { return factors.size(); }
	public Factor getFactor(int i) { return factors.get(i); }
	public int getNumRegions() { return regions.size(); }
	public Region getRegion(int i) { return regions.get(i); }
	public LocalBindingSummary getBinding(int i) { return summaries.get(regions.get(i)); }
	public GlobalBindingSummary getGlobalSummary() { return globalSummary; }
	public Collection<Region> getAllRegions() { return summaries.keySet(); }
	
	public void updateScoredHypothesis(ScoredHypothesis shyp, Collection<Region> rs) {

        switch(shyp.getHypothesis().getOptimizationClass()) {
        case 1:
            BindingHypothesis.Conditional c = (BindingHypothesis.Conditional)(shyp.getHypothesis());
            Factor f1 = ((BindingHypothesis.SimpleBinding)c.getAntecedent()).getFactor();
            Factor f2 = ((BindingHypothesis.SimpleBinding)c.getConsequent()).getFactor();
            
            Set<Region> inconsistent = regionTools.intersection(boundRegions.get(f1), 
                    regionTools.subtract(regionSet, boundRegions.get(f2)));
            
            for(Region r : rs) { 
                if(summaries.containsKey(r)) { 
                    boolean supported = !inconsistent.contains(r);
                    shyp.addRegion(r, supported);
                }
            }
            
            break;

        default:
            for(Region r : rs) { 
                if(summaries.containsKey(r)) { 
                    boolean supported = shyp.getHypothesis().isSupportedBy(summaries.get(r));
                    shyp.addRegion(r, supported);
                }
            }
            break;
        }
	}
	
	public ScoredHypothesis scoreHypothesis(BindingHypothesis hyp, Collection<Region> rs) {
        //System.out.println("Checking \"" + hyp + "\" against " + rs.size() + " regions.");
		ScoredHypothesis sh = new ScoredHypothesis(hyp);
		updateScoredHypothesis(sh, rs);
		return sh;
	}
	
	public ScoredHypothesis scoreHypothesis(BindingHypothesis hyp) {
		ScoredHypothesis sh = new ScoredHypothesis(hyp);
		for(Region r : summaries.keySet()) { 
			boolean supported = hyp.isSupportedBy(summaries.get(r));
			sh.addRegion(r, supported);
		}
		return sh;
	}

	public int getFalsifyingCount(BindingHypothesis hyp) { 
		int c = 0;
		for(Region r : summaries.keySet()) { 
			if(!hyp.isSupportedBy(summaries.get(r))) { 
				c += 1;
			}
		}
		return c;
	}
}
