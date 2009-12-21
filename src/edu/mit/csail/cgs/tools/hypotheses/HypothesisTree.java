package edu.mit.csail.cgs.tools.hypotheses;

import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.SetTools;

public class HypothesisTree {
	
	private static SetTools<Region> tools;
    private static int branching;
	
	static { 
		tools = new SetTools<Region>();
        branching = 10;
	}

	private BindingExplorer explorer;
	private Set<Region> testableRegions, inconsistentRegions;
	private ScoredHypothesis hypothesis;
	private Map<BindingHypothesis,HypothesisTree> children;
	private HypothesisTree parent;
	
	public HypothesisTree(BindingExplorer exp) { 
		explorer = exp;
		testableRegions = new HashSet<Region>(explorer.getAllRegions());
		hypothesis = null;
		inconsistentRegions = testableRegions;
		children = null;
		parent = null;
	}
	
	public HypothesisTree(BindingExplorer exp, BindingHypothesis hyp, Collection<Region> testable, HypothesisTree treeParent) { 
		explorer = exp;
		testableRegions = new HashSet<Region>(testable);
		hypothesis = explorer.scoreHypothesis(hyp, testableRegions);
		inconsistentRegions = tools.subtract(testableRegions, hypothesis.getSupportingRegions());
		children = null;
		parent = treeParent;
	}
    
    public void printTree() { 
        printTree(System.out, 0);
    }
    
    public void printTree(PrintStream ps) { 
        printTree(ps, 0);
    }
    
    public void printTree(PrintStream ps, int depth) { 
        for(int i = 0; i < depth; i++) { ps.print("  "); }
        ps.println(hypothesis != null ? hypothesis : "<root>");
        
        if(children != null) { 
            TreeSet<ScoredHypothesis> shyps = new TreeSet<ScoredHypothesis>();
            for(BindingHypothesis hyp : children.keySet()) { 
                shyps.add(children.get(hyp).getScoredHypothesis()); 
            }
            int i = 0;
            for(ScoredHypothesis shyp : shyps) {
                ps.print(i + ": ");
                children.get(shyp.getHypothesis()).printTree(ps, depth+1);
                i++;
            }
        }
    }
    
    public void addChildrenAtDepth(Collection<BindingHypothesis> hyps, int depth) { 
        if(depth > 0) {
            if(children == null) { addChildren(hyps); } 
            for(BindingHypothesis hyp : children.keySet()) {
                children.get(hyp).addChildrenAtDepth(hyps, depth-1);
            }
        } else { 
            addChildren(hyps);
        }
    }
	
	public void addChildren(Collection<BindingHypothesis> hyps) { 
		for(BindingHypothesis hyp : hyps) { 
			addChild(hyp);
		}
	}
	
	public void addChild(BindingHypothesis hyp) {
		if(!hasTestedHypothesis(hyp)) { 
			if(children == null) { 
				children = new LinkedHashMap<BindingHypothesis,HypothesisTree>();
			}
			
			if(!children.containsKey(hyp)) {
                
                HypothesisTree newTree = new HypothesisTree(explorer, hyp, inconsistentRegions, this);
                
                if(children.size() >= branching) {
                    BindingHypothesis worst = findWorstChild();
                    
                    if(children.get(worst).getScoredHypothesis().getInconsistentScore() >
                        newTree.getScoredHypothesis().getInconsistentScore()) { 

                        System.out.println("\t* Pruning " + worst + ", replaced with " + hyp);
                        children.remove(worst);
                        children.put(hyp, newTree);
                    }
                } else { 
                    System.out.println("\t* Adding " + hyp);
                    children.put(hyp, newTree);
                }
			}
		}
	}
    
    private BindingHypothesis findWorstChild() { 
        if(children == null) { return null; }
        int maxScore = -1;
        BindingHypothesis worst = null;
        for(BindingHypothesis hyp : children.keySet()) { 
            ScoredHypothesis shyp = children.get(hyp).getScoredHypothesis();
            if(shyp.getInconsistentScore() >= maxScore) {
                maxScore = shyp.getInconsistentScore();
                worst = hyp;
            }
        }
        return worst;
    }
    
	public boolean hasTestedHypothesis(BindingHypothesis hyp) {
		return hypothesis != null && 
			(hypothesis.getHypothesis().equals(hyp) || 
			 parent.hasTestedHypothesis(hyp));
	}
	
	public ScoredHypothesis getScoredHypothesis() { return hypothesis; }
	public HypothesisTree getParent() { return parent; }
	
	public int findOptimalScore() { 
		if(children == null) { 
			return hypothesis != null ? hypothesis.getInconsistentScore() :
				testableRegions.size();
		} else { 
			int minScore = inconsistentRegions.size()+1;
			for(BindingHypothesis hyp : children.keySet()) { 
				int score = children.get(hyp).findOptimalScore();
				if(score < minScore) { 
					minScore = score;
				}
			}
			return minScore;
		}
	}
	
	public LinkedList<BindingHypothesis> findOptimalPath() { 
		if(children == null) { 
			LinkedList<BindingHypothesis> path = new LinkedList<BindingHypothesis>();
			if(hypothesis != null) { 
				path.addLast(hypothesis.getHypothesis());
			}
			return path;
		} else { 
			int minScore = inconsistentRegions.size()+1;
			BindingHypothesis mhyp = null;
			for(BindingHypothesis hyp : children.keySet()) { 
				int score = children.get(hyp).findOptimalScore();
				if(score < minScore) { 
					minScore = score;
					mhyp = hyp;
				}
			}
			
			LinkedList<BindingHypothesis> minpath = 
				children.get(mhyp).findOptimalPath();
			if(hypothesis != null) { 
				minpath.addFirst(hypothesis.getHypothesis());
			}
			return minpath;
		}
	}
}
