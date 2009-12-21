/*
 * Created on Jan 22, 2008
 */
package edu.mit.csail.cgs.utils.strings;

import java.io.PrintStream;
import java.util.*;
import java.util.logging.*;

/**
 * @author Timothy Danford
 * 
 * An implementation of suffix trees, and of Ukkonen's Algorithm for 
 * building a suffix tree in linear time, adapted from Chapter 6 of 
 * Gusfield, "Algorithms on Strings, Trees, and Sequences."  
 */
public class UkkonenSuffixTree {
	
	public static void main(String[] args) { 
        String target = args[0];
        String query = args[1];
		UkkonenSuffixTree suffixTree = new UkkonenSuffixTree();
		
		suffixTree.addString(target);
        
        System.out.println("Test1 Tree:");
        suffixTree.print(System.out);
        System.out.println(); System.out.flush();

        System.err.println();
        System.err.flush();
        
        assert suffixTree.check();
        
        Set<StringSuffix> matches = suffixTree.matchString(query);
        System.out.println("Full matches are " + matches);
        matches = suffixTree.matchStringPartial(query);
        System.out.println("Partial matches are " + matches);

        assert suffixTree.check();
	}

    private Vector<TreeString> strings;
    private Vector<TreeEdge> totalStringEdges;
    private TreeNode root;
    private char terminal;
    
    private UkkonenState extState;  // this is a field, because TreeEdge depends on it.
    
    private Logger logger;
    private boolean isLogging;
    private Level minLogLevel;
    
    public UkkonenSuffixTree() {
    	logger = Logger.getLogger("edu.mit.csail.cgs.projects.chipseq.assembler.UkkonenSuffixTree");
    	isLogging = true;
    	minLogLevel = Level.SEVERE;
    	logger.setFilter(new LoggingFilter());
    	logger.addHandler(new LoggingHandler(System.err, false));
    	logger.setUseParentHandlers(false);
    	logger.setLevel(Level.SEVERE);
    	
    	logger.log(Level.FINE, "Logger setup complete.");
    	
        strings = new Vector<TreeString>();
        terminal = '$';
        root = new TreeNode(null);
        totalStringEdges = new Vector<TreeEdge>();
        extState = null;
    }
    
    public TreeString getString(int i) { return strings.get(i); }
    public int size() { return strings.size(); }
    public boolean isTerminal(char c) { return c==terminal; }
    
    public void print(PrintStream ps) { 
    	root.print(0, ps);
    }
    
    public void addString(String str) {
        char[] array = str.toCharArray();
        strings.add(new TreeString(strings.size(), array));
        
        logger.log(Level.INFO, String.format("Adding string \"%s\"", str));

        ukkonenExtendSuffixTree(strings.size()-1);
        //naiveExtendSuffixTree(strings.size()-1);
    }
    
    public Set<StringSuffix> matchString(String str) {  
        char[] array = str.toCharArray();
        EdgeMatch m = findEdge(root, array, 0, array.length, false);
        
        if(m.completedMatch()) { 
            return collectSuffixes(m.lastEdge.tailNode);
        } else { 
            return new TreeSet<StringSuffix>();
        }
    }

    public Set<StringSuffix> matchStringPartial(String str) {  
        char[] array = str.toCharArray();
        EdgeMatch m = findEdge(root, array, 0, array.length, false);

        if (m.lastEdge != null) {
            return collectSuffixes(m.lastEdge.tailNode);
        } else {
            return new TreeSet<StringSuffix>();
        }
    }
    
    public boolean check() { 
        return root.check();
    }

    /** Internal Methods ************************************************************/
    
    private EdgeMatch findEdge(TreeNode currentNode, char[] array,  
    		int start, int end, boolean skipcount) {
    	EdgeMatch em = new EdgeMatch(array, start, end);
    	em.matchFrom(currentNode, skipcount);
    	return em;
    }
    
    private EdgeMatch findEdge(TreeNode currentNode, TreeString string,  
    		int start, int end, boolean skipcount) {
    	EdgeMatch em = new EdgeMatch(string, start, end);
    	em.matchFrom(currentNode, skipcount);
    	return em;
    }
    
    private void naiveExtendSuffixTree(int arrayIdx) {
    	TreeString string = strings.get(arrayIdx);
    	
    	// the array.length-1 constraint, instead of array.length, is because
    	// we assume that the terminal character has already been added to the
    	// string, and we don't want to *just* add the suffix that is that 
    	// character.
    	for(int i = 0; i <= string.length(); i++) {
    		logger.log(Level.FINEST, String.format("Naive Extension: \"%s\"",
    				string.substring(i, string.length()+1)));
    		
    		naiveExtendSuffix(string, i);
    	}
    }
    
    private void naiveExtendSuffix(TreeString string, int start) {
    	EdgeMatch em = findEdge(root, string, start, string.length(), false);
    	StringSuffix stringSuffix = new StringSuffix(string, start);
    	
    	TreeEdge leafEdge = null;
    	if(em.completedMatch()) {
    		leafEdge = em.lastEdge;
    	} else { 
    		if(em.lastEdge == null) { 
    			leafEdge = new TreeEdge(string, start, string.length(), root);
    			root.addEdge(leafEdge);
    		} else { 
				leafEdge = new TreeEdge(string, em.matchedTo, string.length(), em.lastEdge.tailNode);
    			if(em.inEdgeMiddle()) { 
    				int offset = em.lastMatchLength();
    				em.lastEdge.split(offset);
    			}
				em.lastEdge.tailNode.addEdge(leafEdge);
    		}
    	}
    	
		leafEdge.tailNode.suffixes.add(stringSuffix);
    }
    
    private void ukkonenExtendSuffixTree(int arrayIdx) {
    	logger.entering("UkkonenSuffixTree", "ukkonenExtendSuffixTree");
    	logger.log(Level.FINEST, String.format("Ukkonen Algorithm String #%d", arrayIdx));
    	
    	TreeString string = strings.get(arrayIdx);
    	extState = new UkkonenState(string);
    	
    	logger.log(Level.FINEST, String.format(
    			"Ukkonen: (%d,%d)", extState.nextPhaseStart, extState.string.length()));
    	
        for(int phase = extState.nextPhaseStart; phase < extState.string.length(); phase++) {
        	ukkonenSPA(phase);
        	
            System.err.println(String.format("Phase %d results: ", phase));
            print(System.err); System.err.println(); System.err.flush();
        }
        
        logger.log(Level.FINEST, String.format("Finishing edges: %d", extState.lastE));
        extState.finishFinalEdges();

        System.err.println(String.format("Finished results: "));
        print(System.err); System.err.println(); System.err.flush();

    	logger.exiting("UkkonenSuffixTree", "ukkonenExtendSuffixTree");
    }

    /* ukkonenSPA(i) performs phase i of Ukkonen's algorithm.  This 
     * means that we're making sure that array[0,i] (note the inclusivity!) 
     * is a part of the current suffix tree.
     * 
     * Original Description: pg. 106 of Gusfield
     */
    private void ukkonenSPA(int i) {
    	logger.entering("UkkonenSuffixTree", "ukkonenSPA");
        logger.log(Level.FINEST, String.format("i=%d", i));
        
        assert i >= 0;
    	
        /*
         * SPA Step 1:
         * "Increment index e to i+1"
         *  
         * The equivalent of Gusfield's i+1 is, in our situation, just i. 
         * However, the coordinates are inclusive in Gusfield, 
         * and exclusive in our case (along the tree edges).  Therefore, 
         * lastE should be updated to be i+1, exactly. 
         */
        extState.lastE = i+1;
        logger.log(Level.FINEST, String.format("e=%d", extState.lastE));
        
        /*
         * SPA Step 2:
         * "Explicitly compute successive extensions, using the SEA algorithm, 
         * starting at j_i + 1 until reaching the first extension j* where rule3
         * applies or until all extensions are done in this phase."
         *  
         * extState.nextExtStart encodes the (j_i)+1 value.  We start there, and 
         * iterate forward until all extensions have been performed, or until 
         * ukkonenSEA returns false (ukkonenSEA returns a true if rule 1 or rule 2 
         * applies in its extension).
         * 
         * We extend until j==i, because the last extension of each phase is 
         * the extension that *just* adds the new character into the tree.
         */

        logger.log(Level.FINEST, String.format("jstart=%d", extState.nextExtStart));
        boolean keepExtending = true;
        int j = extState.nextExtStart;
        
        while(keepExtending && j <= i) {
        	if(ukkonenSEA(i, j)) { 
        		j++;
                
                // we don't want to just put in the terminal character.
                if(i == extState.string.length()-1 && j == i) { 
                    keepExtending = false;
                }
        	} else { 
        		keepExtending = false;
        	}
        	
        	System.out.println(String.format("Phase %d, Extension %d tree: ", i, j));
        	print(System.out);
        	System.out.println(); System.out.flush(); 
        	System.err.println(); System.err.flush();
        }

        /*
         * SPA Step 3:
         * "Set j_{i+1} to j*-1, to prepare for the next phase."
         */
        extState.nextExtStart = j;
        
        logger.log(Level.FINEST, String.format("j*=%d", extState.nextExtStart));
    	logger.exiting("UkkonenSuffixTree", "ukkonenSPA");
    }
    
    /*
     * ukkonenSEA(i, j) performs extension j of phase i of Ukkonen's algorithm.
     * This means that we're making sure that array[j,i] (note the inclusivity!) 
     * is a part of the current suffix tree.
     * 
     * Original Description: pg. 100 of Gusfield
     */
    private boolean ukkonenSEA(int i, int j) { 
    	logger.exiting("UkkonenSuffixTree", "ukkonenSEA");
        logger.log(Level.FINEST, String.format("j=%d", j));
        
        assert j <= i;

    	boolean rule3 = false;
        
        TreeNode newRule2Node = null;
        
        EdgeMatch m = extState.matcher;
        char lastChar = extState.string.getChar(i);
        boolean lastCharIsTerminal = isTerminal(lastChar);

        /*
         * SEA Step 1:  
         * "Find the first node v at or above the end of S[j-1,i] that either
         * has a suffix link from it or is the root.  This requires walking up 
         * at most one edge from the end of S[j-1,i] in the current tree.  Let
         * \gamma (possibly empty) denote the string between v and the 
         * end of S[j-1,i]." 
         */
        
        /*
         * SEA Step 2:
         * "If v is not the root, traverse the suffix link from v to node
         * s(v) and then walk down from s(v) following the path for string
         * gamma.  If v is the root, then follow the path for S[j,i] from the 
         * root (as in the naive algorithm)." 
         */

        int gammaEnd = i;
        int gammaStart = gammaEnd - extState.gammaLength;
        
        if(extState.nextNode == null || extState.nextNode.isRoot()) {
        	String beta = extState.string.substring(j,i);
            logger.log(Level.FINEST, String.format("beta: %d,%d <%s>%c", j, i, beta, lastChar));

            m.reset(j, i);
            m.matchFrom(root, true);
        } else { 
            logger.log(Level.FINEST, String.format("gammaLength:%d", extState.gammaLength));
            String gamma = extState.string.substring(gammaStart, gammaEnd);
            logger.log(Level.FINEST, String.format("gamma: %d,%d <%s>%c", gammaStart, gammaEnd, gamma, lastChar));

            m.reset(gammaStart, gammaEnd);
            m.matchFrom(extState.nextNode, true);
        }

        /*
         * SEA Step 3:
         * "Using the extension rules, ensure that the string S[j,i]S(i+1) is 
         * in the tree."
         * 
         * In our coordinates, this is array[j,i)+array[i]
         * \beta = array[j,i)
         * 
         * Rule 1: the path \beta ends at a leaf.  (we shouldn't see this case).
         * Rule 2: the path \beta is not continued by array[i].  That is, \beta
         *         ends either at a node (in which case, no child of the node 
         *         starts with array[i]), or in an edge (in which case, the edge 
         *         doesn't continue with array[i]). Either way, we create a new
         *         edge that is labeled with array[i] (coordinates: [i,i+1) ).
         * Rule 3: \beta+array[i] is already in the tree -- either \beta ends in 
         *         an edge that continues with array[i], or at a node that has 
         *         a child under array[i].  Either way, return false (break!).
         */
        
        TreeEdge newEdge = null;

        if(m.lastEdge == null) {
        	logger.log(Level.FINEST, String.format("Found root."));
        	// the \beta string matched to the root (was empty).  So we need
        	// to simply check the children of the root.
        	
        	boolean foundLastChar = !lastCharIsTerminal ? 
        			root.childEdges.containsKey(lastChar) : 
        			root.terminalEdges.containsKey(extState.string.getIndex());
        	
        	if(foundLastChar) { 
        		// Rule 3
        		rule3 = true;
        		logger.log(Level.FINEST, "Rule #3, Root");
        		
        		extState.nextNode = null;
        		extState.gammaLength = 0;
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
        	} else { 
        		// Rule 2
        		logger.log(Level.FINEST, "Rule #2, Root");
        		newEdge = new TreeEdge(extState.string, i, null, root);
        		root.addEdge(newEdge);

        		extState.nextNode = null;
        		extState.gammaLength = 0;
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
        	}

        } else if(m.inEdgeMiddle()) { 
        	int offset = m.lastMatchLength();
        	logger.log(Level.FINEST, String.format("Found edge middle: %d", offset));
        	
        	boolean foundLastChar = !lastCharIsTerminal ? 
        			m.lastEdge.getChar(offset) == lastChar : 
        			(m.lastEdge.string.getIndex() == extState.string.getIndex() && 
        			 offset == m.lastEdge.length()-1);
        	
        	logger.log(Level.FINEST, String.format("foundLastChar: %s", foundLastChar));

        	if(foundLastChar) { 
        		// Rule 3
        		rule3 = true;
        		logger.log(Level.FINEST, "Rule #3, Edge");

        		extState.nextNode = m.lastEdge.headNode;
        		//extState.gammaLength = m.lastMatchLength() + 1;
        		extState.gammaLength = m.lastMatchLength() + (j==i ? 1 : 0);
                
        		assert extState.gammaLength >= 0;
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
        	} else { 
        		// Rule 2
        		logger.log(Level.FINEST, "Rule #2, Edge");
                
        		TreeEdge newLowerEdge = m.lastEdge.split(offset);
                extState.edgesWithE.add(newLowerEdge);
                
        		newEdge = new TreeEdge(extState.string, i, null, m.lastEdge.tailNode);
        		m.lastEdge.tailNode.addEdge(newEdge);
        		newRule2Node = m.lastEdge.tailNode;

        		extState.nextNode = m.lastEdge.headNode;
        		//extState.gammaLength = m.lastEdge.length() + 1;
        		extState.gammaLength = m.lastEdge.length() + (j==i ? 1 : 0);
        		
                assert extState.gammaLength >= 0;
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
        	}

            if(extState.nextNode.suffixLink == null && !extState.nextNode.isRoot()) {
            	logger.log(Level.FINEST, String.format("Walking up edge: %d", 
            			extState.nextNode.parentEdge.length()));
                extState.gammaLength += extState.nextNode.parentEdge.length();
                extState.nextNode = extState.nextNode.parentEdge.headNode;
            }

        } else { 
        	logger.log(Level.FINEST, String.format("Found node."));
        	
        	boolean foundLastChar = !lastCharIsTerminal ?  
        			m.lastEdge.tailNode.childEdges.containsKey(lastChar) : 
        			m.lastEdge.tailNode.terminalEdges.containsKey(extState.string.getIndex());

        	logger.log(Level.FINEST, String.format("foundLastChar: %s", foundLastChar));

        	if(foundLastChar) { 
        		// Rule 3
        		rule3 = true;
        		logger.log(Level.FINEST, "Rule #3, Node");

        		extState.nextNode = m.lastEdge.headNode;
        		//extState.gammaLength = m.lastEdge.length() + 1;
        		extState.gammaLength = m.lastEdge.length() + (j==i ? 1 : 0);
                
        		assert extState.gammaLength >= 0;
                
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
        	} else { 
        		// Rule 2
        		logger.log(Level.FINEST, "Rule #2, Node");
        		newEdge = new TreeEdge(extState.string, i, null, m.lastEdge.tailNode);
        		m.lastEdge.tailNode.addEdge(newEdge);

        		extState.nextNode = m.lastEdge.headNode;
        		//extState.gammaLength = m.lastEdge.length() + 1;
        		extState.gammaLength = m.lastEdge.length() + (j==i ? 1 : 0);
        		
        		logger.log(Level.FINEST, String.format("nextGamma: %d", extState.gammaLength));
                assert extState.gammaLength >= 0;
        	}
            
            if(extState.nextNode.suffixLink == null && !extState.nextNode.isRoot()) {
            	logger.log(Level.FINEST, String.format("Walking up edge: %d", 
            			extState.nextNode.parentEdge.length()));
                extState.gammaLength += extState.nextNode.parentEdge.length();
                extState.nextNode = extState.nextNode.parentEdge.headNode;
            }
        }
        
        if(extState.nextNode != null) {
        	logger.log(Level.FINEST, "Following suffix link.");
            extState.nextNode = extState.nextNode.suffixLink;
        } else { 
        	logger.log(Level.FINEST, "Suffix link not found.");
        }
        
        if(newEdge != null) { 
            newEdge.tailNode.suffixes.add(extState.currentSuffix);
            extState.nextSuffix();
            
        	extState.edgesWithE.add(newEdge);
        	logger.log(Level.FINEST, String.format("Added suffix: %d", j));
        }
        
        /*
         * SEA Step 4:
         * "If a new internal node w was created in extension j-1 (by extension rule 2) 
         * then by Lemma 6.1.1 string alpha must end at node s(w), the end node for the 
         * suffix link from w.  Create the suffix link (w, s(w)) from w to s(w)."
         * 
         * This wording is confusing -- is there a typo in Gusfield?  I'm not sure where
         * the 'w' comes from.  
         */
        if(extState.rule2Node != null) {
        	if(m.lastEdge != null) { 
        		extState.rule2Node.suffixLink = m.lastEdge.tailNode;
            	logger.log(Level.FINEST, "Adding suffix link --> internal node.");
        	} else {
        		extState.rule2Node.suffixLink = root;
            	logger.log(Level.FINEST, "Adding suffix link --> root.");
        	}
        }

        /*
         * Update any state that will be needed in the next extension.
         */
        extState.rule2Node = newRule2Node;
        
    	logger.exiting("UkkonenSuffixTree", "ukkonenSEA");
    	
    	// "Rule 3 is a show stopper" means that, if we encounter rule 3, 
        // we *don't* continue.
        return !rule3;
    }
    
    private Set<StringSuffix> collectSuffixes(TreeNode tn) { 
        TreeSet<StringSuffix> set = new TreeSet<StringSuffix>();
        tn.collectSuffixes(set);
        return set;
    }

    /** Helper Classes **************************************************************/
    
    public class TreeString {
    	
    	private int index;
    	private char[] array;
    	
    	private TreeString(int idx, char[] a) { 
    		index = idx;
    		array = a;
    	}
    	
    	public int getIndex() { return index; }
    	public int length() { return array.length+1; }
    	public char getChar(int i) { return i < array.length ? array[i] : terminal; }
    	
    	public boolean matches(int offset, TreeString str, int strOffset) { 
			assert str != null;
			assert offset >= 0;
			assert offset <= array.length;
			assert strOffset >= 0;
			assert strOffset <= str.array.length;

    		if(offset==array.length || strOffset == str.array.length) { 
    			return index==str.index && offset==array.length && strOffset==str.array.length;
    		} else { 
    			return array[offset] == str.array[strOffset];
    		}
    	}
    	
    	public boolean matches(int offset, char[] str, int strOffset) { 
    		if(offset == array.length || strOffset >= str.length) { 
    			return false;
    		} else { 
    			return array[offset]==str[strOffset];
    		}
    	}
    	
    	public String substring(int start, int end) { 
    		StringBuilder sb = new StringBuilder();
    		for(int i = start; i < end; i++) { 
    			sb.append(getChar(i));
    		}
    		return sb.toString();
    	}
    	
    	public int hashCode() { 
    		int code = 17;
    		code += index; code *= 37;
    		return code;
    	}
    	
    	public boolean equals(Object o) { 
    		if(!(o instanceof TreeString)) { return false; }
    		TreeString ts = (TreeString)o;
    		return ts.index == index;
    	}
    	
    	public String toString() { 
    		return String.format("#%d:%s", index, new String(array));
    	}
    }

    public class StringSuffix implements Comparable<StringSuffix> { 

    	private int offset;
    	private TreeString string;
        
        private StringSuffix(TreeString ts, int off) {
        	string = ts;
            offset = off;
        }
        
        public int getStringIndex() { return string.getIndex(); }
        public int getOffset() { return offset; }
        
        public String getSuffixString() { 
            StringBuilder sb = new StringBuilder();
            for(int i = offset; i <= string.length(); i++) { 
                sb.append(string.getChar(i));
            }
            return sb.toString();
        }
        
        public String toString() { 
        	return String.format("(#%d,+%d)", getStringIndex(), offset); 
        }
        
        public int hashCode() { 
            int code = 17;
            code += string.hashCode(); code *= 37;
            code += offset; code *= 37;
            return code;
        }
        
        public int compareTo(StringSuffix ss) {
        	int stringID = getStringIndex();
            if(stringID < ss.getStringIndex()) { return -1; }
            if(stringID > ss.getStringIndex()) { return 1; }
            if(offset < ss.offset) { return -1; }
            if(offset > ss.offset) { return 1; }
            return 0;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof StringSuffix)) { return false; }
            StringSuffix ss = (StringSuffix)o;
            if(!ss.string.equals(string)) { return false; }
            return offset==ss.offset;
        }
    }    

    /** Internal Classes ************************************************************/
    
    private class TreeNode { 

        public TreeEdge parentEdge;
        public TreeNode suffixLink;
        public Set<StringSuffix> suffixes;

        private Map<Character,TreeEdge> childEdges;
        private Map<Integer,TreeEdge> terminalEdges;

        public TreeNode(TreeEdge p) { 
            parentEdge = p;
            suffixLink = null;
            suffixes = new TreeSet<StringSuffix>();

            childEdges = new TreeMap<Character,TreeEdge>();
            terminalEdges = new TreeMap<Integer,TreeEdge>();
        }
        
        // Walks the tree, checking to make sure that we haven't violated any 
        // logical constraints.  This is a method for debugging. 
        public boolean check() { 
            String path = pathLabel();
            if(suffixLink != null) { 
                String suffixPath = suffixLink.pathLabel();
                if(!path.substring(1, path.length()).equals(suffixPath)) {
                    logger.log(Level.SEVERE, String.format("Suffix Link for node (%s) didn't match: %s", path, suffixPath));
                    return false;
                }
            }
            
            for(char c : childEdges.keySet()) { 
                TreeEdge e = childEdges.get(c);
                if(!e.check()) { 
                    return false;
                }
            }
            
            for(int k : terminalEdges.keySet()) { 
            	TreeEdge e = terminalEdges.get(k);
            	if(!e.check()) { 
            		return false;
            	}
            }
            
            return true;
        }
        
        public String pathLabel() { 
            if(parentEdge != null) { 
                StringBuilder sb = new StringBuilder(parentEdge.headNode.pathLabel());
                for(int i = 0; i < parentEdge.length(); i++) { 
                    sb.append(parentEdge.getChar(i));
                }
                return sb.toString();
            } else { 
                return "";
            }
        }
        
        public void print(int indent, PrintStream ps) { 
        	if(isLeaf()) {
        		printSuffixes(ps);
        		ps.println();
        	} else { 
        		int i = 0;
        		for(char c : childEdges.keySet()) { 
        			TreeEdge edge = childEdges.get(c);
        			edge.print(indent, i!=0, ps);
        			i++;
        		}
        		for(int k : terminalEdges.keySet()) { 
        			TreeEdge edge = terminalEdges.get(k);
        			edge.print(indent, i!=0, ps);
        			i++;
        		}
        	}
        }
        
        public void printSuffixes(PrintStream ps) { 
    		ps.print(" [");
    		int i = 0;
    		for(StringSuffix ss : suffixes) {
    			ps.print((i == 0 ? "" : ",") + ss.toString());
    			i++;
    		}
    		ps.print("]");        	
        }
        
        public boolean isRoot() { return parentEdge == null; }
        public boolean isLeaf() { return childEdges.isEmpty(); }
        
        public void collectSuffixes(Set<StringSuffix> suffices) { 
            suffices.addAll(suffixes);
            for(char key : childEdges.keySet()) { 
                childEdges.get(key).tailNode.collectSuffixes(suffices);
            }
            for(int key : terminalEdges.keySet()) { 
            	terminalEdges.get(key).tailNode.collectSuffixes(suffices);
            }
        }
        
        public void addEdge(TreeEdge e) {
        	if(e.start() == e.string.length()-1) { 
        		int key = e.string.getIndex();
        		if(terminalEdges.containsKey(key)) { throw new IllegalArgumentException(); }
        		
        		e.headNode = this;
        		terminalEdges.put(key, e);
        	} else {
        		char key = e.getChar(0);
                if(childEdges.containsKey(key)) { throw new IllegalArgumentException(); }
                e.headNode = this;
                childEdges.put(key, e);
        	}
        }
    }

    /*
     * This is the state that we need to carry through, from extension to 
     * extension and from phase to phase, of the Ukkonen extension algo-
     * rithm.  Particularly important is the "lastE" field, which is 
     * globally referenced by a certain subset of the TreeEdges (the 'leaf'
     * edges) during the execution of the Ukkonen Algorithm.  As a result,
     * there is a global field ('extState') that is maintained during any 
     * run of the algorithm.
     */
    private class UkkonenState {

    	public TreeString string;
    	
        // this is set for the very first phase, when we're 
        // adding a string that shares a prefix with a string
        // already in the tree.
        public int nextPhaseStart;
        
        // these are updated once a phase.
        public int lastE;
        public int nextExtStart;

        // these are updated each extension (potentially).
        public EdgeMatch matcher;
        public TreeNode rule2Node;
        public TreeNode nextNode;
        public int gammaLength;
        public LinkedList<TreeEdge> edgesWithE;
        public StringSuffix currentSuffix;
        
        public UkkonenState(TreeString str) {
        	string = str;
        	
            lastE = 0;
            edgesWithE = new LinkedList<TreeEdge>();
            nextPhaseStart = 0;
            nextExtStart = 0;
            matcher = null;
            nextNode = root;
            gammaLength = 0;
            rule2Node = null;
            
            if(string.getIndex() > 0) { 
                matcher = findEdge(root, string, 0, string.length(), false);
                nextPhaseStart = matcher.matchedTo;
                nextExtStart = 0;
                lastE = matcher.matchedTo;
                
                logger.log(Level.FINEST, String.format(
                        "String %s can start at phase %d (E:%d)", 
                        string.toString(), nextPhaseStart, lastE));
            } else { 
            	matcher = new EdgeMatch(string, 0, string.length());
            }

            currentSuffix = new StringSuffix(string, 0);
        }
        
        public void nextSuffix() { 
            currentSuffix = new StringSuffix(string, currentSuffix.getOffset() + 1);
        }
        
        public void finishFinalEdges() {
            for(TreeEdge edge : edgesWithE) {
                if(edge.isEndSymbolic()) { 
                    edge.setEnd(lastE);
                }
            }
            edgesWithE.clear();
        }
    }

    /*
     * We can't just match into the tree in a naive way -- we need to 
     * store some state about how *far* we matched into the tree, where
     * we ended up, etc.   As a result, this EdgeMatch class implements
     * tree descent -- it walks as far down into the tree as possible 
     * and then (upon return from the matchFrom() method) contains all the 
     * state we need about where that match ended.
     */
    private class EdgeMatch {

    	// only one of these is non-null.  
    	public char[] array;
    	public TreeString string;
    	
        public int start, end;
        
        public TreeEdge lastEdge;
        public int matchingFrom, matchedTo;
        
        public EdgeMatch(char[] a, int st, int ed) { 
        	array = a;
        	string = null;
			assert string != null || array != null;
        	
        	start = st;
        	end = ed;

        	lastEdge = null;
        	matchingFrom = matchedTo = -1;
        }
        
        public EdgeMatch(TreeString s, int st, int ed) { 
        	array = null;
        	string = s;
			assert string != null || array != null;
        	
        	start = st;
        	end = ed;

        	lastEdge = null;
        	matchingFrom = matchedTo = -1;
        }
        
        public char getChar(int i) { 
        	return array != null ? array[i] : string.getChar(i);
        }
        
        public void reset(int st, int ed) { 
            start = st;
            end = ed;
            lastEdge = null;
            matchingFrom = matchedTo = -1;
        }
        
        public String matchingString() {
        	if(array != null) { 
        		return new String(array, start, end-start);
        	} else { 
        		return string.substring(start,end);
        	}
        }
        
        public String matchedString() {
        	if(array != null) { 
        		return new String(array, start, matchedTo-start);
        	} else { 
        		return string.substring(start,matchedTo);
        	}
        }
        
        public String currentMatchString() {
        	if(array != null) { 
        		return new String(array, matchingFrom, end-matchingFrom);
        	} else { 
        		return string.substring(matchingFrom,end);
        	}
        }
        
        public String toString() { 
        	return String.format("EdgeMatch (%d,%d-->%d,%d) '%s' in '%s'", 
                    start, matchingFrom, matchedTo, end, matchedString(), matchingString());
        }

        public void matchFrom(TreeNode startNode, boolean skipcount) {
			assert string != null || array != null;
			assert lastEdge == null;

        	matchingFrom = start;
        	matchedTo = start;
        	
        	// a base case.
        	if(end-start <= 0) { return; }

        	char nextChar = getChar(matchingFrom);

			assert startNode != null;

            System.err.println("startNode is " + startNode + " and nextChar is " + nextChar);
        	
        	if(startNode.childEdges.containsKey(nextChar)) { 
        		lastEdge = startNode.childEdges.get(nextChar);
        	} else if(skipcount) {
        		System.err.println("Failure Node: ");
        		startNode.print(0, System.err);
        		System.err.println(); System.err.flush();
        		throw new IllegalArgumentException();
        	} 
        	
        	boolean keepMatching = lastEdge != null;
        	
        	while(keepMatching) {
        		int remaining = end-matchingFrom;
        		int matching = skipcount ? 
        				Math.min(lastEdge.length(), remaining) :
						(string != null ? lastEdge.countMatches(string, matchingFrom, end) : 
										  lastEdge.countMatches(array, matchingFrom, end));
        		matchedTo = matchingFrom + matching;
        		
        		if(matching < lastEdge.length()) {
        			// we either matched all the way into the middle of the 
        			// current edge, or we found a mismatch in the current edge.
        			// either way, update the matchedTo variable and we're done.
        			keepMatching = false;
        		} else { 
        			if(matching < remaining) {
        				nextChar = getChar(matchedTo);
        				
        				if(lastEdge.tailNode.childEdges.containsKey(nextChar)) {
        					lastEdge = lastEdge.tailNode.childEdges.get(nextChar);
        					matchingFrom += matching;
        					
        				} else if(skipcount) {  
                            System.err.println("ERROR TREE: ");
                            print(System.err); System.err.println(); System.err.flush();
                            
                            String err = String.format("[%s] node:'%s' (next: %c)", toString(), startNode.pathLabel(), nextChar);
        					throw new IllegalArgumentException(err);
        				} else { 
        					keepMatching = false;
        				}
        			} else { 
        				keepMatching = false;
        			}
        		}
        	}
        }
        
        public int lastMatchLength() { return matchedTo - matchingFrom; }
        public int matchLength() { return matchedTo - start; }
        
        public boolean inEdgeMiddle() { 
        	return lastEdge != null && lastMatchLength() < lastEdge.length();
        }
        
        public boolean completedMatch() { return lastEdge != null && matchedTo == end; }
    }
    
    private class TreeEdge { 

		// This is a pointer to the overall string of which this edge is a substring.
    	public TreeString string;  

		// The nodes between which this edge exists in the tree.
        public TreeNode headNode, tailNode;
        
        // These are tricky -- 'end' can be null, and when it is, it takes
        // the value of extState.lastE.  This is part of the Ukkonen extension
        // algorithm, and as a result if you want to find the coordinates of 
        // a TreeEdge, you should always call start() and end(), the methods.
        private Integer start, end;  // coordinates: [start, end)
        
        public TreeEdge(TreeString str, Integer st, Integer ed, TreeNode h) {
            assert st >= 0;
            assert ed == null ? st < extState.lastE : st < ed;
            
            string = str;
            start = st;
            end = ed;
            headNode = h;
            tailNode = new TreeNode(this);
        }
        
        public TreeEdge(TreeString str, Integer st, Integer ed, TreeNode h, TreeNode t) { 
            assert st >= 0;
            assert ed == null ? st < extState.lastE : st < ed;

            string = str;
            start = st;
            end = ed;
            headNode = h;
            tailNode = t;
        }
        
        public boolean check() { 
            if(end == null) {
                logger.log(Level.SEVERE, String.format("Edge %s still has a [null] for it's end-value.", edgeLabel()));
                return false;
            }
            
            return tailNode.check();
        }
        
        public String edgeLabel() { 
            return String.format("%s+%c", headNode.pathLabel(), getChar(0));
        }
        
        public void print(int indent, boolean printIndent, PrintStream ps) {
        	char nodeChar = headNode.suffixLink != null ? '*' : '-';
        	for(int i = 0; printIndent && i < indent; i++) { ps.print(" "); }
            //String coords = String.format("(%d,%d)", start, end);
            String coords = "";
        	ps.print(String.format("%c%s%s", nodeChar, coords, getSubstring()));
        	tailNode.print(indent+length()+coords.length()+1, ps);
        }
        
        public TreeEdge split(int offset) {
            if(offset <= 0 || offset >= length()) { 
            	throw new IllegalArgumentException(String.format(
            			"Illegal edge split offset %d (length %d)",
            			offset, length()));
            }
            
            TreeNode newInternal = new TreeNode(this);
            TreeEdge newLowerEdge = new TreeEdge(string, start+offset, end, newInternal, tailNode);
            tailNode = newInternal;
            newInternal.addEdge(newLowerEdge);
            
            end = start+offset;
            
            return newLowerEdge;
        }
        
        public char getChar(int offset) {
        	return string.getChar(start+offset);
        }

		public int countMatches(char[] array, int startMatch, int endMatch) { 
        	int ei = start;
        	int mi = startMatch;
        	int eend = end();
        	
        	for(; ei < eend && mi < endMatch; ei++, mi++) {
				assert mi >= 0; 
				assert mi < array.length;

				if(array[mi] != string.getChar(ei)) { 
					return ei-start;
				}
        	}
        	return ei-start;
		}
        
        public int countMatches(TreeString str, int startMatch, int endMatch) {
        	int ei = start;
        	int mi = startMatch;
        	int eend = end();
        	
        	for(; ei < eend && mi < endMatch; ei++, mi++) {
        		if(!string.matches(ei, str, mi)) { 
        			return ei-start;
        		}
        	}
        	return ei-start;
        }
        
        public void setEnd(int e) {
            assert e > start;
            end = e;
        }
        
        public boolean isEndSymbolic() { return end == null; }
        
        public boolean isTerminal() { 
        	return start==string.length();
        }
        
        public int end() { return end != null ? end : extState.lastE; }
        public int start() { return start; }
        public int length() { return end() - start(); }
        
        public String getSubstring() { return string.substring(start(), end()); }
        
        public boolean isLeafEdge() { 
            return tailNode.isLeaf();
        }
    }
    
    private class LoggingHandler extends Handler {
    	
    	private PrintStream logstream;
    	private boolean closeStream;
    	
    	public LoggingHandler(PrintStream ps) { 
    		logstream = ps;
    		closeStream = false;
    	}
    	
    	public LoggingHandler(PrintStream ps, boolean cs) { 
    		logstream = ps;
    		closeStream = cs;
    	}
    	
    	public void setCloseStream(boolean cs) { closeStream = cs; }

    	public void close() throws SecurityException {
    		if(closeStream) { 
    			logstream.close();
    		}
		}

		public void flush() {
			logstream.flush();
		}

		public void publish(LogRecord rec) {
			String msg = String.format("%s: %s", rec.getLevel().toString(), rec.getMessage());
			logstream.println(msg);
		}  
    }

    /*
     * The Level.intValue() method apparently doesn't do what we want?  So 
     * I'll write my own.  This lets us filter by "greater than" relations
     * on Level's, with the natural (to me) ordering.
     */
    private class LoggingFilter implements Filter {
    	
    	public LoggingFilter() { 
    		
    	}
    	
    	public boolean isLoggable(LogRecord rec) { 
    		return isLogging && logValue(rec.getLevel()) >= logValue(minLogLevel);
    	}
    	
    	private int logValue(Level l) { 
    		if(l.equals(Level.FINEST)) { 
    			return 0;
    		} else if(l.equals(Level.FINE)) { 
    			return 1; 
    		} else if(l.equals(Level.INFO)) { 
    			return 2;
    		} else if(l.equals(Level.WARNING)) { 
    			return 3;
    		} else if(l.equals(Level.SEVERE)) { 
    			return 4;
    		}
    		return -1;
    	}
    }
}

