/*
 * Created on Sep 7, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.BIND;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class BindInteractionSet {
    
    public static void main(String[] args) { 
        try { 
            File ints = new File(args[0]), refs = new File(args[1]), labels = new File(args[2]);
            File names = new File(args[3]);
			int taxon = Integer.parseInt(args[4]);
            printPairwiseInteractions(ints, refs, labels, names, taxon);
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
    }
    
    public static void printPairwiseInteractions(File ints, File refs, File labels, File nameFile, int taxon)
        throws IOException { 
        
        BindInteractionSet bis = loadFullInteractionSet(ints, refs, labels);
        LineFile namef = new LineFile(nameFile);
        Vector<String> names = new Vector<String>(namef.getAllLines());

        Collection<Pair<String,String>> pairs = bis.findConnectedNames(names, taxon);
        for(Pair<String,String> p : pairs) { 
            System.out.println(p.getFirst() + "\t" + p.getLast());
        }
    }

	public static Collection<BindInteractor> filterByType(Collection<BindInteractor> original, 
		String type) { 
		
		LinkedList<BindInteractor> lst = new LinkedList<BindInteractor>();
		for(BindInteractor bi : original) { 
			if(bi.getType().equals(type)) { 
				lst.addLast(bi);
			}
		}
		return lst;
	}

	public static Collection<BindInteractor> filterByTaxon(Collection<BindInteractor> original, int taxon) { 
		
		LinkedList<BindInteractor> lst = new LinkedList<BindInteractor>();
		for(BindInteractor bi : original) { 
			if(bi.getTaxon() == taxon) { 
				lst.addLast(bi);
			}
		}
		return lst;
	}

    public static BindInteractionSet loadFullInteractionSet(File ints, File refs, File labels) 
        throws IOException {
        BindInteractionSet bis = new BindInteractionSet(ints);
        bis.loadInteractorLabels(labels);
        bis.loadInteractionReferences(refs);
        return bis;
    }
    
    private Map<Integer,BindInteraction> interactions;
    private Map<String,BindInteractor> interactors;
    private Map<BindInteractor,Set<BindInteraction>> interactor2interactions;

    public BindInteractionSet(File f) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        LinkedList<FlatFileRecord> recs = new LinkedList<FlatFileRecord>();
        while((line = br.readLine())!=null) { 
            line = line.trim();
            if(line.length() > 0) { 
                FlatFileRecord rec = new FlatFileRecord(line);
                recs.addLast(rec);
            }
        }
        br.close();

        interactions = new HashMap<Integer,BindInteraction>();
        interactors = new HashMap<String,BindInteractor>();
        interactor2interactions = new HashMap<BindInteractor,Set<BindInteraction>>();
        
        for(FlatFileRecord r : recs) { 
            BindInteractor a = r.createA(), b = r.createB();
            if(!interactors.containsKey(a.getKey())) { interactors.put(a.getKey(), a); }
            if(!interactors.containsKey(b.getKey())) { interactors.put(b.getKey(), b); }
            a = interactors.get(a.getKey());
            b = interactors.get(b.getKey());
            BindInteraction inter = new BindInteraction(r.rgid, r.bid, a, b, r.AB);
            interactions.put(inter.getBID(), inter);
            if(!interactor2interactions.containsKey(a)) { 
                interactor2interactions.put(a, new HashSet<BindInteraction>());
            }
            if(!interactor2interactions.containsKey(b)) { 
                interactor2interactions.put(b, new HashSet<BindInteraction>());
            }
            interactor2interactions.get(a).add(inter);
            interactor2interactions.get(b).add(inter);
        }        
    }
    
    public void loadInteractorLabels(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        while((line = br.readLine()) != null) { 
            if(line.length() > 0) { 
                int start = 0, end = line.indexOf("\t");
                int bid = Integer.parseInt(line.substring(start, end));
                BindInteraction inter = interactions.get(bid);
                BindInteractor a = inter.getA(), b = inter.getB();
				//System.out.println("Naming (" + bid + ") : " + inter);
                
                start = end+1;
                end = line.indexOf("\t", start);
				if(end==-1) { end = line.length(); }
                if(start < line.length() && end - start > 1) {
                    String shortName = line.substring(start, end);
                    a.addShortName(shortName);
					//System.out.println("\tA_Short: " + shortName);
                }

                start = end+1;
                end = line.indexOf("\t", start);
				if(end==-1) { end = line.length(); }
                if(start < line.length() && end - start > 1) {
                    String shortName = line.substring(start, end);
                    b.addShortName(shortName);
					//System.out.println("\tB_Short: " + shortName);
                }

                start = end+1;
                end = line.indexOf("\t", start);
				if(end==-1) { end = line.length(); }
                if(start < line.length() && end-start > 1) { 
                    String[] array = line.substring(start, end).split("\\|");
                    for(int i = 0; i < array.length; i++) { 
                        a.addOtherName(array[i]);
						//System.out.println("\tA_Other: " + array[i]);
                    }
                }

                start = end+1;
                end = line.indexOf("\t", start);
				if(end==-1) { end = line.length(); }
                if(start < line.length() && end-start > 1) { 
                    String[] array = line.substring(start, end).split("\\|");
                    for(int i = 0; i < array.length; i++) { 
                        b.addOtherName(array[i]);
						//System.out.println("\tB_Other: " + array[i]);
                    }
                }
            }
        }
        br.close();
    }
    
    public void loadInteractionReferences(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                String[] array = line.split("\t");
                int rgid = Integer.parseInt(array[0]), bid = Integer.parseInt(array[1]);
                int pmid = Integer.parseInt(array[2]);
				int methodIndex = -1;
				if(array.length >= 4) { 
					String methodName = array[3];
					methodIndex = BindReference.methodName2Index.get(methodName);
				}
                BindReference bref = new BindReference(pmid, methodIndex);
                int interactionKey = bid;
                if(interactions.containsKey(interactionKey)) { 
                    BindInteraction bi = interactions.get(interactionKey);
                    bi.addReference(bref);
                }
            }
        }
        
        br.close();
    }
    
    public Collection<BindInteraction> getAllInteractions() { return interactions.values(); }
    public Collection<BindInteraction> getInteractions(BindInteractor actor) {
        return interactor2interactions.get(actor);
    }
    
    public BindInteraction lookupInteraction(int bid) { 
        return interactions.get(bid);
    }
    
    public Collection<BindInteraction> lookupRedundantInteractions(int rgid) { 
        Set<BindInteraction> lst = new HashSet<BindInteraction>();
        for(int bid : interactions.keySet()) { 
            if(interactions.get(bid).getRGID() == rgid) { 
                lst.add(interactions.get(bid));
            }
        }
        return lst;
    }
    
    public int size() { return interactions.size(); }
    public int getNumInteractions() { return interactions.size(); }
    public int getNumInteractors() { return interactors.size(); }
    
    public boolean isInteractingPair(BindInteractor a, BindInteractor b) { 
        Collection<BindInteraction> inters = interactor2interactions.get(a);
        for(BindInteraction bi : inters) { 
            if(bi.getA().equals(a) && bi.getB().equals(b)) { 
                return true;
            }
        }
        return false;
    }
    
    public Collection<Pair<String,String>> findConnectedNames(Vector<String> names, int taxon) { 
        Map<BindInteractor,String> named = new HashMap<BindInteractor,String>();
        for(String name : names) { 
            Collection<BindInteractor> values = findByName(name);
			values = filterByType(values, "protein");
			values = filterByTaxon(values, taxon);
            for(BindInteractor actor : values) { 
                if(named.containsKey(actor)) { 
                    throw new IllegalArgumentException("2+ names (\"" + name + 
                            "\" and \"" + named.get(actor) + "\") map to the same " +
                            "interactor \"" + actor.toString() + "\"");
                }
                named.put(actor, name);
            }
        }
        
        LinkedList<Pair<String,String>> pairs = new LinkedList<Pair<String,String>>();
        Vector<BindInteractor> ordered = new Vector<BindInteractor>(named.keySet());
        
        for(int i = 0; i < ordered.size(); i++) {
            BindInteractor a = ordered.get(i);
            for(int j = i + 1; j < ordered.size(); j++) {
                BindInteractor b = ordered.get(j);
                if(isInteractingPair(a, b)) { 
                    String na = named.get(a), nb = named.get(b);
                    Pair<String,String> p = new Pair<String,String>(na, nb);
                    pairs.addLast(p);
                }
            }
        }
        
        return pairs;
    }
    
    public Collection<BindInteractor> findByName(String n) { 
        Set<BindInteractor> ints = new HashSet<BindInteractor>();
        for(String k : interactors.keySet()) { 
            if(interactors.get(k).containsName(n)) { 
                ints.add(interactors.get(k));
            }
        }
        return ints;
    }
    
    public Collection<BindInteractor> findByNameFragment(String n) { 
        Set<BindInteractor> ints = new HashSet<BindInteractor>();
        for(String k : interactors.keySet()) { 
            if(interactors.get(k).findName(n)) { 
                ints.add(interactors.get(k));
            }
        }
        return ints;
    }    

    public static class FlatFileRecord { 
        public int rgid, bid;
        public String A_type, B_type;
        public String A_db, B_db;
        public String A_acc, B_acc;
        public int A_id, B_id;
        public int A_tax, B_tax;
        public String AB;
        
        public FlatFileRecord(String line) { 
            String[] array = line.split("\t");
            rgid = Integer.parseInt(array[0]);
            bid = Integer.parseInt(array[1]);
            int i = 2;
            A_type = array[i++];
            A_db = array[i++];
            A_acc = array[i++];
            A_id = Integer.parseInt(array[i++]);
            A_tax = Integer.parseInt(array[i++]);
            B_type = array[i++];
            B_db = array[i++];
            B_acc = array[i++];
            B_id = Integer.parseInt(array[i++]);
            B_tax = Integer.parseInt(array[i++]);
            AB = array[i++];
        }
        
        public BindInteractor createA() { 
            return new BindInteractor(A_id, A_acc, A_db, A_type, A_tax);
        }
        
        public BindInteractor createB() { 
            return new BindInteractor(B_id, B_acc, B_db, B_type, B_tax);
        }
        
        public int hashCode() { 
            int code = 17;
            code += rgid; code *= 37;
            code += bid; code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof FlatFileRecord)) { return false; }
            FlatFileRecord r = (FlatFileRecord)o;
            if(rgid != r.rgid) { return false; }
            if(bid != r.bid) { return false; }
            return true;
        }
        
        public String toString() { 
            return "(" + rgid + "/" + bid + "): " + A_acc + "/" + A_id + " --> " + B_acc + "/" + 
            B_id + ")";
        }
    }

}
