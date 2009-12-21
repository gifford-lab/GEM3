/*
 * Created on Nov 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;

public class ProvisionalGOTerm implements Comparable<ProvisionalGOTerm> { 
    
    private OboObject base;
    private String id, name;
    private HashSet<ProvisionalGOTerm> parents;
    private HashSet<ProvisionalGOTerm> children;
    private int insertPriority;
    
    public ProvisionalGOTerm(OboObject b) { 
        base = b;
        insertPriority = 0;
        Collection<String> nvs = b.getValues("name");
        Collection<String> ivs = b.getValues("id");
        if(nvs.size() != 1) { base.debugPrint(System.err); throw new IllegalArgumentException(); }
        if(ivs.size() != 1) { base.debugPrint(System.err); throw new IllegalArgumentException(); }
        Iterator<String> itr = nvs.iterator();
        name = itr.next();
        itr = ivs.iterator();
        id = itr.next();
        parents = new HashSet<ProvisionalGOTerm>();
        children = new HashSet<ProvisionalGOTerm>();
    }
    
    public String getName() { return name; }
    public String getID() { return id; }
    public OboObject getObject() { return base; }
    public int getInsertPriority() { return insertPriority; }
    
    public Collection<ProvisionalGOTerm> getParents() { return parents; }
    
    public boolean isChild(ProvisionalGOTerm parent) { 
        return parents.contains(parent);
    }
    
    public boolean isDescendant(ProvisionalGOTerm p) {
        if(equals(p)) { return false; }
        if(isChild(p)) { return true; }
        for(ProvisionalGOTerm par : parents) { 
            if(par.isDescendant(p)) { return true; }
        }
        return false;
    }
    
    public int compareTo(ProvisionalGOTerm p) { 
        if(insertPriority > p.insertPriority) { return -1; }
        if(insertPriority < p.insertPriority) { return 1; }
        return id.compareTo(p.id);
    }
    
    private void updateInsertPriority(int p) { 
        updateInsertPriority(new LinkedList<String>(), p);
    }
    
    private void updateInsertPriority(LinkedList<String> path, int p) {
        
        if(path.contains(id)) {
            /*
            System.out.print("Cycle! : ");
            for(String pid : path) { System.out.print(" " + pid); }
            System.out.println(" --> " + id);
            */
        } else { 
            if(p > insertPriority) { 
                insertPriority = p;
                path.addLast(id);
                for(ProvisionalGOTerm par : parents) { 
                    par.updateInsertPriority(path, p + 1);
                }
                path.removeLast();
            }
        }
    }
    
    public void addParent(ProvisionalGOTerm p) { 
        parents.add(p);
        p.children.add(this);
        p.updateInsertPriority(insertPriority + 1);
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append(id);
        sb.append(" [");
        for(ProvisionalGOTerm gt : parents) { 
            sb.append(" " + gt.id);
        }
        sb.append(" ]");
        return sb.toString();
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ProvisionalGOTerm)) { return false; }
        ProvisionalGOTerm g = (ProvisionalGOTerm)o;
        if(!name.equals(g.name)) { return false; }
        if(!id.equals(g.id)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += name.hashCode(); code *= 37;
        code += id.hashCode(); code *= 37;
        return code;
    }
}

