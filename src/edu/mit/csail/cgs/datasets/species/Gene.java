package edu.mit.csail.cgs.datasets.species;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;

/** A <code>Gene</code> is a <code>NamedStrandedRegion</code> that stores
 * an id, a set of aliases, and a source.  
 *
 * The distinction between name and id can be fuzzy, but name is usually
 * a biological name (eg, CDC15) and id is usually a systematic name (eg YAR019C).  
 * The name and the id may be the same.
 */
public class Gene extends NamedStrandedRegion {
    
    private int DBID;
    private String id;
    private Set<String> aliases;
    private String source;
    private String rep;
    
    public Gene(Genome g, String c, int start, int end, String name, String id, String src) {
        super(g,c,start,end,name,' ');
        this.id = id;
        aliases = new HashSet<String>();
        source = src;
        DBID = -1;
        rep = null;
    }

    public Gene(Genome g, String c, int start, int end, String name, String id, char str, String src) {
        super(g,c,start,end,name,str);
        this.id = id;
        aliases = new HashSet<String>();
        source = src;
        DBID = -1;
        rep = null;
    }

    public Gene(Genome g, String c, int start, int end, String name, String id, char str, String src, int dbid) {
        super(g,c,start,end,name,str);
        this.id = id;
        aliases = new HashSet<String>();
        source = src;
        DBID = dbid;
        rep = null;
    }

    public int getDBID() { return DBID; }
    public String getID () {return id;}
    
    public String toString() {
        if(rep == null) { rep = getStringRep(); }
        return rep;
    }
    
    public String getStringRep() { 
        StringBuilder sb = new StringBuilder();
        sb.append(getID());
        Collection<String> otherNames = getNonIDNames();
        if(otherNames.size() > 0) {
            sb.append(" (");
            boolean first = true;
            for(String n : otherNames) {
                sb.append((first ? "" : ", ") + n);
                first = false;
            }
            sb.append(")");
        }
        return sb.toString();
    }
    
    public String getSource() { return source; }
    /** Returns the transcription start site.  This is the result of <code>getStart()</code>
     *  if the gene is on the + strand and the result of <code>getEnd()</code> if the
     *  gene is on the - strand.
     */
    public int getTSS() { return getFivePrime();}
    
    public Collection<String> getAliases() { return aliases; }
    /** Returns a collections of Strings that includes the name, id, and aliases. */
    public Collection<String> getAllNames() { 
        TreeSet<String> names = new TreeSet<String>();
        names.addAll(aliases);
        names.add(id);
        names.add(getName());
        return names;
    }
    
    public Collection<String> getNonIDNames() { 
        TreeSet<String> names = new TreeSet<String>();
        names.addAll(aliases);
        names.add(getName());
        if(names.contains(getID())) { names.remove(getID()); }
        return names;
    }
    
    public void addAlias(String a) {
        if(a != null && !a.equals(getName()) && !a.equals(id)) { 
            aliases.add(a); 
        }
    }
    
    public boolean hasName(String a) { 
        return getName().equals(a) || id.equals(a) || aliases.contains(a);
    }
    
    public boolean sharesName(Gene g) { 
        if(g.hasName(getName()) || g.hasName(id)) { return true; }
        for(String a : aliases) { 
            if(g.hasName(a)) { return true; }
        }
        return false;
    }
    
    public int hashCode() { 
        int code = super.hashCode();
        code += id.hashCode(); code *= 37;
        code += getStrand(); code *= 37;
        code += source.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof Gene)) { return false; }
        Gene g = (Gene)o;
        if(!super.equals(g)) { return false; }
        if(!id.equals(g.id)) { return false; }
        if(getStrand() != g.getStrand()) { return false; }
        if(!source.equals(g.source)) { return false; }
        return true;
    }
}
