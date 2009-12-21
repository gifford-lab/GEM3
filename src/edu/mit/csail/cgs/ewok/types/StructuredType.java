package edu.mit.csail.cgs.ewok.types;

import java.util.*;

public class StructuredType implements EchoType {

	private TreeMap<String,EchoType> subtypes;
	private String name;
	private String totalName;
    
    public StructuredType(String n) { 
        subtypes = new TreeMap<String,EchoType>();
        name = n;
        buildTotalName();
    }
	
	public StructuredType(String n, Collection<EchoType> ts) { 
		name = n;
		subtypes = new TreeMap<String,EchoType>();
        
        int i = 0;
        for(EchoType t : ts) {
            subtypes.put(String.valueOf(i), t);
            i++;
        }
        
		buildTotalName();
	}
	
    public StructuredType(String n, String tn, EchoType t) { 
        name = n;
        subtypes = new TreeMap<String,EchoType>();
        subtypes.put(tn, t);
        buildTotalName();
    }
    
    public StructuredType(String n, String tn1, EchoType t1, String tn2, EchoType t2) { 
        name = n;
        subtypes = new TreeMap<String,EchoType>();
        subtypes.put(tn1, t1);
        subtypes.put(tn2, t2);
        buildTotalName();
    }
    
    public StructuredType(StructuredType base, String tn, EchoType t) { 
        name = base.name;
        subtypes = new TreeMap<String,EchoType>(base.subtypes);
        if(subtypes.containsKey(tn)) { throw new IllegalArgumentException(tn); }
        subtypes.put(tn, t);
        buildTotalName();
    }
    
	private void buildTotalName() {
        StringBuilder tnb = new StringBuilder();
		tnb.append(name + "( ");
        TreeSet<String> keys = new TreeSet<String>(subtypes.keySet());
        for(String key : keys) { 
            tnb.append(key + ":" + subtypes.get(key).getName() + " ");
        }
		tnb.append(")");
        totalName = tnb.toString();
	}
	
	public String toString() { return getName(); }
	public String getName() { return totalName; }
    public Set<String> getTypeKeys() { return subtypes.keySet(); }
    public EchoType getTypeValue(String k) { return subtypes.get(k); }

	public boolean isSubType(EchoType et) {
        if(et instanceof GeneralType) { return true; }
		if(!(et instanceof StructuredType)) { return false; }
        
		StructuredType st = (StructuredType)et;
		if(!name.equals(st.name)) { return false; }
        
        for(String key : st.subtypes.keySet()) { 
            if(!subtypes.containsKey(key)) { return false; } 
            if(!subtypes.get(key).isSubType(st.subtypes.get(key))) { return false; }
        }        
		return true;
	}

	public boolean isParameterSubstitute(EchoType paramType) {
		return isSubType(paramType);
	}

	public boolean isReturnSubstitute(EchoType retType) {
		return retType.isSubType(this);
	}
	
	public int hashCode() { 
		int code = 17;
		code += name.hashCode(); code *= 37;
        for(String key : subtypes.keySet()) { code += key.hashCode(); }
        code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof StructuredType)) { return false; }
		StructuredType st = (StructuredType)o;
		if(!name.equals(st.name)) { return false; }
		if(st.subtypes.size() != subtypes.size()) { return false; }
        
        for(String key : subtypes.keySet()) { 
            if(!st.subtypes.containsKey(key)) { return false; }
            if(!subtypes.get(key).equals(st.subtypes.get(key))) { return false; }
        }
		return true;
	}
}
