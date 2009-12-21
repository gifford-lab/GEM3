package edu.mit.csail.cgs.ewok.types;

public class ClassType implements EchoType {
	
	private Class internalClass;
	private String name;
	
	public ClassType(Class c) { 
		internalClass = c;
		String cn = c.getName();
		String[] array = cn.split("\\.");
		name = array[array.length-1];
	}
	
	public Class getInternalClass() { return internalClass; }

	public String getName() {
		return name;
	}

	public boolean isSubType(EchoType et) {
        if(et instanceof GeneralType) { return true; }
		if(!(et instanceof ClassType)) { return false; }
		Class oc = ((ClassType)et).getInternalClass();
		return isSubclass(internalClass, oc);
	}
	
    public static boolean isSubclass(Class subc, Class superc) {
        return superc.isAssignableFrom(subc);
    }

	public boolean isParameterSubstitute(EchoType paramType) {
		return isSubType(paramType);
	}

	public boolean isReturnSubstitute(EchoType retType) {
		return retType.isSubType(this);
	}
	
	public int hashCode() { 
		int code = 17;
		code += internalClass.hashCode(); code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ClassType)) { return false; }
		ClassType ct = (ClassType)o;
		if(!ct.internalClass.equals(internalClass)) { return false; }
		return true;
	}
}
