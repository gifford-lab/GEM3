package edu.mit.csail.cgs.ewok.types;

public interface EchoType {
	
	public static final EchoType OBJECT_CLASS = new GeneralType();
	
	public String getName();
	public boolean isSubType(EchoType et);
	public boolean isParameterSubstitute(EchoType paramType);
	public boolean isReturnSubstitute(EchoType retType);
}
