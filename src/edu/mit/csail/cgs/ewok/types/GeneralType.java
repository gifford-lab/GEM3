/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.types;

public class GeneralType implements EchoType {
    
    public GeneralType() {}

    public String getName() {
        return "OBJECT";
    }
    
    public boolean isParameterSubstitute(EchoType paramType) {
        return isSubType(paramType);
    }

    public boolean isReturnSubstitute(EchoType retType) {
        return retType.isSubType(this);
    }

    public boolean isSubType(EchoType et) {
        if(et.equals(this)) { return true; }
        return false;
    }

    public boolean equals(Object o) { 
        return o instanceof GeneralType;
    }
    
    public int hashCode() { 
        return 17;
    }

    public String toString() { return getName(); }
}
