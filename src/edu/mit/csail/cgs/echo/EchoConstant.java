package edu.mit.csail.cgs.echo;

import java.util.*;

import edu.mit.csail.cgs.echo.gui.EchoGUIConstant;
import edu.mit.csail.cgs.ewok.types.*;

/**
 * @author Timothy Danford
 */
public class EchoConstant { 

	private SelfDescribingConstant constant;
	private Object value;
    
    private EchoGUIConstant guiPeer; 

	public EchoConstant(SelfDescribingConstant sdc) { 
		constant = sdc;
		value = null;
        guiPeer = null;
	}
    
    public void setPeer(EchoGUIConstant p) { guiPeer = p; }
    public EchoGUIConstant getPeer() { return guiPeer; }
	
	public EchoType getConstantClass() { 
		return constant.getConstantClass(); 
	}

	public void init() { 
		value = null;
	}

	public Object evaluate() { 
		if(value == null) { value = constant.getConstantValue(); }
		return value;
	}

	public SelfDescribingConstant getConstant() { return constant; }
		
	public String[] findMatchingParamNames(Reverb obj) {
		SelfDescribingVerb verb = obj.getVerb();
		Vector<String> p = new Vector<String>();
		
		String[] params = verb.getParameterNames();
		EchoType[] paramClasses = verb.getParameterClasses();
		for(int i = 0; params != null && i < params.length; i++) {
			if(constant.getConstantClass().isSubType(paramClasses[i])) { 
				p.add(params[i]);
			}
		}
		
		if(p.size() > 0) { 
			return p.toArray(new String[p.size()]);
		} else { 
			return null;
		}
	}	
}
