/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import java.util.*;
import edu.mit.csail.cgs.echo.*;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.utils.Factory;

public class GUIVerbSelector {

	private Class verbClass;
    private SelfDescribingVerbFactory verbfact;
    private String name;
    private String param;
    
    public GUIVerbSelector(GUIVerbSelector s, String p) { 
        verbClass = s.verbClass;
        verbfact = s.verbfact;
        name = s.name;
        param = p;
    }
    
    public GUIVerbSelector(SelfDescribingVerb v) { 
		verbClass = v.getClass();
		verbfact = new SelfDescribingVerbFactory(verbClass);
        String classname = verbClass.getName();
        String[] array = classname.split("\\.");
        name = array.length > 0 ? array[array.length-1] : classname;
        param = "";
    }
    
    public GUIVerbSelector(Class c) { 
		verbClass = c;
		verbfact = new SelfDescribingVerbFactory(verbClass);
        String classname = verbClass.getName();
        String[] array = classname.split("\\.");
        name = array.length > 0 ? array[array.length-1] : classname;
        param = "";
    }
    
    public SelfDescribingVerb getVerb() { 
		return verbfact.createObject(param); 
	}
    
    public String getParam() { return param; }
    public String getClassName() { return verbClass.getName(); }
    public String getName() { return name; }
    public String toString() { return name; }
    
    public int hashCode() { 
        int code = 17;
        code += getClassName().hashCode(); code *= 37;
        code += param.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof GUIVerbSelector)) { return false; }
        GUIVerbSelector s = (GUIVerbSelector)o;
        if(!param.equals(s.param)) { return false; }
        return s.getClassName().equals(getClassName());
    }
}
