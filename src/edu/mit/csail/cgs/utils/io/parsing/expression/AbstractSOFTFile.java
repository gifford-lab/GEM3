package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;

public class AbstractSOFTFile {

	protected SOFTAttributes attrs;
	
	public AbstractSOFTFile() { 
		attrs = new SOFTAttributes();
	}
	
	public SOFTAttributes getAttributes() { return attrs; }
	
	public void addAttribute(String key, String value) { 
		attrs.addKeyValue(key, value);
	}
}
