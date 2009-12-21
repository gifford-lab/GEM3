/*
 * Author: tdanford
 * Date: Jun 19, 2008
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;

public class TextLayout {
	
	public TextLayout() { 
	}
	
	public Vector<String> paragraphLayout(String text, int lineLength) { 
		String[] wa = text.split("\\s+");
		Vector<String> words = new Vector<String>();
		for(int i = 0; i < wa.length; i++) { 
			words.add(wa[i]);
		}
		Vector<Vector<String>> laidOutWords = paragraphVectorLayout(words, lineLength);
		Vector<String> layout = new Vector<String>();
		for(int i = 0; i < laidOutWords.size(); i++) { 
			StringBuilder sb = new StringBuilder();
			Vector<String> ws = laidOutWords.get(i);
			for(int j = 0; j < ws.size(); j++) { 
				if(sb.length() > 0) { sb.append(" "); }
				sb.append(ws.get(j));
			}
			layout.add(sb.toString());
		}
		return layout;
	}

	public Vector<Vector<String>> paragraphVectorLayout(Vector<String> words, int lineLength) { 
		Vector<Vector<String>> lines = new Vector<Vector<String>>();
		
		int lineOffset = 0;
		Vector<String> currentLine = new Vector<String>();
		lines.add(currentLine);
		for(int wi = 0; wi < words.size(); wi++) { 
			String word = words.get(wi);
			if(lineOffset + word.length() > lineLength) { 
				lineOffset = 0;
				currentLine = new Vector<String>();
				lines.add(currentLine);
			}
			
			currentLine.add(word);
			lineOffset += word.length() + 1;
		}
		
		return lines;
	}
}
