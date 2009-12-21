package edu.mit.csail.cgs.viz.colors;

import java.util.*;
import java.awt.*;

public class ColorRange {
	
	public static int VOLUME=256*3;
	public static int AREA=256*2;
	public static int WIDTH=256;

	private Color[] array;
	
	public ColorRange(int size) { 
		array = new Color[size];
		int blockSize = VOLUME/(size+1);
		int[] c = new int[3];
		
		for(int i = 0; i < size; i++) { 
			int offset = (i+1) * blockSize;
			offsetToColor(offset, c);
			array[i] = new Color(c[0], c[1], c[2]);
		}
	}
	
	public int size() { return array.length; }
	public Color getColor(int i) { return array[i]; }
	
	private void offsetToColor(int offset, int[] array) {
		array[2] = offset % WIDTH;
		array[1] = ((offset-array[2]) / WIDTH) % WIDTH;
		array[0] = ((offset-array[1]) / AREA) % WIDTH;
	}
}
