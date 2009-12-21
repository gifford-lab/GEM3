package edu.mit.csail.cgs.viz.colors;

import java.awt.Color;
import java.util.*;

/* contains a list of colors and returns them when asked in
   such a way that it repeats a color as infrequently as possible */

public class ColorSet {

    private Map<String,Color> assigned;
    
    public static Color translateColorName(String n) {
        n = n.toLowerCase();
        if(n.equals("blue")) { return Color.blue; }
        if(n.equals("green")) { return Color.green; }
        if(n.equals("red") ) { return Color.red; }
        if(n.equals("cyan")) { return Color.cyan; }
        if(n.equals("orange")) { return Color.orange; }
        if(n.equals("magenta")) { return Color.magenta; }
        if(n.equals("lightgray")) { return Color.LIGHT_GRAY; }
        if(n.equals("darkgray")) { return Color.DARK_GRAY; }
        if(n.equals("pink")) { return Color.pink; }
        if(n.equals("yellow")) { return Color.yellow; }
        return null;
    }
    
    private static Color[] defaultColors = {Color.BLUE,
                                            Color.GREEN,
                                            Color.RED,
                                            Color.CYAN,
                                            Color.ORANGE,
                                            Color.MAGENTA,
                                            Color.PINK.darker(),
                                            new Color(180,0,230),
                                            Color.DARK_GRAY,
                                            Color.YELLOW};
    
    private Vector<Color> colors;
    private int pos;
    private Random rand;

    public ColorSet() {
        pos = 0;
        rand = new Random();
        colors = new Vector<Color>();
        assigned = new HashMap<String,Color>();
        for(int i = 0; i < defaultColors.length; i++) { colors.add(defaultColors[i]); }
    }
    
    public Color sampleColor() { 
    	int r = rand.nextInt(256);
    	int g = rand.nextInt(256);
    	int b = rand.nextInt(256);
    	return new Color(r, g, b);
    }
    
    public void addColor() { colors.add(sampleColor()); }
    
    public void ensureColorCount(int c) { 
    	while(colors.size() < c) { addColor(); }
    }

    public Color getColor() {
        Color c = colors.get(pos++);
        pos = pos % colors.size();
        return c;
    }
    public Color getColor(String key) {
        if (assigned.containsKey(key)) {
            return assigned.get(key);
        }
        Color c = getColor();
        assigned.put(key,c);
        return c;
    }

    public Color[] getColors() {
    	return (Color[])colors.toArray(new Color[colors.size()]);
    }

    public int colorCount() { return colors.size(); }
    public Color colorAt(int i) { return colors.get(i); }

    public void reset() {
        pos = 0;
        assigned.clear();
    }

}
