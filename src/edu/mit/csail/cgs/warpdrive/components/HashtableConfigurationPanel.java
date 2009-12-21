package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.awt.*;
import javax.swing.*;

public class HashtableConfigurationPanel extends JPanel {
    public static String configuringKey = ".configuring";

    private Hashtable<String,Object> table;
    private java.util.List<String> order;
    private Hashtable<String,JTextField> textfields, doublefields, integerfields;
    private Hashtable<String,JColorChooser> colorfields;
    private Hashtable<String,JCheckBox> checkfields;

    public HashtableConfigurationPanel(Hashtable<String,Object> table) {
        super();
        this.table = table;
        this.order = null;
    }
    public HashtableConfigurationPanel(Hashtable<String,Object> table, java.util.List<String> order) {
        super();
        this.table = table;
        this.order = order;
    }
    public void init() {
        textfields = new Hashtable<String,JTextField>();
        checkfields = new Hashtable<String,JCheckBox>();
        colorfields = new Hashtable<String,JColorChooser>();
        doublefields = new Hashtable<String,JTextField>();
        integerfields = new Hashtable<String,JTextField>();
        
        java.util.List<String> keys;
        if (order == null) {
            keys = new ArrayList<String>();
            keys.addAll(table.keySet());
            Collections.sort(keys);
        } else {
            keys = order;
            for (String s : table.keySet()) {
                if (!keys.contains(s)) {
                    keys.add(s);
                }
            }
        }
        GridBagLayout gridbag = new GridBagLayout();
        setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();        
        constraints.weightx = 1.0;
        for (String k : keys) {
            if (k == configuringKey) {continue;}
            JLabel label = new JLabel(k);            
            constraints.fill = GridBagConstraints.BOTH;
            label.setHorizontalAlignment(SwingConstants.RIGHT);
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            gridbag.setConstraints(label,constraints);
            JComponent comp = null;
            add(label);
            if (table.get(k) instanceof String) {
                JTextField field = new JTextField(table.get(k).toString());
                textfields.put(k,field);
                comp = field;
            } else if (table.get(k) instanceof Boolean) {
                JCheckBox field = new JCheckBox("",(Boolean)table.get(k));
                checkfields.put(k,field);
                comp = field;
            } else if (table.get(k) instanceof Color) {
                Color c = (Color)table.get(k);
                //                ColorSelection selection = 
                //                    new ColorSelection(c);
                JColorChooser selection = new JColorChooser(c);
                comp = selection;
                colorfields.put(k,selection);
            } else if (table.get(k) instanceof Double) {
                JTextField field = new JTextField(table.get(k).toString());
                comp = field;
                doublefields.put(k,field);
            } else if (table.get(k) instanceof Integer) { 
                JTextField field = new JTextField(table.get(k).toString());
                comp = field;
                integerfields.put(k, field);
            }
            if (comp != null) {
                constraints.gridwidth = GridBagConstraints.REMAINDER;
                gridbag.setConstraints(comp,constraints);
                add(comp);
            }
        }
    }

    public void parse() {
        for (String k : textfields.keySet()) {
            table.put(k,textfields.get(k).getText());
        }
        for (String k : doublefields.keySet()) {
            table.put(k,Double.parseDouble(doublefields.get(k).getText()));
        }
        for (String k : integerfields.keySet()) {
            table.put(k,Integer.parseInt(integerfields.get(k).getText()));
        }
        for (String k : checkfields.keySet()) {
            table.put(k,checkfields.get(k).isSelected());
        }
        for(String k : colorfields.keySet()) { 
            Color c = colorfields.get(k).getColor();
            if(c != null) { table.put(k,c); }
        }
        
        
        clearConfiguring();
    }

    public void clearConfiguring() {
        table.remove(configuringKey);
    }
}
