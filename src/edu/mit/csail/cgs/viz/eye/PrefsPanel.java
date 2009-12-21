package edu.mit.csail.cgs.viz.eye;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.models.*;


public class PrefsPanel<M extends Model> extends JPanel {
		
    private M value;
    private ModelFieldAnalysis analysis;
    private Map<Field,JTextField> stringFields, integerFields, doubleFields;
    private Map<Field,JRadioButton[]> booleanFields;
    private Map<Field,JColorChooser> colorFields;
		
    public PrefsPanel(M mdl) { 
        value = mdl;
        analysis = new ModelFieldAnalysis(mdl.getClass());
        stringFields = new HashMap<Field,JTextField>();
        integerFields = new HashMap<Field,JTextField>();
        doubleFields = new HashMap<Field,JTextField>();
        booleanFields = new HashMap<Field,JRadioButton[]>();
        colorFields = new HashMap<Field,JColorChooser>();
        init();
    }
		
    public int numFields() { 
        return stringFields.size() + booleanFields.size() + 
				doubleFields.size() + integerFields.size() + colorFields.size();
    }
    public Dimension getPreferredSize() {
        return new Dimension(600, numFields() * 30 + 250 * colorFields.size() + 50);
    }
	
    public void saveToModel() { 
        for(Field f : stringFields.keySet()) { 
            JTextField cmp = stringFields.get(f);
            String fval = cmp.getText();
            try {
                f.set(value, fval);
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
        }
        for(Field f : integerFields.keySet()) { 
            JTextField cmp = integerFields.get(f);
            try {
                Integer fval = Integer.parseInt(cmp.getText());
                f.set(value, fval);
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            } catch(NumberFormatException nfe) { 
                System.err.println(String.format("%s is an invalid Integer", cmp.getText()));
            }
        }
        for(Field f : doubleFields.keySet()) { 
            JTextField cmp = doubleFields.get(f);
            try {
                Double fval = Double.parseDouble(cmp.getText());
                f.set(value, fval);
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            } catch(NumberFormatException nfe) { 
                System.err.println(String.format("%s is an invalid Double", cmp.getText()));
            }
        }
        for(Field f : booleanFields.keySet()) { 
            JRadioButton[] array = booleanFields.get(f);
            try {
                Boolean fval = array[0].isSelected();
                f.set(value, fval);
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
        }
        for (Field f : colorFields.keySet()) {
            JColorChooser cmp = colorFields.get(f);
            try {
                f.set(value, cmp.getColor());
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
        }
    }
		
    private void init() {
        Vector<Field> sfs = analysis.findTypedFields(String.class);
        Vector<Field> ifs = analysis.findTypedFields(Integer.class);
        Vector<Field> dfs = analysis.findTypedFields(Double.class);
        Vector<Field> bfs = analysis.findTypedFields(Boolean.class);
        Vector<Field> cfs = analysis.findTypedFields(Color.class);
			
        int fieldCount = sfs.size() + ifs.size() + dfs.size() + bfs.size();
        GridBagLayout gridbag = new GridBagLayout();
        setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();        
        constraints.weightx = 1.0;
        constraints.fill = GridBagConstraints.BOTH;

        for(Field f : sfs) { 
            String name = f.getName();
            String fvalue = (String)analysis.get(f.getName(), value);
            JLabel label = new JLabel(String.format("%s:", name));
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            JTextField component = new JTextField();
            if(fvalue != null) { component.setText(fvalue); }
            stringFields.put(f, component);
            gridbag.setConstraints(label,constraints);        
            add(label);        
            constraints.gridwidth = GridBagConstraints.REMAINDER;
            gridbag.setConstraints(component, constraints);
            add(component);
        }

        for(Field f : ifs) { 
            String name = f.getName();
            Integer fvalue = (Integer)analysis.get(f.getName(), value);
            JTextField component = new JTextField();
            if(fvalue != null) { component.setText(String.valueOf(fvalue)); }
            integerFields.put(f, component);
            JLabel label = new JLabel(String.format("%s:", name));
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            gridbag.setConstraints(label,constraints);        
            add(label);        
            constraints.gridwidth = GridBagConstraints.REMAINDER;
            gridbag.setConstraints(component, constraints);
            add(component);
        }
        for(Field f : dfs) { 
            String name = f.getName();
            Double fvalue = (Double)analysis.get(f.getName(), value);
            JTextField component = new JTextField();
            if(fvalue != null) { component.setText(String.valueOf(fvalue)); }
            doubleFields.put(f, component);
            JLabel label = new JLabel(String.format("%s:", name));
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            gridbag.setConstraints(label,constraints);        
            add(label);        
            constraints.gridwidth = GridBagConstraints.REMAINDER;
            gridbag.setConstraints(component, constraints);
            add(component);
        }
			
        for(Field f : bfs) { 
            String name = f.getName();
            Boolean fvalue = (Boolean)analysis.get(f.getName(), value);
            ButtonGroup group = new ButtonGroup();
            JRadioButton[] array = new JRadioButton[] { 
                new JRadioButton("True"),
                new JRadioButton("False")
            };
            group.add(array[0]);
            group.add(array[1]);
            if(fvalue == null || fvalue) { 
                array[0].setSelected(true);
            } else { 
                array[1].setSelected(true);
            }
				
            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout());
            buttonPanel.add(array[0]);
            buttonPanel.add(array[1]);
				
            booleanFields.put(f, array);
				
            JLabel label = new JLabel(String.format("%s:", name));
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            gridbag.setConstraints(label,constraints);        
            add(label);        
            constraints.gridwidth = GridBagConstraints.REMAINDER;
            gridbag.setConstraints(buttonPanel, constraints);
            add(buttonPanel);
        }
        for (Field f : cfs) {
            String name = f.getName();
            Color fvalue = (Color)analysis.get(f.getName(),value);
            JColorChooser chooser = new JColorChooser();
            if (fvalue != null) {
                chooser.setColor(fvalue);
            }
            colorFields.put(f, chooser);
            JLabel label = new JLabel(String.format("%s:",name));
            constraints.gridwidth = GridBagConstraints.RELATIVE;
            gridbag.setConstraints(label,constraints);        
            add(label);        
            constraints.gridwidth = GridBagConstraints.REMAINDER;
            gridbag.setConstraints(chooser, constraints);
            add(chooser);
        }

    }		
    public M getModelValue() { return value; }
}