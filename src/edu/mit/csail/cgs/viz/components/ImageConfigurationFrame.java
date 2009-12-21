package edu.mit.csail.cgs.viz.components;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;


public class ImageConfigurationFrame extends JFrame implements ActionListener {
    private AbstractPaintable parent;
    private JCheckBox rasterbox;
    private JTextField widthfield, heightfield;
    private JButton okbutton, cancelbutton;

    public ImageConfigurationFrame(AbstractPaintable p) {
        parent = p;
        JLabel boxlabel = new JLabel("Configure Save-as-image");
	    rasterbox = new JCheckBox("Raster Image?",parent.getImageRaster());
	    JLabel widthlabel = new JLabel("Width");
	    JLabel heightlabel = new JLabel("Height");
	    widthfield = new JTextField(Integer.toString(parent.sImageWidth));
	    heightfield = new JTextField(Integer.toString(parent.sImageHeight));
	    okbutton = new JButton("OK");
	    cancelbutton = new JButton("Cancel");
	    okbutton.addActionListener(this);
	    cancelbutton.addActionListener(this);
	    
	    JPanel toppanel = new JPanel();
	    toppanel.setLayout(new BorderLayout());
	    
	    JPanel buttonpanel = new JPanel();
	    buttonpanel.add(okbutton);
	    buttonpanel.add(cancelbutton);
	    toppanel.add(buttonpanel,BorderLayout.SOUTH);
	
	    JPanel infopanel = new JPanel();
	    infopanel.setLayout(new BorderLayout());
	    infopanel.add(rasterbox,BorderLayout.NORTH);
	    JPanel textpanel = new JPanel();
	    textpanel.setLayout(new GridLayout(2,2));
	    textpanel.add(widthlabel);
	    textpanel.add(widthfield);
	    textpanel.add(heightlabel);
	    textpanel.add(heightfield);
	    infopanel.add(textpanel,BorderLayout.CENTER);
	    
	    toppanel.add(infopanel,BorderLayout.CENTER);
	
	    getContentPane().add(toppanel);
	    setMinimumSize(new Dimension(150,150));
	    setSize(getPreferredSize());
	    pack();
	    setVisible(true);
	}
	
	public void actionPerformed (ActionEvent e) {
	    if (e.getSource() == okbutton) {
	        parent.setImageRaster(rasterbox.isSelected());
	        try {
	            parent.sImageWidth = Integer.parseInt(widthfield.getText());
	        } catch (NumberFormatException ex) {
	        }
	        try {
	            parent.sImageHeight = Integer.parseInt(heightfield.getText());
	        } catch (NumberFormatException ex) {
	        }
	
	
	        this.dispose();
	    } else if (e.getSource() == cancelbutton) {
	        this.dispose();
	    }
	}     
}	

