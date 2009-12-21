package edu.mit.csail.cgs.metagenes.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

public class OptionsFrame extends JFrame{
	
	private JTextField fontSizeField = new JTextField(10);
	private JTextField maxColorField = new JTextField(10);
	private JTextField minColorField = new JTextField(10);
	private JTextField lineWeightField = new JTextField(10);
	private JLabel profileMax, profileMin;
	private JButton fontSizeEntry;
	private JButton maxColorEntry, minColorEntry;
	private JButton lineWeightEntry;
	private JButton colorChangeEntry;
	private JButton finished;
	private final ProfilePanel panel;
	private final ProfileLinePanel linepanel;
	private JColorChooser colorPick;

	public OptionsFrame(final ProfilePanel p, final ProfileLinePanel lp){
		panel = p;
		linepanel=lp;

		//Initialize buttons
		fontSizeEntry = new JButton(createFontUpdateAction());
		maxColorEntry = new JButton(createMaxColorUpdateAction());
		minColorEntry = new JButton(createMinColorUpdateAction());
		lineWeightEntry = new JButton(createLineWeightUpdateAction());
		colorChangeEntry = new JButton(createColorUpdateAction());
		finished= new JButton(createFinishAction());
		
		//Font size panel
		JPanel fp = new JPanel();
		fp.setLayout(new BoxLayout(fp, BoxLayout.LINE_AXIS));
		JLabel message = new JLabel("Font Size:");
		fp.add(message);
		fontSizeEntry.setMaximumSize(new Dimension(100, 30));
		fp.add(fontSizeField);
		fp.add(fontSizeEntry);
		fp.setSize(new Dimension(300,30));
		fp.setMaximumSize(new Dimension(300,30));
		fp.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
		
		//Color chooser
		JPanel cp = new JPanel();
		cp.setLayout(new BoxLayout(cp, BoxLayout.PAGE_AXIS));
		colorPick = new JColorChooser(panel.getPeakColor());
		cp.add(colorPick);
		cp.add(colorChangeEntry);
		
		//Meta-peak Style
		JRadioButton rad_histo = new JRadioButton("Histogram");
		JRadioButton rad_line = new JRadioButton("Line");
		rad_histo.addActionListener(createRadioUpdateAction());
		rad_histo.setActionCommand("Histo");
		rad_line.addActionListener(createRadioUpdateAction());
		rad_line.setActionCommand("Line");
		rad_line.setSelected(true);
		ButtonGroup styleMenu = new ButtonGroup();
		styleMenu.add(rad_histo);
		styleMenu.add(rad_line);
		JPanel sp = new JPanel();
		sp.setLayout(new BoxLayout(sp, BoxLayout.PAGE_AXIS));
		JLabel sLabel = new JLabel("Meta-Peak Style:");
		sp.add(sLabel);
		sp.add(rad_histo);
		sp.add(rad_line);
		
		//Line profile max/min color panel
		JPanel lmp = new JPanel();
		lmp.setLayout(new BoxLayout(lmp, BoxLayout.LINE_AXIS));
		profileMax = new JLabel(String.format("Value for line max color:", linepanel.getMaxColorVal()));
		lmp.add(profileMax);
		maxColorEntry.setMaximumSize(new Dimension(200, 30));
		lmp.add(maxColorField);
		lmp.add(maxColorEntry);
		lmp.setSize(new Dimension(400,30));
		lmp.setMaximumSize(new Dimension(400,30));
		lmp.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
		JPanel lmp2 = new JPanel();
		lmp2.setLayout(new BoxLayout(lmp2, BoxLayout.LINE_AXIS));
		profileMin = new JLabel(String.format("Value for line min color:", linepanel.getMinColorVal()));
		lmp2.add(profileMin);
		minColorEntry.setMaximumSize(new Dimension(200, 30));
		lmp2.add(minColorField);
		lmp2.add(minColorEntry);
		lmp2.setSize(new Dimension(400,30));
		lmp2.setMaximumSize(new Dimension(400,30));
		lmp2.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
		
		//Line weight panel
		JPanel lwp = new JPanel();
		lwp.setLayout(new BoxLayout(lwp, BoxLayout.LINE_AXIS));
		JLabel lineWeightMessage = new JLabel("Line Profile Weight");
		lwp.add(lineWeightMessage);
		lineWeightEntry.setMaximumSize(new Dimension(100, 30));
		lwp.add(lineWeightField);
		lwp.add(lineWeightEntry);
		lwp.setSize(new Dimension(300,30));
		lwp.setMaximumSize(new Dimension(300,30));
		lwp.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
		
		
		//Container for everything
		Container c = (Container)getContentPane();
		c.setLayout(new BoxLayout(c, BoxLayout.PAGE_AXIS));
		c.add(fp);
		c.add(Box.createRigidArea(new Dimension(0, 10)));
		c.add(cp);
		c.add(Box.createRigidArea(new Dimension(0, 10)));
		c.add(sp);
		c.add(Box.createRigidArea(new Dimension(0, 10)));
		c.add(lmp);
		c.add(lmp2);
		c.add(Box.createRigidArea(new Dimension(0, 10)));
		c.add(lwp);
		c.add(Box.createRigidArea(new Dimension(0, 10)));
		c.add(finished);
		//c.setVisible(true);
		
		//setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	}
	
	private Action createFontUpdateAction(){
		return new AbstractAction("OK") { 
			public void actionPerformed(ActionEvent e) { 
				String tmp = fontSizeField.getText();
				int fSize;
				if(tmp.length()>0){
					fSize = Integer.valueOf(tmp);
				}else{
					fSize=12;
				}
				panel.updateFontSize(fSize);
				linepanel.updateFontSize(fSize);
				OptionsFrame.this.setVisible(false);
			} 
		};
	}
	private Action createLineWeightUpdateAction(){
		return new AbstractAction("OK") { 
			public void actionPerformed(ActionEvent e) { 
				String tmp = lineWeightField.getText();
				int w;
				if(tmp.length()>0){
					w = Integer.valueOf(tmp);
				}else{
					w=1;
				}
				linepanel.updateLineWeight(w);
				OptionsFrame.this.setVisible(false);
			} 
		};
	}
	private Action createColorUpdateAction(){
		return new AbstractAction("Update Color") { 
			public void actionPerformed(ActionEvent e) { 
				Color c = colorPick.getColor();
				if(c==null){
					c=Color.blue; 
				}
				panel.updateColor(c);
				linepanel.updateColor(c);
				OptionsFrame.this.setVisible(false);
			} 
		};
	}
	private Action createRadioUpdateAction(){
		return new AbstractAction() { 
			public void actionPerformed(ActionEvent e) { 
				panel.setStyle(e.getActionCommand());
			} 
		};
	}
	private Action createMaxColorUpdateAction(){
		return new AbstractAction("OK") { 
			public void actionPerformed(ActionEvent e) { 
				String tmp = maxColorField.getText();
				double max;
				if(tmp.length()>0){
					max = Double.valueOf(tmp);
				}else{
					max=linepanel.getMaxColorVal();
				}
				linepanel.setMaxColorVal(max);
				OptionsFrame.this.setVisible(false);
			} 
		};
	}
	private Action createMinColorUpdateAction(){
		return new AbstractAction("OK") { 
			public void actionPerformed(ActionEvent e) { 
				String tmp = minColorField.getText();
				double min;
				if(tmp.length()>0){
					min = Double.valueOf(tmp);
				}else{
					min=linepanel.getMinColorVal();
				}
				linepanel.setMinColorVal(min);
				OptionsFrame.this.setVisible(false);
			} 
		};
	}
	private Action createFinishAction(){
		return new AbstractAction("Done") {
			public void actionPerformed(ActionEvent e) {
				Color c = colorPick.getColor();
				if(c!=null){
					panel.updateColor(c);
					linepanel.updateColor(c);
				}
				String tmp = fontSizeField.getText();
				int fSize;
				if(tmp.length()>0){
					fSize = Integer.valueOf(tmp);
					panel.updateFontSize(fSize);
					linepanel.updateFontSize(fSize);
				}
				tmp = maxColorField.getText();
				double max=linepanel.getMaxColorVal();
				if(tmp.length()>0){
					max = Double.valueOf(tmp);
					linepanel.setMaxColorVal(max);
				}
				tmp = minColorField.getText();
				double min=linepanel.getMinColorVal();
				if(tmp.length()>0){
					min = Double.valueOf(tmp);
					linepanel.setMinColorVal(min);
				}
				tmp = lineWeightField.getText();
				int w=1;
				if(tmp.length()>0){
					w = Integer.valueOf(tmp);
					linepanel.updateLineWeight(w);
				}
				OptionsFrame.this.setVisible(false);
			}
		};
	}
	public void startup() {
		// We want to call setVisible() and *then* pack() -- but once setVisible() has 
		// been called on a Swing component, you shouldn't call any other methods of that 
		// component except *from the Swing thread*.  Therefore, this hack.  
		Runnable r = new Runnable() {
			public void run() { 
				//OptionsFrame.this.setPreferredSize(new Dimension(400,500));
				OptionsFrame.this.setLocation(getX() + 150, getY() + 50);
				//OptionsFrame.this.setVisible(true);
				OptionsFrame.this.pack();
			}
		};
		SwingUtilities.invokeLater(r);
	}
}
