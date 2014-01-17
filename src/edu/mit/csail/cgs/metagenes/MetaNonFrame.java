package edu.mit.csail.cgs.metagenes;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.FlowLayout;
import java.io.File;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.metagenes.swing.ProfileLinePanel;
import edu.mit.csail.cgs.metagenes.swing.ProfilePanel;
import edu.mit.csail.cgs.viz.paintable.PaintableScale;

public class MetaNonFrame{
	private Genome genome;
	private BinningParameters params;
	private MetaProfile profile;
	private MetaProfileHandler handler;
	private MetaUtils utils;
	private PaintableScale peakScale, lineScale;
	private ProfileLinePanel linePanel;
	private ProfilePanel panel;
	
	public MetaNonFrame(Genome g, BinningParameters bps, PointProfiler pp, boolean normalizedMeta) {
		peakScale = new PaintableScale(0.0, 1.0);
		lineScale = new PaintableScale(0.0, 1.0);
		
		genome = g;
		params = bps;
		handler = new MetaProfileHandler("MetaProfile", params, pp, normalizedMeta);
		profile = handler.getProfile();
		linePanel = new ProfileLinePanel(params, lineScale);
		profile.addProfileListener(linePanel);
		utils = new MetaUtils(genome);
		
		panel = new ProfilePanel(profile, peakScale);
	}
	public MetaNonFrame(Genome g, BinningParameters bps, PointProfiler pp, boolean normalizedMeta, double peakMax) {
		peakScale = new PaintableScale(0.0, peakMax);
		lineScale = new PaintableScale(0.0, 1.0);
		
		genome = g;
		params = bps;
		handler = new MetaProfileHandler("MetaProfile", params, pp, normalizedMeta);
		profile = handler.getProfile();
		linePanel = new ProfileLinePanel(params, lineScale);
		profile.addProfileListener(linePanel);
		utils = new MetaUtils(genome);
		
		panel = new ProfilePanel(profile, peakScale);
	}
	public void setColor(Color c){
		panel.updateColor(c);
		linePanel.updateColor(c);
	}
	public void setLinePanelColorQuanta(double [] q){
		linePanel.setLineColorQuanta(q);
	}
	public void setDrawColorBar(boolean c){
		linePanel.setDrawColorBar(c);
	}
	public void saveImages(String root){
		try {
			System.out.println("Saving images with root name: "+root);
			panel.saveImage(new File(root+"_profile.png"), 1250, 700);
			linePanel.saveImage(new File(root+"_lines.png"), linePanel.getPanelWidth(), linePanel.getPanelLength());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public void savePointsToFile(String root){
		String fileName = String.format("%s.points.txt", root);
		profile.saveToFile(fileName);
	}
	public MetaProfileHandler getHandler() { return handler; }
	public MetaUtils getUtils(){return utils;}
	public void clusterLinePanel(){linePanel.cluster();}
	public void setLineMin(double m){linePanel.setMinColorVal(m);}
	public void setLineMax(double m){linePanel.setMaxColorVal(m);}
	public void setLineThick(int t){linePanel.updateLineWeight(t);}
}
