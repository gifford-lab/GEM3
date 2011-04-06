package edu.mit.csail.cgs.metagenes;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.metagenes.swing.MultiProfilePanel;
import edu.mit.csail.cgs.viz.paintable.PaintableScale;

public class MetaNonFrameMultiSet{
	private Genome genome;
	private BinningParameters params;
	private List<Profile> profiles = new ArrayList<Profile>();
	private List<MetaProfileHandler> handlers = new ArrayList<MetaProfileHandler>();
	private MetaUtils utils;
	private PaintableScale peakScale;
	private MultiProfilePanel panel;
	
	public MetaNonFrameMultiSet(List<String> setNames, Genome g, BinningParameters bps, PointProfiler pp, boolean normalizedMeta) {
		peakScale = new PaintableScale(0.0, 1.0);
		
		genome = g;
		params = bps;
		for(int i=0; i<setNames.size(); i++){
			MetaProfileHandler handler = new MetaProfileHandler(setNames.get(i), params, pp, normalizedMeta);
			handlers.add(handler);
			profiles.add(handler.getProfile());
		}
		utils = new MetaUtils(genome);
		
		panel = new MultiProfilePanel(profiles, peakScale);
	}
		
	public void saveImage(String root){
		try {
			System.out.println("Saving images with root name: "+root);
			panel.saveImage(new File(root+"_profile.png"), 1200, 700);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public void savePointsToFile(String root){
		String fileName = String.format("%s.points.txt", root);
		try {
			FileWriter fout = new FileWriter(fileName);
			int start = (-1*(params.getWindowSize()/2))+params.getBinSize()/2;
			int step = params.getWindowSize()/params.getNumBins();
			
			fout.write("Pos");
			for(Profile q : profiles)
				fout.write("\t"+q.getName());
			fout.write("\n");
			int i=0;
			for(int  k= start; k<params.getWindowSize(); k+=step){
				fout.write(k);
				for(Profile q : profiles)
					fout.write("\t"+q.value(i));
				fout.write("\n");
				i++;
			}			
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public List<MetaProfileHandler> getHandlers() { return handlers; }
	public MetaUtils getUtils(){return utils;}
}
