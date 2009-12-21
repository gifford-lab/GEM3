package edu.mit.csail.cgs.viz.preferences;

import java.io.File;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;

public class PreferencesTest implements PreferencesListener {
	
	public static void main(String[] args) {
		System.out.println("Testing...");
		
		PreferencesTest test = new PreferencesTest();
		try {
			test.model.setValue("Genome", Organism.findGenome("mm8"));
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		test.model.setValue("FirstString", "foo");
		test.model.setValue("SecondInt", 3);
		test.model.setValue("ThirdBoolean", true);
		test.model.setValue("FourthBoolean", false);
		test.model.setValue("FifthFile", new File("test.txt"));
		
		test.showDialog();
	}
	
	private PreferencesModel.Default model;

	public PreferencesTest() { 
		model = new PreferencesModel.Default();
		model.addListener(this);
	}
	
	public void showDialog() { 
		PreferencesDialog dlg = new PreferencesDialog(model);
	}

	public void preferencesUpdateCanceled(PreferencesEvent evt) {
		System.out.println("Canceled.");
	}

	public void preferencesUpdated(PreferencesEvent evt) {
		System.out.println("Updated!");
		for(String key : model.getKeys()) { 
			System.out.println("\t" + key + " ==> " + model.getValue(key));
		}
	}
}
