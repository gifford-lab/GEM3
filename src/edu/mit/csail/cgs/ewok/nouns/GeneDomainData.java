package edu.mit.csail.cgs.ewok.nouns;

import java.util.*;
import java.text.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.ExonicGene;
import edu.mit.csail.cgs.datasets.species.Gene;

public class GeneDomainData extends RegionDomainData {
	
	private static NumberFormat nf;
	
	static { 
		nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(1);
		nf.setMinimumFractionDigits(1);
	}

	private Gene gene;
	
	private Region window;
	private RegionDomainData windowData;
	
	private Vector<Region> exons;
	private Vector<RegionDomainData> exonData;
	
	public GeneDomainData(Gene g, int win) { 
		super(g);
		gene = g;
		window = new Region(gene.getGenome(), gene.getChrom(),
				gene.getStart()-win, gene.getEnd()+win);
		windowData = new RegionDomainData(window);
		
		exons = new Vector<Region>();
		exonData = new Vector<RegionDomainData>();
		
		if(g instanceof ExonicGene) {
			ExonicGene eg = (ExonicGene)g;
			Iterator<Region> exitr = eg.getExons();
			while(exitr.hasNext()) {
				Region ex = exitr.next();
				exons.add(ex);
				exonData.add(new RegionDomainData(ex));
			}
		}
	}
	
	public void addDomain(Region dom) { 
		super.addDomain(dom);
		windowData.addDomain(dom);
		for(RegionDomainData rdd : exonData) { 
			rdd.addDomain(dom);
		}
	}
	
	public RegionDomainData getWindowData() { return windowData; }
	public RegionDomainData getExonData(int i) { return exonData.get(i); }
	public int getNumExons() { return exons.size(); }
	public Region getExon(int i) { return exons.get(i); }
	public Region getWindow() { return window; }
	public Gene getGene() { return gene; }
	
	public boolean isTSSCovered() {
		Point tss = new Point(gene.getGenome(), gene.getChrom(), gene.getTSS());
		for(Region unc : uncoveredRegions) { 
			if(unc.contains(tss)) { return false; }
		}
		return true;
	}
	
	public void printData() { 
		double cover = (double)getCoveredBP() / (double)getBP();
		double wincover = (double)windowData.getCoveredBP() / (double)windowData.getBP();
		
		System.out.println("Gene: " + gene.getID());
		System.out.println("\tLocation: " + gene.getLocationString());
		System.out.println("\t# Domains: " + getNumDomains());
		for(Region dom : domains) { 
			System.out.println("\t\t" + dom.getLocationString());
		}
		
		System.out.println("\t% Coverage: " + nf.format(cover));
		System.out.println("\t% Window Coverage: " + nf.format(wincover));

		if(exons.size() >= 1) { 
			System.out.println("\t# Exons: " + getNumExons());
			double ex_cover = (double)(exonData.get(0).getCoveredBP()) / 
				(double)(exonData.get(0).getBP());
			System.out.println("\t% 1st Exon Coverage: " + nf.format(ex_cover));
		}
		
		System.out.println();
	}
}
