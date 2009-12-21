package edu.mit.csail.cgs.warpdrive.components;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.SQLException;
import java.util.*;

import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSlider;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYDataset;

import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDataset;
import edu.mit.csail.cgs.datasets.chipchip.ExptNameVersion;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.warpdrive.WarpOptions;

public class ChipChipTwoArrayScatterFrame extends JFrame {

	JFreeChart chart;
	Color color;

	public ChipChipTwoArrayScatterFrame() throws SQLException, NotFoundException {
		final ChipChipTwoArrayScatterFrame thisFrame = this;
		JMenuBar menubar = new JMenuBar();
		JMenu filemenu, colormenu;
		filemenu = new JMenu("File");
		JMenuItem item;
		filemenu.add((item = new JMenuItem("Close")));
		item.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				thisFrame.dispose();
			}});
		filemenu.add((item = new JMenuItem("Exit")));
		item.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.exit(0);
			}});
		menubar.add(filemenu);

		colormenu = new JMenu("Color");
		colormenu.add((item = new JMenuItem("Set Color")));
		item.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				thisFrame.chooseColor();
			}});
		menubar.add(colormenu);

		setJMenuBar(menubar);
	}

	public void addChart(JFreeChart chart) {
		this.chart = chart;
		ChartPanel panel = new ChartPanel(chart);

		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(panel,BorderLayout.CENTER);
		setSize(1000,1000);  
	}

	public void chooseColor() {
		final ChipChipTwoArrayScatterFrame thisFrame = this;
		final JFrame frame = new JFrame("Plot Color");
		frame.getContentPane().setLayout(new BorderLayout());

		final JSlider slider = new JSlider(0,
				255,
				color.getAlpha());
		frame.getContentPane().add(slider,
				BorderLayout.NORTH);


		final JColorChooser chooser = new JColorChooser(color);
		frame.getContentPane().add(chooser,
				BorderLayout.CENTER);

		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new GridBagLayout());
		Dimension buttonSize = new Dimension(40,25);
		JButton ok = new JButton("OK");
		JButton cancel = new JButton("Cancel");
		ok.setMaximumSize(buttonSize);
		cancel.setMaximumSize(buttonSize);
		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color t = chooser.getColor();
				color = new Color(t.getRed(),
						t.getGreen(),
						t.getBlue(),
						slider.getValue());
				Plot plot= chart.getPlot();
				XYItemRenderer renderer = ((XYPlot)plot).getRenderer();
				renderer.setSeriesPaint(0, color);
				frame.dispose();
			}});
		cancel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frame.dispose();
			}});
		buttonPanel.add(ok);
		buttonPanel.add(cancel);

		frame.getContentPane().add(buttonPanel,
				BorderLayout.SOUTH);
		frame.setSize(500,500);
		frame.setVisible(true);
	}

	public JFreeChart generateChart(SQLData data1,
			SQLData data2,
			String key,
			Genome genome,
			double sample) throws NotFoundException {
		ArrayList<Double> ratios1 = new ArrayList<Double>();
		ArrayList<Double> ratios2 = new ArrayList<Double>();
		ChromosomeGenerator<Genome> chromgen = new ChromosomeGenerator<Genome>();
		Iterator<Region> chromiter = chromgen.execute(genome);

		Map<Integer,Double> values1 = new HashMap<Integer,Double>();
		Map<Integer,Double> values2 = new HashMap<Integer,Double>();

		while (chromiter.hasNext()) {
			Region chromregion = chromiter.next();
			data1.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
			data2.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());

			for(int i = 0; i < data1.getCount(); i++) {
				double sum = 0.0;
				int count = 0;
				if (data1.getReplicates(i) > 1) {
					//System.err.print("More than one replicate in exp1: "+chromregion.getChrom()+" "+i+" expt ids:");
					for (int j=0; j<data1.getReplicates(i); j++) {
						//System.err.print(" "+data1.getExptID(i, j));
					}
					//System.err.println();
					continue;
				}
				for(int j = 0; j < data1.getReplicates(i); j++) { 
					sum += data1.getRatio(i, j);
					count += 1;
				}

				values1.put(data1.getPos(i), sum / (double)count);
			}

			for(int i = 0; i < data2.getCount(); i++) {
				double sum = 0.0;
				int count = 0;
				if (data2.getReplicates(i) > 1) {
					//System.err.print("More than one replicate in exp2: "+chromregion.getChrom()+" "+i+" expt ids:");
					for (int j=0; j<data2.getReplicates(i); j++) {
						//System.err.print(" "+data2.getExptID(i, j));
					}
					//System.err.println();
					continue;
				}
				for(int j = 0; j < data2.getReplicates(i); j++) { 
					sum += data2.getRatio(i, j);
					count += 1;
				}

				values2.put(data2.getPos(i), sum / (double)count);
			}
		}

		SetTools<Integer> tools = new SetTools<Integer>();
		Set<Integer> commonPos = tools.intersection(values1.keySet(), values2.keySet());
		System.err.println(commonPos.size()+" probes in common");
		DefaultXYDataset dataset = new DefaultXYDataset();
		double[][] values = new double[2][commonPos.size()];
		int posind = 0;
		for(int pos : commonPos) {
			double v1 = values1.get(pos);
			double v2 = values2.get(pos);
			//values[0][posind] = 0.5 * (Math.log(v1) + Math.log(v2)); // A
			//values[1][posind] = Math.log(v1) - Math.log(v2); // M
			values[0][posind] = Math.log(v2);
			values[1][posind] = Math.log(v1);
			posind++;
		}

		dataset.addSeries(key, values);
		chart = ChartFactory.createScatterPlot(key,"V2","V1",dataset,PlotOrientation.VERTICAL,false,false,false);
		chart.setAntiAlias(true);

		XYPlot plot= chart.getXYPlot();
		
		ValueAxis dAxis = plot.getDomainAxis();
		ValueAxis rAxis = plot.getRangeAxis();
		
		double maxbound = Math.max(dAxis.getUpperBound(), rAxis.getUpperBound());
		double minbound = Math.min(dAxis.getLowerBound(), rAxis.getLowerBound());
		
		dAxis.setUpperBound(maxbound);
		rAxis.setUpperBound(maxbound);
		dAxis.setLowerBound(minbound);
		rAxis.setLowerBound(minbound);
		
		plot.addAnnotation(new XYLineAnnotation(minbound, minbound, maxbound, maxbound));

		color = new Color((float)0, 
				(float)0,
				(float)1.0,
				(float).1);

		XYItemRenderer renderer = ((XYPlot)plot).getRenderer();
		renderer.setSeriesPaint(0, color);
//		renderer.setBaseStroke(new BasicStroke((float).5));
//		renderer.setSeriesStroke(0,new BasicStroke((float).5));
		renderer.setSeriesShape(0,new Rectangle(1,1));
		
//		plot.addRangeMarker(new ValueMarker(0.0));

		return chart;
	}

	public static void main(String args[]) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		String line;
		String[] split;
		ArrayList<ExptNameVersion> expts = new ArrayList<ExptNameVersion>();
		while((line = br.readLine()) != null) {
			split = line.split("\t");
			expts.add(new ExptNameVersion(split[0],split[1],split[2]));
		}
		
		Genome g = new Genome("Mus musculus","mm8");
		ChipChipDataset agilent = new ChipChipDataset(g);
		for (int i=0; i<expts.size(); i++) {
			for (int j=i+1; j<expts.size(); j++) {
				ChipChipData data1 = agilent.getData(expts.get(i));
				ChipChipData data2 = agilent.getData(expts.get(j));
				if (data1 instanceof SQLData && data2 instanceof SQLData) {
					ChipChipTwoArrayScatterFrame f = new ChipChipTwoArrayScatterFrame();
					JFreeChart chart = f.generateChart((SQLData) data1, (SQLData) data2,
							expts.get(i).toString() + " " + expts.get(j).toString(),
							g,
							.2);
					ChartUtilities.saveChartAsJPEG(new File(args[1]+"/scatplot"+i+"_"+j+".jpg"), chart, 1000, 1000);
					System.err.println("image saved to: "+args[1]+"/scatplot"+i+"_"+j+".jpg");
					System.err.println("Created for " + expts.get(i).toString() + " and " + expts.get(j).toString());
				}
			}
		}
		/*
		WarpOptions options = WarpOptions.parseCL(args);
		Genome g = new Genome(options.species,options.genome);
		ChipChipDataset agilent = new ChipChipDataset(g);
		ExptNameVersion env1 = options.agilentdata.get(0);
		ExptNameVersion env2 = options.agilentdata.get(1);
		ChipChipData data1 = agilent.getData(env1);
		ChipChipData data2 = agilent.getData(env2);
		if (data1 instanceof SQLData && data2 instanceof SQLData) {
			ChipChipTwoArrayScatterFrame f = new ChipChipTwoArrayScatterFrame();
			JFreeChart chart = f.generateChart((SQLData) data1, (SQLData) data2,
					env1.toString() + " " + env2.toString(),
					g,
					.2);
			if (options.saveimage) {
				ChartUtilities.saveChartAsJPEG(new File(options.filename), chart, 1000, 1000);
				System.err.println("image saved to: "+options.filename);
			} else {
				f.addChart(chart);
				f.setVisible(true);
			}
			System.err.println("Created for " + env1.toString() + " and " + env2.toString());
		}
		*/
	}

}
