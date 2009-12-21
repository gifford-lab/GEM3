package edu.mit.csail.cgs.viz;
import edu.mit.csail.cgs.utils.*;

import java.lang.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.imageio.*;

public class FunctionDisplay { 

    public static void main(String[] args) { 
	try { 
	    int start = -1;
	    int stop = -1;
	    if(args.length > 1) { 
		start = Integer.parseInt(args[1]);
		stop = Integer.parseInt(args[2]);
	    }

	    FunctionArrayDisplayFrame fadf = 
		new FunctionArrayDisplayFrame(new File(args[0]), start, stop);
	} catch(IOException ie) { 
	    ie.printStackTrace(System.err);
	}
    }

    public static class FunctionArrayDisplayFrame extends JFrame { 

	private FunctionArrayDisplay fArrayPanel;
	private String[] fInfoLines;
	private boolean fShowStats;

	public FunctionArrayDisplayFrame(String name, FunctionDisplay fd, boolean showStats) { 
	    super(name);
	    LinkedList lst = new LinkedList();
	    lst.addLast(fd);

	    FunctionArrayDisplay fad = 
		new FunctionArrayDisplay(1, 1, lst);
	    fShowStats = showStats;
	    fad.setShowFunctionStats(fShowStats);
	    fArrayPanel = fad;
	    fInfoLines = new String[0];
	    Container c = (Container)getContentPane();
	    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	    c.setLayout(new BorderLayout());
	    c.add(fad, BorderLayout.CENTER);

	    setJMenuBar(createMenuBar());
	    
	    setSize(300, 300);
	    setVisible(true);
	}

	public FunctionArrayDisplayFrame(File input, int start, int stop) throws IOException { 
	    super("Function Display Array");
	    BufferedReader br = new BufferedReader(new FileReader(input));
	    String line;
	    StringTokenizer st;
	    line = br.readLine();
	    st = new StringTokenizer(line);
	    int lines = Integer.parseInt(st.nextToken());
	    int rows = Integer.parseInt(st.nextToken());
	    int cols = Integer.parseInt(st.nextToken());
	    System.out.println("(" + lines + "," + rows + "," + cols + ")");

	    if(start != -1 && stop != -1) { 
		int diff = stop - start + 1;
		if(diff <= cols) { 
		    rows = 1; cols = diff;
		} else { 
		    rows = (diff / cols) + 1;
		}
	    }
	    
	    line = br.readLine();
	    st = new StringTokenizer(line, ":");
	    fInfoLines = new String[st.countTokens()];
	    for(int i = 0; i < fInfoLines.length; i++) { 
		fInfoLines[i] = st.nextToken();
	    }

	    Color[] funcColors = 
		{Color.red, Color.blue, Color.red, Color.blue, Color.black, Color.green};
	    Color axesColor = null;
	    LinkedList lst = new LinkedList();
	    
	    for(int i = 0; i < lines; i++) { 
		line = br.readLine();
		if(start != -1 && stop != -1 && (i < start || i > stop)) { 
		    continue;
		} 

		st = new StringTokenizer(line, ",");

		String name = st.nextToken();
		if(name.trim().equals("*")) { 
		    axesColor = Color.black;
		    name = st.nextToken();
		} else { 
		    axesColor = null;
		}

		double lx = Double.parseDouble(st.nextToken());
		double hx = Double.parseDouble(st.nextToken());
		double ly = Double.parseDouble(st.nextToken());
		double hy = Double.parseDouble(st.nextToken());
		FunctionDisplay fd = 
		    new FunctionDisplay(name, 10, 10, false, axesColor, lx, hx, ly, hy);
		int ci = 0;
		while(st.hasMoreTokens()) { 
		    String nextTok = st.nextToken();
		    if(nextTok.startsWith("pts") || nextTok.startsWith("lines")) {
			StringTokenizer ptTok = new StringTokenizer(nextTok);
			ptTok.nextToken();
			double[] vals = new double[ptTok.countTokens()];
			for(int p = 0; p < vals.length; p++) { 
			    vals[p] += Double.parseDouble(ptTok.nextToken());
			}
			
			fd.addPoints(vals, funcColors[ci], nextTok.startsWith("lines"));
		    } else {
			fd.addFunction(createFunction(nextTok), 
				       funcColors[ci]);
		    }
		    if(ci < funcColors.length-1) { ci++; }
		}
		lst.addLast(fd);
	    }
	    br.close();

	    FunctionArrayDisplay fad = new FunctionArrayDisplay(rows, cols, lst);
	    fArrayPanel = fad;
	    fShowStats = false;
	    fad.setShowFunctionStats(fShowStats);

	    Container c = (Container)getContentPane();
	    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    c.setLayout(new BorderLayout());
	    c.add(fad, BorderLayout.CENTER);

	    setJMenuBar(createMenuBar());
	    	    
	    setSize(500, 500);
	    setVisible(true);
	}

	private JMenuBar createMenuBar() { 
	    JMenuBar jmb = new JMenuBar();
	    JMenu menu; JMenuItem item;

	    jmb.add(menu = new JMenu("File"));
	    menu.add(item = new JMenuItem(createShowInfoAction()));
	    menu.add(item = new JMenuItem(createFlipShowStatsAction()));
	    menu.add(item = new JMenuItem(fArrayPanel.createSaveImageAction()));
	    menu.add(item = new JMenuItem("Exit"));
	    item.addActionListener(new ActionListener() { 
		    public void actionPerformed(ActionEvent e) { 
			System.exit(0);
		    }
		});

	    return jmb;
	}

	public void setShowStats(boolean v) { 
	    fShowStats = v;
	    fArrayPanel.setShowFunctionStats(v);
	}

	private Action createFlipShowStatsAction() { 
	    return new AbstractAction("Toggle Show Stats") { 
		    public void actionPerformed(ActionEvent e) { 
			setShowStats(!fShowStats);
		    }
		};
	}

	private Action createShowInfoAction() { 
	    return new AbstractAction("Show Info...") { 
		    public void actionPerformed(ActionEvent e) { 
			StringSetDialog dlg = new StringSetDialog(fInfoLines, false);
			dlg.setVisible(true);
		    }
		};
	}
    }

    public static class FunctionArrayDisplay extends JPanel { 
	
	private int fRows, fCols;
	private int fWSpacing, fHSpacing;
	private int fFDWidth, fFDHeight;
	private FunctionDisplay[][] fArray;
	private ComponentListener fResizeListener;

	public FunctionArrayDisplay(int rows, int cols, 
				    Collection fdCollect) { 
	    super();
	    fRows = rows; fCols = cols;
	    fWSpacing = 5;
	    fHSpacing = 5;

	    fArray = new FunctionDisplay[fRows][fCols];
	    Iterator itr = fdCollect.iterator();
	    for(int i = 0; i < fRows; i++) { 
		for(int j = 0; j < fCols; j++) { 
		    if(itr.hasNext()) { 
			fArray[i][j] = (FunctionDisplay)itr.next();
		    } else { 
			fArray[i][j] = null;
		    }
		}
	    }
	    
	    recalcDims();
	    fResizeListener = new ComponentAdapter() { 
		    public void componentResized(ComponentEvent ce) { 
			recalcDims();
			repaint();
		    }
		};
	    addComponentListener(fResizeListener);

	    addMouseListener(new MouseAdapter() { 
		    public void mouseClicked(MouseEvent e) { 
			int row = e.getY() / (fHSpacing + fFDHeight);
			int col = e.getX() / (fWSpacing + fFDWidth);
			if(fArray[row][col] != null) { 
			    String name = fArray[row][col].getName();
			    FunctionDisplay fd = fArray[row][col].cloneDisplay();
			    FunctionArrayDisplayFrame fadf = 
				new FunctionArrayDisplayFrame(name, fd, true);
			}
		    }
		});
	}

	public void setShowFunctionStats(boolean v) { 
	    for(int i = 0; i < fRows; i++) { 
		for(int j = 0; j < fCols; j++) { 
		    if(fArray[i][j] != null) { 
			fArray[i][j].setShowFunctionStats(v);
		    }
		}
	    }
	    
	    repaint();
	}

	public Action createSaveImageAction() { 
	    return new AbstractAction("Save As Image...") { 
		    public void actionPerformed(ActionEvent e) { 
			String pwdName = System.getProperty("user.dir");
			JFileChooser chooser;
			if(pwdName != null) { 
			    chooser = new JFileChooser(new File(pwdName));
			} else {
			    chooser = new JFileChooser();
			}
			
			int v = chooser.showOpenDialog(FunctionArrayDisplay.this);
			if(v == JFileChooser.APPROVE_OPTION) { 
			    File f = chooser.getSelectedFile();
			    try {
				saveImage(f);
			    } catch(IOException ie) {
				ie.printStackTrace(System.err);
			    }
			}
			
		    }
		};
	}

	public void saveImage(File f) throws IOException { 
	    int w = getWidth();
	    int h = getHeight();
	    BufferedImage im = 
		new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics g = im.getGraphics();
	    g.setColor(Color.white);
	    g.fillRect(0, 0, w, h);
	    paintComponent(g);
	    ImageIO.write(im, "jpg", f);
	}
	
	protected void paintComponent(Graphics g) { 
	    super.paintComponent(g);
	    for(int i = 0; i < fRows; i++) { 
		for(int j = 0; j < fCols; j++) { 
		    int llx = fWSpacing + (fFDWidth + fWSpacing) * j;
		    int lly = fFDHeight + fHSpacing + 
			(fFDHeight + fHSpacing) * i;
		    if(fArray[i][j] != null) { 
			fArray[i][j].draw(g, llx, lly, llx + fFDWidth, lly - fFDHeight);
		    }
		}
	    }
	}

	private void recalcDims() { 
	    int w = getWidth();
	    int h = getHeight();
	    fFDWidth = (w - (fWSpacing * (fCols+1))) / fCols;
	    fFDHeight = (h - (fHSpacing * (fRows+1))) / fRows;
	    for(int i = 0; i < fRows; i++) { 
		for(int j = 0; j < fCols; j++) { 
		    if(fArray[i][j] != null) { 
			fArray[i][j].setSize(fFDWidth, fFDHeight);
		    }
		}
	    }
	}
    }
    
    public static RealValuedFunction createFunction(String desc) { 
	StringTokenizer st = new StringTokenizer(desc);
	String first = st.nextToken();

	if(first.equals("gaussian")) { 
	    double mean = 0.0; 
	    double var = 1.0;
	    if(st.hasMoreTokens()) { 
		mean = Double.parseDouble(st.nextToken());
	    }
	    if(st.hasMoreTokens()) { 
		var = Double.parseDouble(st.nextToken());
	    }
	    
	    return new GaussianFunction(mean, var);
	}

	return null;
    }

    private String fName;
    private int fHeight, fWidth;
    private double fLoX, fHiX, fLoY, fHiY;
    private double fXScale, fYScale;  // x-units per pixel
    private boolean fShowFunctionStats;

    private Vector fFuncs;
    private Vector fColors;
    private Color fAxesColor;

    private Vector fRendered;

    private Vector fPoints;
    private Vector fPointsColors;

    public FunctionDisplay(String name, 
			   int w, int h, boolean showStats,
			   Color axesColor,
			   double lx, double hx, 
			   double ly, double hy) { 
	fName = name;
	fShowFunctionStats = showStats;
	fWidth = w; fHeight = h;

	fFuncs = new Vector();
	fColors = new Vector();
	fRendered = new Vector();
	fPoints = new Vector();
	fPointsColors = new Vector();

	fAxesColor = axesColor;

	setExtents(lx, hx, ly, hy);
    }

    public FunctionDisplay cloneDisplay() { 
	FunctionDisplay fd = new FunctionDisplay(fName, fWidth, fHeight,
						 fShowFunctionStats, fAxesColor,
						 fLoX, fHiX, fLoY, fHiY);
	for(int i = 0; i < fFuncs.size(); i++) { 
	    fd.addFunction((RealValuedFunction)fFuncs.get(i), (Color)fColors.get(i));
	}
	
	for(int i = 0; i < fPoints.size(); i++) { 
	    RenderedPoints rp = (RenderedPoints)fPoints.get(i);
	    fd.addPoints(rp.getValues(), (Color)fPointsColors.get(i), rp.getLine());
	}

	return fd;
    }

    public String getName() { return fName; }
    public void setShowFunctionStats(boolean v) { fShowFunctionStats = v; }
    
    public void addFunction(RealValuedFunction f, Color c) { 
	fFuncs.add(f);
	fColors.add(c);
	RenderedFunction rf = new RenderedFunction(f);
	fRendered.add(rf);
    }

    public void addPoints(double[] vals, Color c, boolean line) { 
	fPoints.add(new RenderedPoints(vals, line));
	fPointsColors.add(c);
    }

    public void setSize(int width, int height) { 
	fHeight = height; fWidth = width;
	fXScale = (fHiX - fLoX) / (double)fWidth;
	fYScale = (fHiY - fLoY) / (double)fHeight;
	for(int i = 0; i < fRendered.size(); i++) { 
	    ((RenderedFunction)fRendered.get(i)).resampleFunction();
	}

	for(int i = 0; i < fPoints.size(); i++) { 
	    ((RenderedPoints)fPoints.get(i)).resamplePoints();
	}
    }

    public void setExtents(double loX, double hiX, 
			   double loY, double hiY) { 
	fLoX = loX; fHiX = hiX;
	fLoY = loY; fHiY = hiY;
	fXScale = (fHiX - fLoX) / (double)fWidth;
	fYScale = (fHiY - fLoY) / (double)fHeight;

	for(int i = 0; i < fRendered.size(); i++) { 
	    ((RenderedFunction)fRendered.get(i)).resampleFunction();
	}
	for(int i = 0; i < fPoints.size(); i++) { 
	    ((RenderedPoints)fPoints.get(i)).resamplePoints();
	}
    }

    public void draw(Graphics g, int llx, int lly, int urx, int ury) { 
	drawAxes(g, llx, lly);
	for(int i = 0; i< fFuncs.size(); i++) { 
	    g.setColor((Color)(fColors.get(i)));
	    ((RenderedFunction)fRendered.get(i)).renderFunction(g, llx, lly);
	}

	for(int i = 0; i < fPoints.size(); i++) { 
	    g.setColor((Color)fPointsColors.get(i));
	    ((RenderedPoints)fPoints.get(i)).renderPoints(g, llx, lly, urx, ury);
	}
    }

    public void drawAxes(Graphics g, int llx, int lly) { 
	Graphics2D g2 = (Graphics2D)g;
	Stroke oldStroke = g2.getStroke();

	// draws the frame
	if(fAxesColor != null) { 
	    g2.setColor(fAxesColor);
	    g2.drawRect(llx, lly-fHeight, fWidth, fHeight);
	    g2.setStroke(new BasicStroke((float)3.0));
	}
	g2.setStroke(oldStroke);

	// draws the origin-lines, if visible
	g2.setColor(Color.gray);
	if(fLoX < 0.0 && fHiX > 0.0) { 
	    int zeroX = (int)Math.round((-fLoX) / fXScale);
	    g2.drawLine(llx + zeroX, lly - fHeight, llx + zeroX, lly);
	}

	if(fLoY < 0.0 && fHiY > 0.0) { 
	    int zeroY = (int)Math.round((-fLoY) / fYScale);
	    g2.drawLine(llx, lly - zeroY, llx + fWidth, lly - zeroY);
	}
    }

    public static class GaussianFunction implements RealValuedFunction { 

	private String fName;
	private double fMean, fVar;
	private double fNorm;

	public GaussianFunction() { 
	    double m = 0.0; 
	    double v = 1.0;
	    fName = "Gaussian(" + m + "," + v + ")";
	    fMean = m; fVar = v;
	    fNorm = 1.0 / Math.sqrt(2.0 * Math.PI * fVar);
	}

	public GaussianFunction(double m, double v) { 
	    fName = "Gaussian(" + m + "," + v + ")";
	    fMean = m; fVar = v;
	    fNorm = 1.0 / Math.sqrt(2.0 * Math.PI * fVar);
	}

	public String getName() { return fName; }
	public double getMean() { return fMean; }
	public double getVar() { return fVar; }

	public double eval(double x) { 
	    double diff = x - fMean;
	    double expt = - (diff * diff) / (2.0 * fVar);
	    return fNorm * Math.exp(expt);
	}
    }

    private class RenderedPoints { 

	private double[] fValues;
	private int[] fOffsets;
	private int fDim;
	private boolean fLines;
	
	public RenderedPoints(double[] vals) { 
	    fValues = (double[])vals.clone();
	    fOffsets = new int[vals.length];
	    fDim = 2;
	    fLines = false;
	}

	public RenderedPoints(double[] vals, boolean lines) { 
	    fValues = (double[])vals.clone();
	    fOffsets = new int[vals.length];
	    fDim = 2;
	    fLines = lines;
	}

	public double[] getValues() { return fValues; }
	public boolean getLine() { return fLines; }

	public void resamplePoints() { 
	    for(int i = 0; i < fValues.length; i++) { 
		fOffsets[i] = (int)Math.round((fValues[i] - fLoX) / fXScale);
		if(fOffsets[i] < 0 || fOffsets[i] >= fWidth) { fOffsets[i] = -1; }
	    }
	}

	public void renderPoints(Graphics g, int llx, int lly, int urx, int ury) {
	    int zeroY = 0;
	    if(fLoY < 0.0 && fHiY > 0.0) { 
		zeroY = (int)Math.round(-fLoY / fYScale);
	    }

	    for(int i = 0; i < fOffsets.length; i++) { 
		if(fOffsets[i] >= 0) { 
		    if(fLines) { 
			g.drawRect(llx + fOffsets[i], ury, fDim, lly - ury);
		    } else { 
			g.drawRect(llx + fOffsets[i] - fDim, lly - zeroY - fDim, 
				   fDim * 2, fDim * 2);
		    }
		}
	    }
	}
    }

    private class RenderedFunction { 

	private RealValuedFunction fBaseFunction;
	private int[] fPixYArray;
	
	public RenderedFunction(RealValuedFunction base) { 
	    fBaseFunction = base;
	    resampleFunction();
	}

	public void resampleFunction() { 
	    int array_length = fWidth;
	    if(array_length < 0) { array_length = 0; }
	    fPixYArray = new int[array_length];
	    for(int i = 0; i < fPixYArray.length; i++) { 
		double x = fLoX + ((double)i * fXScale);
		double y = fBaseFunction.eval(x);
		fPixYArray[i] = (int)Math.round((y - fLoY) / fYScale);
		if(fPixYArray[i] < 0 || fPixYArray[i] >= fHeight) { 
		    fPixYArray[i] = -1;
		}
	    }
	}

	public void renderFunction(Graphics g, int llx, int lly) { 
	    int lx, ly, cx, cy;
	    lx = ly = cx = cy = -1;
	    for(int i = 0; i < fPixYArray.length; i++) { 
		cx = i; cy = fPixYArray[i];
		if(ly != -1 && cy != -1) { 
		    g.drawLine(lx + llx, lly - ly, 
			       cx + llx, lly - cy);
		}

		lx = cx; ly = cy;
	    }

	    if(fShowFunctionStats && fBaseFunction instanceof GaussianFunction) { 
		GaussianFunction gf = (GaussianFunction)fBaseFunction;
		double mean = gf.getMean();
		double var = gf.getVar();
		double heightAtMean = 1.0 / Math.sqrt(Math.PI * 2.0 * var);
		int mx = (int)Math.round((mean - fLoX) / fXScale);
		int my = (int)Math.round((heightAtMean - fLoY) / fYScale);
		int dim = 2;
		g.drawRect(llx + mx - dim, lly - my - dim, dim * 2, dim * 2);
		g.drawString("(" + mean + ")", llx + mx + 2, lly - my);
	    }
	}
    }
}
