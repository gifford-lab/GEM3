/*
 * Created on Aug 28, 2007
 */
package edu.mit.csail.cgs.tools.sequence;

import java.util.*;
import java.util.regex.*;
import java.util.logging.*;
import java.io.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.text.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * @author tdanford
 */
public class SequenceViewer extends JFrame implements CaretListener {
    
    public static void main(String[] args) { 
        SequenceViewer sv = new SequenceViewer();
    }
    
    private static Pattern regionPattern = Pattern.compile("([\\w\\d]+):([\\w\\d]+):(\\d+)-(\\d+)");

    private Logger logger;
    private Level logLevel;
    
    private SequenceGenerator seqgen;
    private Region currentRegion;
    private String currentSequence;
    private boolean currentStrand;

    private JTextPane sequenceArea;
    private StyledDocument seqDoc;
    private JTextField locField, matchField;
    private JButton retrieveSequence, matchSequence;
    private JLabel coordinateLabel;
    
    private JRadioButton forwardStrand, reverseStrand;
    private ButtonGroup strandGroup;
    
    public SequenceViewer() { 
        super("SequenceViewer");
        init("");
    }
    public SequenceViewer(String initialLocation) {
        super("SequenceViewer");
        init(initialLocation);
        fetchSequence(initialLocation);
    }
    public SequenceViewer(Region initialLocation) {
        super("SequenceViewer");
        String initial = getRegionString(initialLocation);
        init(initial);
        fetchSequence(initial);
    }
    

    private void init(String location) {
        logger = Logger.getLogger("edu.mit.csail.cgs.tools.sequence.SequenceViewer");
        logLevel = Level.INFO;
        logger.log(logLevel, "Starting SequenceViewer..");
        
        seqgen = new SequenceGenerator();
        currentRegion = null;
        currentSequence = null;
        currentStrand = true;
        
        sequenceArea = new JTextPane();
        sequenceArea.addCaretListener(this);
        seqDoc = sequenceArea.getStyledDocument();
        initStyles();
        
        locField = new JTextField(location);
        retrieveSequence = new JButton("Fetch");
        matchField = new JTextField();
        matchSequence = new JButton("Match");
        coordinateLabel = new JLabel(" ");
        
        forwardStrand = new JRadioButton("Forward");
        reverseStrand = new JRadioButton("Reverse");
        strandGroup = new ButtonGroup();
        strandGroup.add(forwardStrand);
        strandGroup.add(reverseStrand);
        forwardStrand.setSelected(true);
        
        JPanel strandPanel = new JPanel();
        strandPanel.setLayout(new BorderLayout());
        strandPanel.add(forwardStrand, BorderLayout.NORTH);
        strandPanel.add(reverseStrand, BorderLayout.SOUTH);
        
        JPanel locPanel = new JPanel();
        locPanel.setLayout(new BorderLayout());
        locPanel.add(locField, BorderLayout.CENTER);
        locPanel.add(retrieveSequence, BorderLayout.EAST);

        JPanel matchPanel = new JPanel();
        matchPanel.setLayout(new BorderLayout());
        matchPanel.add(matchField, BorderLayout.CENTER);
        matchPanel.add(matchSequence, BorderLayout.EAST);

        JPanel inputPanel = new JPanel();
        inputPanel.setLayout(new BorderLayout());
        inputPanel.add(strandPanel, BorderLayout.NORTH);
        inputPanel.add(locPanel, BorderLayout.SOUTH);

        JPanel seqPanel = new JPanel();
        seqPanel.setLayout(new BorderLayout());
        seqPanel.add(new JScrollPane(sequenceArea), BorderLayout.CENTER);
        seqPanel.setBorder(new TitledBorder("Chromosomal Sequence"));
        seqPanel.add(matchPanel, BorderLayout.NORTH);
        seqPanel.add(coordinateLabel, BorderLayout.SOUTH);
        
        Container c = (Container)getContentPane();
        c.setLayout(new BorderLayout());
        c.add(seqPanel, BorderLayout.CENTER);
        c.add(inputPanel, BorderLayout.SOUTH);
        
        setJMenuBar(createMenubar());
        
        retrieveSequence.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                fetchSequence(locField.getText());
            }
        });
        
        matchSequence.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                findAndMark(matchField.getText());
            }
        });
        
        forwardStrand.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent e) { 
        		setStrand(true);
        	}
        });
        
        reverseStrand.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent e) { 
        		setStrand(false);
        	}
        });
        
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(400, 600);
        setVisible(true);
    }
    
    public void setStrand(boolean value) {
    	if(value != currentStrand) {
    		currentSequence = reverseComplement(currentSequence);
    		currentStrand = value;
    		updateTextArea();
    	}
    }
    
    private JMenuBar createMenubar() { 
        JMenuBar bar = new JMenuBar();
        JMenu menu = null;
        JMenuItem item = null;
        final JFrame thisframe = this;
        
        bar.add(menu = new JMenu("File"));
        menu.add(item = new JMenuItem("Close"));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) { thisframe.dispose(); }
        });
        menu.add(item = new JMenuItem("Exit"));
        item.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { System.exit(0); }
        });
        
        bar.add(menu = new JMenu("Edit"));
        menu.add(item = new JMenuItem("Clear"));
        item.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                clearMarkings();
            }
        });
        
        return bar;
    }
    
    private void recordCoordinates(int dot, int mark) { 
        int start = Math.min(dot, mark), end = Math.max(dot, mark);
        start += currentRegion.getStart(); 
        end += currentRegion.getStart();
        
        String msg = start + "," + end;
        if(start == end) { 
            msg = String.valueOf(start);
        }
        SwingUtilities.invokeLater(new CoordinateSetter(msg));
    }

    public void caretUpdate(CaretEvent evt) {
        recordCoordinates(evt.getDot(), evt.getMark());
    }

    private void initStyles() { 
        //Initialize some styles.
        Style def = StyleContext.getDefaultStyleContext().
                        getStyle(StyleContext.DEFAULT_STYLE);

        Style regular = seqDoc.addStyle("regular", def);
        StyleConstants.setFontFamily(def, "Courier");
        StyleConstants.setFontSize(def, 14);

        Style s = seqDoc.addStyle("italic", regular);
        StyleConstants.setItalic(s, true);

        s = seqDoc.addStyle("bold", regular);
        StyleConstants.setBold(s, true);

        s = seqDoc.addStyle("red", regular);
        StyleConstants.setForeground(s, Color.red);
    }

    public String getRegionString(Region r) {
        return String.format("%s:%s:%d-%d",
                             r.getGenome().getVersion(),
                             r.getChrom(),
                             r.getStart(),
                             r.getEnd());
    }
    
    public void fetchSequence(String seqstring) { 
        logger.log(logLevel, "Retrieving sequence for " + seqstring);
        Matcher m = regionPattern.matcher(seqstring);
        if(m.matches()) { 
            String genomeName = m.group(1);
            String chromName = m.group(2);
            int start = Integer.parseInt(m.group(3));
            int end = Integer.parseInt(m.group(4));
            try {
                Genome genome = Organism.findGenome(genomeName);
                Region r = new Region(genome, chromName, start, end);
                String seq = seqgen.execute(r);
                seq = seq.toUpperCase();

                logger.log(logLevel, "Sequence: [" + seq + "]");
                
                currentRegion = r;
                currentSequence = seq;
                currentStrand = true;
                
                updateTextArea();
                
            } catch (NotFoundException e) {
                e.printStackTrace();
                logger.log(Level.WARNING, "Couldn't find genome \"" + genomeName + "\"");
            }
        } else { 
            logger.log(Level.WARNING, "Couldn't parse \"" + seqstring + "\"");
        }
    }
    
    public void findAndMark(String substr) { 
        Pattern p = Pattern.compile(substr);
        Matcher m = p.matcher(currentSequence);
        while(m.find()) { 
            int st = m.start(); 
            int end = m.end();
            SwingUtilities.invokeLater(new TextMarker(st, end, "red", false));
        }
    }
    
    private String reverseComplement(String str) { 
    	StringBuffer sb = new StringBuffer();
    	for(int i = str.length()-1; i >= 0; i--) { 
    		sb.append(complement(str.charAt(i)));
    	}
    	return sb.toString();
    }
    
    private char complement(char c) { 
    	switch(c) { 
    	case 'A': return 'T';
    	case 'a': return 't';
    	case 'T': return 'A';
    	case 't': return 'a';
    	case 'G': return 'C';
    	case 'g': return 'c';
    	case 'C': return 'G';
    	case 'c': return 'g';
    	default: return c;
    	}
    }
    
    public void clearMarkings() { 
        SwingUtilities.invokeLater(new TextMarker(0, currentSequence.length(), "regular", true));
    }
    
    private void updateTextArea() { 
        SwingUtilities.invokeLater(new Updater());
    }
    
    private class Updater implements Runnable { 
        public void run() { 
            try {
                seqDoc.remove(0, seqDoc.getLength());
                seqDoc.insertString(0, currentSequence, seqDoc.getStyle("regular"));
            } catch (BadLocationException e) {
                e.printStackTrace();
            }
        }
    }
    
    private class CoordinateSetter implements Runnable { 
        private String lbl; 
        public CoordinateSetter(String l) { lbl = l; }
        public void run() { 
            coordinateLabel.setText(lbl);
        }
    }
    
    private class TextMarker implements Runnable {
        private String stylename;
        private int start, end;
        private boolean rep;
        
        public TextMarker(int st, int ed, String s, boolean r) { 
            start = st; end = ed; stylename = s; rep = r; 
        }
        
        public void run() { 
            int length = end-start;
            logger.log(logLevel, "Marking " + start + "-" + end + " as " + stylename);
            seqDoc.setCharacterAttributes(start, end-start, seqDoc.getStyle(stylename), rep);
        }
    }

}
