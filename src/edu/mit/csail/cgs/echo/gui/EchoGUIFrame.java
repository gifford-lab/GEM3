/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

import edu.mit.csail.cgs.ewok.types.SelfDescribingConstant;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingExpander;
import edu.mit.csail.cgs.ewok.verbs.motifs.StringToGCCounts;

import edu.mit.csail.cgs.datasets.locators.BayesLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.echo.EchoConstant;
import edu.mit.csail.cgs.echo.Reverb;
import edu.mit.csail.cgs.echo.components.*;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.viz.components.SelectionEvent;

public class EchoGUIFrame extends JFrame {
    
    public static void main(String[] args) { 
        EchoGUIFrame f = new EchoGUIFrame();
        
        f.registerVerb(TabbedFileRegionSink.class);
        f.registerVerb(CountingSink.class);
        f.registerVerb(ListingSink.class);
        f.registerVerb(ChromRegionWrapper.class);
        f.registerVerb(RefGeneGenerator.class);
        f.registerVerb(BayesBindingGenerator.class);
        //f.registerVerb(BindingEventAnnotator.class);
        f.registerVerb(BindingExpander.class);
        f.registerVerb(UniqueFilter.class);
        f.registerVerb(GeneToPromoter.class);
        f.registerVerb(PairCombiner.class);
        f.registerVerb(RegionBreaker.class);
        f.registerVerb(RegionExpander.class);
        f.registerVerb(ConstantGenerator.class);
		f.registerVerb(StoringSink.class);
		f.registerVerb(DifferenceSink.class);
		f.registerVerb(ToStringMapper.class);
		f.registerVerb(CountingMapper.class);
        f.registerVerb(FirstElementMapper.class);
        f.registerVerb(IteratorIdentityMapper.class);
        f.registerVerb(IteratorToCollectionMapper.class);
        f.registerVerb(CollectionToIteratorMapper.class);
        f.registerVerb(FindClosest.class);
        f.registerVerb(MatchClosest.class);
        f.registerVerb(StrandedToStart.class);
        f.registerVerb(NonNullPairFilter.class);
        f.registerVerb(PairFlipper.class);
        f.registerVerb(PairFirst.class);
        f.registerVerb(PointToRegion.class);
        f.registerVerb(SequenceGenerator.class);
        f.registerVerb(StringToGCCounts.class);
        f.registerVerb(CountFractionSink.class);

		f.setVisible(true);
		f.pack();
    }
    
    private EchoGUI gui;
    private GUIVerbSelectionPanel verbSelection;
    private GUIConstantCreationPanel constCreation;
    private JButton runButton;
    private JMenuBar menuBar;
    
    public EchoGUIFrame() {
        super("Echo");
        
        gui = new EchoGUI();
        verbSelection = new GUIVerbSelectionPanel();
        constCreation = new GUIConstantCreationPanel();
        
        verbSelection.addEventListener(new Listener<SelectionEvent<GUIVerbSelector>>() {
            public void eventRegistered(SelectionEvent<GUIVerbSelector> e) {
                GUIVerbSelector sel = e.getFirstValue();
                SelfDescribingVerb v = sel.getVerb();
				addVerb(sel.getName(), v);
            } 
        });
        
        constCreation.addEventListener(new Listener<CreationEvent<SelfDescribingConstant>>() {
            public void eventRegistered(CreationEvent<SelfDescribingConstant> e) {
                SelfDescribingConstant k = e.getValue();
                if(k != null) { addConstant(k.toString(), k); }
            } 
        });
        
        runButton = new JButton("Run Program");
        
        Container c = (Container)getContentPane();
        c.setLayout(new BorderLayout());
        
        c.add(gui, BorderLayout.CENTER);
        c.add(verbSelection, BorderLayout.WEST);
        c.add(constCreation, BorderLayout.EAST);
        verbSelection.add(runButton, BorderLayout.SOUTH);
        
        runButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                gui.runBase();
            }
        });
        
        menuBar = createMenu();
        setJMenuBar(menuBar);
        
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
    
    private JMenuBar createMenu() { 
        JMenuBar bar = new JMenuBar();
        JMenu menu;
        JMenuItem item;
        
        bar.add(menu = new JMenu("File"));
        menu.add(item = new JMenuItem("Exit"));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) { 
		System.exit(0);
                //EchoGUIFrame.this.dispose();
            }
        });
        
        bar.add(menu = new JMenu("Edit"));
        menu.add(item = new JMenuItem("Clear"));
        item.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                gui.clear();
            }
        });
        
        return bar;
    }
    
    public void registerVerb(Class c) { 
    	verbSelection.addVerb(new GUIVerbSelector(c));
    }
    
	private void addVerb(String name, SelfDescribingVerb v) { 
		if(v instanceof Sink || v instanceof MultiSink) { 
			addEchoSink(name, new Reverb(v)); 
		} else if(v instanceof Generator) { 
			addEchoGenerator(name, new Reverb(v)); 
		} else { 
			addEchoReverb(name, new Reverb(v));
		}
	}

	private void addConstant(String name, SelfDescribingConstant c) { 
		addEchoConstant(name, new EchoConstant(c));
	}
    
    private void addEchoConstant(String name, EchoConstant inp) { 
        Rectangle rect = gui.getRandomRect(30, 30);
        EchoGUIConstant c = new EchoGUIConstant(gui, name, rect, inp);
        gui.addEchoComponent(c);
    }

    private void addEchoReverb(String name, Reverb inp) { 
        Rectangle rect = gui.getRandomRect(30, 30);
        EchoGUIReverb c = new EchoGUIReverb(gui, name, rect, inp);
        gui.addEchoComponent(c);
    }
    
    private void addEchoSink(String name, Reverb inp) { 
        Rectangle rect = gui.getRandomRect(30, 30);
        EchoGUISink c = new EchoGUISink(gui, name, rect, inp);
        gui.addEchoComponent(c);
    }
    
    private void addEchoGenerator(String name, Reverb inp) { 
        Rectangle rect = gui.getRandomRect(30, 30);
        EchoGUIGenerator c = new EchoGUIGenerator(gui, name, rect, inp);
        gui.addEchoComponent(c);
    }
        

}
