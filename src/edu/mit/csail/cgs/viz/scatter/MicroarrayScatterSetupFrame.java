package edu.mit.csail.cgs.viz.scatter;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import edu.mit.csail.cgs.datasets.species.Genome;

/* Frame to configure a microarray data scatter plot.  This class is just a wrapper around the 
   MicroarrayScatterSetupPane.

   Command line usage:
   java MicroarrayScatterSetupFrame [--log] [--mvsa] [--species 'Mus musculus;mm8 [--exptone 'name;version;replicate'] [--expttwo 'name;version;replicate']]

*/

public class MicroarrayScatterSetupFrame extends JFrame implements ActionListener {

    private JButton ok, cancel;
    public MicroarrayScatterSetupPane pane;

    public MicroarrayScatterSetupFrame(Genome g) {
        pane = new MicroarrayScatterSetupPane(g);
        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());
        JPanel buttonpanel = new JPanel();
        ok = new JButton("OK");
        cancel = new JButton("Cancel");
        buttonpanel.add(ok);
        buttonpanel.add(cancel);
        ok.addActionListener(this);
        cancel.addActionListener(this);
        panel.add(pane, BorderLayout.CENTER);
        panel.add(buttonpanel, BorderLayout.SOUTH);
        getContentPane().add(panel);
        this.setSize(600,600);
        setVisible(true);
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == ok) {
            JFrame f = new JFrame();
            ScatterPanel p = new ScatterPanel("",
                                              pane.parse(),
                                              .1);
            f.getContentPane().add(p, BorderLayout.CENTER);
            f.setPreferredSize(new Dimension(800,800));
            f.setSize(800,800);
            f.setVisible(true);
            this.dispose();
        } else if (e.getSource() == cancel) {
            this.dispose();
        }


    }

}