package edu.mit.csail.cgs.warpdrive.components;

import edu.mit.csail.cgs.warpdrive.WarpOptions;
import edu.mit.csail.cgs.datasets.species.Genome;

public interface PainterContainer {
    public void addPaintersFromOpts(WarpOptions opts);
    public Genome getGenome();
}
