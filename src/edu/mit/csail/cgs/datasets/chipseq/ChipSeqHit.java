package edu.mit.csail.cgs.datasets.chipseq;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class ChipSeqHit extends StrandedRegion {

  private ChipSeqAlignment align;
  private double weight = 1.0;

  public ChipSeqHit(Genome g, String chrom, int start, int end, char strand, ChipSeqAlignment align, double weight) {
    super(g,chrom,start,end,strand);
    this.align = align;
    this.weight = weight;
  }

  public ChipSeqAlignment getAlignment() { return align; }

  public double getWeight() {
    return weight;
  }


  public void setWeight(double weight) {
    this.weight = weight;
  }

  public String toString() { 
    return String.format("%s:%c", getLocationString(), getStrand());
  }

  public ChipSeqHit extendHit(int ext) { 
    if(getStrand() == '+') { 
      return new ChipSeqHit(getGenome(), getChrom(), getStart(), getEnd() + ext, getStrand(), align, weight);
    } else { 
      return new ChipSeqHit(getGenome(), getChrom(), getStart() - ext, getEnd(), getStrand(), align, weight);
    }
  }
  public ChipSeqHit shiftExtendHit(int ext, int shift) { 
    if(getStrand() == '+') { 
      return new ChipSeqHit(getGenome(), getChrom(), getStart()+shift-(ext/2), getEnd()+shift+(ext/2), getStrand(), align, weight);
    } else { 
      return new ChipSeqHit(getGenome(), getChrom(), getStart()-shift-(ext/2), getEnd()-shift+(ext/2), getStrand(), align, weight);
    }
  }

  public int hashCode() { 
    int code = 17;
    code += super.hashCode(); code *= 37;
    code += align.hashCode(); code *= 37;
    return code; 
  }

  public boolean equals(Object o) { 
    if(!(o instanceof ChipSeqHit)) { return false; }
    ChipSeqHit d = (ChipSeqHit)o;
    if(!super.equals(d)) { return false; }
    if(!align.equals(d.align)) { return false; }
    return true;
  }  
}
