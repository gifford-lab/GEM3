package edu.mit.csail.cgs.datasets.chipseq;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.projects.readdb.SingleHit;

public class ChiaPetLoader extends ChipSeqLoader {

	public ChiaPetLoader() throws SQLException, IOException {
		super();
		// TODO Auto-generated constructor stub
	}

	public ChiaPetLoader(boolean openClient) throws SQLException, IOException {
		super(openClient);
	}

	public List<ChipSeqHit> loadByRegion(ChipSeqAlignment align, Region r) throws IOException {
		try {
			List<ChipSeqHit> tor = convert(client.getPairedHits(Integer.toString(align.getDBID()),
					r.getGenome().getChromID(r.getChrom()),
					true,
					r.getStart(),
					r.getEnd(),
					null,
					null), align, true);
			tor.addAll(convert(client.getPairedHits(Integer.toString(align.getDBID()),
					r.getGenome().getChromID(r.getChrom()),
					false,
					r.getStart(),
					r.getEnd(),
					null,
					null), align, false));
			return tor;
		} catch (ClientException e) {
			throw new IllegalArgumentException(e);
		}
	}

	public List<ChipSeqHit> convert(Collection<PairedHit> input, ChipSeqAlignment align, boolean isLeft) {
		Genome g = align.getGenome();
		List<ChipSeqHit> output = new ArrayList<ChipSeqHit>();
		if (isLeft) {
			for (PairedHit s : input) {
				int start = s.leftPos;
				int end = s.leftStrand ? s.leftPos + s.leftLength : s.leftPos - s.leftLength;
				output.add(new ChiaPetHit(g, g.getChromName(s.leftChrom), Math.min(start,end), Math.max(start,end),
						s.leftStrand ? '+' : '-', align, s.weight));
			}
		} else {
			for (PairedHit s : input) {
				int start = s.rightPos;
				int end = s.rightStrand ? s.rightPos + s.rightLength : s.rightPos - s.rightLength;
				output.add(new ChiaPetHit(g, g.getChromName(s.rightChrom), Math.min(start,end), Math.max(start,end),
						s.rightStrand ? '+' : '-', align, s.weight));
			}
		}
		return output;
	}

	public List<ChipSeqHit> convert(Collection<SingleHit> input, ChipSeqAlignment align) {
		Genome g = align.getGenome();
		List<ChipSeqHit> output = new ArrayList<ChipSeqHit>();
		for (SingleHit s : input) {
			int start = s.pos;
			int end = s.strand ? s.pos + s.length : s.pos - s.length;
			output.add(new ChiaPetHit(g, g.getChromName(s.chrom), Math.min(start,end), Math.max(start,end),
					s.strand ? '+' : '-', align, s.weight));
		}
		return output;
	}

}
