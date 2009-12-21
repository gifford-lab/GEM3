package edu.mit.csail.cgs.utils.parsing.alignment;

import java.io.*;
import java.util.*;

public class SAMParser {

    public static final String VERSION = "$Id$";

    public static final int QNAME = 0,
        FLAGS = 1,
        REFNAME = 2,
        POS = 3,
        MAPQ = 4,
        CIGAR = 5,
        MATEREFNAME = 6,
        MATEPOS = 7,
        INSERTSIZE = 8,
        SEQ = 9,
        QUAL = 10,
        OPT = 11;

    public static final int 
        PAIRED = 0x0001,
        MAPPEDPAIR = 0x0002,
        UNMAPPEDQUERY = 0x0004,
        UNMAPPEDMATE = 0x0008,
        QUERYSTRAND = 0x0010,
        MATESTRAND = 0x0020,
        READFIRST = 0x0040,
        READSECOND = 0x0080,
        NOTPRIMARY = 0x0100,
        FAILEDQC = 0x0200,
        DUPLICATE = 0x0400;

    /**
     * Fills in the SAMRecord with the data from the SAM-formatted line.
     */
    public void parseFromLine(String line,
                              SAMRecord record,
                              boolean wantoptional) {
        String cols[] = line.split("\t");
        record.qname = cols[QNAME];
        record.flags = Integer.parseInt(cols[FLAGS]);
        record.tname = cols[REFNAME];
        record.pos = Integer.parseInt(cols[POS]);
        record.mapq = Integer.parseInt(cols[MAPQ]);
        record.cigar = cols[CIGAR];
        record.matetname = cols[MATEREFNAME];
        record.matepos = Integer.parseInt(cols[MATEPOS]);
        record.insertsize = Integer.parseInt(cols[INSERTSIZE]);
        record.seq = cols[SEQ];
        record.qual = cols[QUAL];
        record.nfields = 0;
        if (wantoptional && cols.length > OPT ) {
            record.nfields = cols.length - OPT;
            if (record.tags == null || record.tags.length < record.nfields) {
                record.tags = new String[record.nfields];
                record.values = new String[record.nfields];
                record.types = new char[record.nfields];
            }

            for (int i = OPT; i < cols.length; i++) {
                String p[] = cols[i].split(":");
                record.tags[i-OPT] = p[0];
                record.types[i-OPT] = p[1].charAt(0);
                record.values[i-OPT] = p[2];
            }
        }
    }

    public List<SAMRecord> parseStream(InputStream is, boolean wantoptional) throws IOException {
        ArrayList<SAMRecord> output = new ArrayList<SAMRecord>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        String line;
        while ((line = reader.readLine()) != null) {
            SAMRecord r = new SAMRecord();
            output.add(r);
            parseFromLine(line, r, wantoptional);
        }
        return output;
    }



}