package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;
import java.io.*;

public abstract class ChipChipPieceFile extends ChipChipFile {

    public ChipChipPieceFile(String filename, String type) {
        super(filename,type);
    }
    public int getPieceCount() {
        int piece = 1;
        while (true) {
            String name = getBaseFilename() + "." + piece + "." + getType();
            File f = new File(name);
            if (!f.exists()) {
                break;
            }
            piece++;
        }
        return piece - 1;
    }
    public String getPieceFilename(int piece) {
        return getBaseFilename() + "." + piece + "." + getType();
    }
    public String getPieceFilename(int piece, String type) {
        return getBaseFilename() + "." + piece + "." + type;
    }


}
