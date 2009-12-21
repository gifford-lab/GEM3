package edu.mit.csail.cgs.warpdrive.components;

import edu.mit.csail.cgs.warpdrive.WarpOptions;
import java.io.File;
import java.io.IOException;
import edu.mit.csail.cgs.utils.NotFoundException;

public class Snapshot {

    
    public static void main(String args[]) {
        String picturename = null;
        int w = 1600, h = 1200;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--picture")) {
                picturename = args[++i];
            }
            if (args[i].equals("--width")) {
                w = Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--height")) {
                h = Integer.parseInt(args[++i]);
            }
        }
        if (picturename == null) {return;}
        try {
            WarpOptions opts = WarpOptions.parseCL(args);
            RegionPanel panel = new RegionPanel(opts);
            File file = new File(picturename);
            panel.computeLayout(0,0,w,h);
            while (!panel.allCanPaint()) {
                Thread.yield();
            }
            panel.saveImage(file,w,h,true);
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
            e.printStackTrace();
        }

        
        System.exit(0);
    }


}