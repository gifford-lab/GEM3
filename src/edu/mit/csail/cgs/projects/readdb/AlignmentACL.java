package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import java.nio.channels.*;

public class AlignmentACL implements Closeable {

    /* set of principals that can read, write, and change acls
       for this directory.  Default ACLs are empty, so client code
       should remember to add at least one user 
    */
    private Set<String> readACL, writeACL, adminACL;

    public AlignmentACL() {
        readACL = new HashSet<String>();
        writeACL = new HashSet<String>();
        adminACL = new HashSet<String>();
    }
    public AlignmentACL(String fname) throws IOException {
        readACL = new HashSet<String>();
        writeACL = new HashSet<String>();
        adminACL = new HashSet<String>();
        readFromFile(fname);
    }
    /* get the ACLs.  These return the underlying sets to allow other
       readdb code to *add* to the acls in addition to reading them.
    */
    protected Set<String> getReadACL() {return readACL;}
    protected Set<String> getWriteACL() {return writeACL;}
    protected Set<String> getAdminACL() {return adminACL;}
    /* writes an ACL file.  locks the file and translates locking problems
       into IOExceptions
    */
    public void writeToFile(String fname) throws IOException {
        File f = new File(fname);
        FileChannel channel = new RandomAccessFile(f,"rw").getChannel();
        FileLock lock = channel.lock();

        PrintStream ps = new PrintStream(fname);
        ps.print("read:");
        for (String s : readACL) {
            ps.print(" " + s);
        }
        ps.println();
        ps.print("write:");
        for (String s : writeACL) {
            ps.print(" " + s);
        }
        ps.println();
        ps.print("admin:");
        for (String s : adminACL) {
            ps.print(" " + s);
        }
        ps.println();
        ps.close();
        
        lock.release();
        channel.close();
    }
    public void readFromFile(String fname) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = null;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split(":");
            String type = pieces[0];
            pieces = pieces[1].replaceAll("^ +","").split(" ");
                Set<String> a = type.equals("read") ? readACL : (type.equals("write") ? writeACL : (type.equals("admin") ? adminACL : null));
                if (a == null) {throw new IOException ("invalid type " + type + " in acl " + fname);}
                for (int i = 0; i < pieces.length; i++) {
                    a.add(pieces[i]);
                }
        }
        reader.close();
    }
    public void close() {}
}