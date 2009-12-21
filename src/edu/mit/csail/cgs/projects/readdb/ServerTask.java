package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.net.*;
import java.util.*;
import java.util.logging.*;
import javax.security.sasl.*;
import javax.security.auth.callback.*;

/** 
 * ServerTask represents a client connection.  Server creates ServerTasks when it receives
 * a connection and passes them to Dispatch.  Dispatch manages a pool of WorkerThreads and
 * assigns them to ServerTasks as the tasks appear to be available.  
 */

public class ServerTask {

    /* instance variables */
    /* Server that we're working for.  Need a reference to it so we can ask it
       for paths and configuration information and such
    */
    private Server server;
    /* instance variables set for each request */
    /* true if the client sent a close command.  If the client sent a close command or
       the thread has detected an error from the client, then shouldClose is true
       and the server will close the connection.
    */
    private boolean shouldClose;
    /* Socket, streams from the socket */
    private Socket socket;
    private BufferedInputStream instream;
    private OutputStream outstream;
    private WritableByteChannel outchannel;
    /* if authenticate was successful, this holds a username.  Null otherwise */
    private String username;
    /* checked to see what byte order to use for raw ints going across the network.
    */
    private ByteOrder myorder, clientorder;
    /* buffer for readLine */
    private int bufferpos;
    private byte[] buffer;
    private static final int MAXPARAMLINES = 100;
    /* other variables maintained across calls to Run but reset between connections */
    private String type;
    private List<String> params;
    private SaslServer sasl;
    private Map<String,String> saslprops;
    private String uname; // temporary, used by authenticate

    public ServerTask(Server serv, Socket s) throws IOException {
        myorder = ByteOrder.nativeOrder();          
        buffer = new byte[8192];
        params = new ArrayList<String>();
        saslprops = new HashMap<String,String>();
        saslprops.put("Sasl.POLICY_NOPLAINTEXT","true");
        saslprops.put("Sasl.POLICY_NOANONYMOUS","true");
        server = serv;
        socket = s;
        shouldClose = false;
        type = null;
        username = null;
        uname = null;
        instream = new BufferedInputStream(socket.getInputStream());
        outstream = socket.getOutputStream();
        outchannel = Channels.newChannel(outstream);
        bufferpos = 0;
        clientorder = myorder;
        sasl = null;
        socket.setTcpNoDelay(true);
//         if (server.debug()) {
//             System.err.println("New ServerTask " + this + " on socket " + socket);
//         }
    }
    public boolean shouldClose() {
//         if (shouldClose && server.debug()) {
//             System.err.println("Should close " + socket + " for " + this);
//         }

        return shouldClose;
    }
    public void close () {
        try {
            if (server.debug()) {
                System.err.println("Closing Socket " + socket + " for " + this);
            }
            if (!socket.isClosed()) {
                socket.close();
            }
        } catch (Exception e) {
            // ignore it
        }
    }
    public boolean inputAvailable() {
        boolean avail = false;
        try {
            if (!socket.isClosed() && socket.isConnected()) {
                avail = instream.available() > 0;
            } else {
//                 if (server.debug()) {
//                     System.err.println("Setting shouldClose 1 " + socket + " for " + this);
//                 }
                shouldClose = true;
            }                
        } catch (IOException e) {
//             if (server.debug()) {
//                 System.err.println("Setting shouldClose 2 " + socket + " for " + this);
//             }
            avail = false;
            shouldClose = true;
        }
        return avail;
    }
    /**
     * main method for the task.  This method is asynchronous- it shouldn't block too long on the client.  It does block
     * on disk reads and such and it does block on the client while waiting, eg, for more hits to store.  It does
     * not block while reading parameters or waiting for the next command.
     *
     * The asynchronous behavior is achieved through readLine(), which returns null if it doesn't have a complete line.
     * Typically run then returns and waits until it's called again, at which point there is hopefully a complete
     * line to deal with.  aside from authenticate(), the basic procedure is
     *   - read the type of request
     *   - read the list of request parameters
     *   - call processRequest()
     *
     * Run CANNOT throw any exceptions in the current model.  If the outer block gets an exception, it swallows it
     * and returns, setting shouldClose = true to indicate that this connection to a client should be closed.
     */
    public void run() {
        try {
//             if (server.debug()) {
//                 System.err.println("ST running with username=" + username + " type=" + type + " params="+ params);
//             }
            if (username == null) {
                if (!authenticate()) {
                    printAuthError();
//                     if (server.debug()) {
//                         System.err.println("Setting shouldClose 3 " + socket + " for " + this);
//                     }
                    shouldClose = true;
                    return;
                }
                if (username == null) { 
                    return ;
                }
                server.getLogger().log(Level.INFO,"ServerTask " + Thread.currentThread() + " authenticated " + username + " from " + socket.getInetAddress());
                printString("authenticated as " + username + "\n");
            }
            if (type == null) {
                type = readLine();
            }
            if (type == null) { 
                return;
            } else {
                while (true) {
                    String p = readLine();
                    if (p == null) { 
                        return; 
                    } else {
                        if (p.equals("ENDREQUEST")) {
                            processRequest();
                            type = null;
                            params.clear();
                            if (outstream != null) { outstream.flush(); }                            
                            break;
                        } else {
                            params.add(p);
                            if (params.size() > MAXPARAMLINES) {
//                                 if (server.debug()) {
//                                     System.err.println("Setting shouldClose 4 " + socket + " for " + this);
//                                 }
                                shouldClose = true;
                                return;
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
//             if (server.debug()) {
//                 System.err.println("Setting shouldClose 5 " + socket + " for " + this);
//             }
            shouldClose = true;
            Lock.releaseLocks(this);
            System.gc();
            System.runFinalization();
            return;
        }
    }

    /** reads and handles a request on the Socket.
     */
    public void processRequest () throws IOException {
        if (type.length() == 0) {
//             if (server.debug()) {
//                 System.err.println("Setting shouldClose 6 " + socket + " for " + this);
//             }
            shouldClose = true;
            printString("Empty request type");
            return;
        }
        if (type.equals("exists")) {
            processExists(params);            
        } else if (type.equals("store")) {
            processStore(params);
        } else if (type.equals("bye")) {
            if (server.debug()) {
                System.err.println("Got bye");
                //                System.err.println("Setting shouldClose 7 " + socket + " for " + this);
            }
            shouldClose = true;
        } else if (type.equals("getchroms")) {
            processGetChroms(params);
        } else if (type.equals("getacl")) {
            processGetACL(params);            
        } else if (type.equals("setacl")) {
            processSetACL(params);
        } else if (type.equals("deletealign")) {
            processDeleteAlignment(params);
        } else if (type.equals("byteorder")) {
            processByteOrder(params);
        } else if (type.equals("addtogroup")) {
            processAddToGroup(params);
        } else if (type.equals("shutdown")) {
            server.getLogger().log(Level.INFO,"Received shutdown from " + username);
            if (server.isAdmin(username)) {
                printOK();
                server.keepRunning(false);
            } else {
                printAuthError();
            }
            shouldClose = true;
        } else {
            processFileRequest(type,params);
        }
        long done = System.currentTimeMillis();
        //        server.getLogger().log(Level.INFO,String.format("%s from %s took %d", type, username, done - afterType));
    }
    /**
     * Handles the subset of requests that deal with a particular
     * file that we expect to exist
     */
    public void processFileRequest(String type,List<String> params)  throws IOException{
        String alignid = params.remove(0);
        String chromid = params.remove(0);
        if (alignid == null || alignid.length() == 0) {
            printString("null or empty alignment " + alignid);
            return;
        }
        if (chromid == null || chromid.length() == 0) {
            printString("null or empty chromosome " + chromid);
            return;
        }

        File directory = new File(server.getAlignmentDir(alignid));
        if (!directory.exists()) {
            printString("No Such Alignment\n");
            return;
        } 
        String headerfname = server.getHeaderFileName(alignid, chromid);
        String aclfname = server.getACLFileName(alignid);
        Lock.readLock(this,headerfname);
        Lock.readLock(this,aclfname);
        Header header;
        AlignmentACL acl = null;        
        try {
            File hf = new File(headerfname);
            if (!hf.exists()) {
                printString("No Such chromosome or alignment\n");
                Lock.readUnLock(this,headerfname);
                Lock.readUnLock(this,aclfname);
                return;
            }
            acl = server.getACL(aclfname);
        } catch (IOException e) {
            // happens if the file doesn't exist or if we can't read it at the OS level
            server.getLogger().log(Level.INFO,String.format("read error on %s, %s : %s",
                                                            aclfname,
                                                            server.getACLFileName(alignid),
                                                            e.toString()));
            printAuthError();
            Lock.readUnLock(this,aclfname);
            Lock.readUnLock(this,headerfname);
            return;
        }
        try {
            header = server.getHeader(headerfname);
        } catch (IOException e) {
            // happens if the file doesn't exist or if we can't read it at the OS level
            server.getLogger().log(Level.INFO,String.format("read error on %s, %s : %s",
                                                            headerfname,
                                                            server.getHeaderFileName(alignid, chromid),
                                                            e.toString()));
            printAuthError();
            Lock.readUnLock(this,aclfname);
            Lock.readUnLock(this,headerfname);
            return;
        }

        if (!authorizeRead(acl)) {
            printAuthError();
            server.getLogger().log(Level.INFO,String.format("%s can't read %s, %s",
                                                            username,
                                                            aclfname,
                                                            server.getACLFileName(alignid)));
            Lock.readUnLock(this,aclfname);
            Lock.readUnLock(this,headerfname);
            return;
        }
        Lock.readUnLock(this,aclfname);
        String hitsfname = server.getHitsFileName(alignid, chromid);
        String weightsfname = server.getWeightsFileName(alignid, chromid);
        Lock.readLock(this,hitsfname);
        Lock.readLock(this,weightsfname);
        Hits hits = server.getHits(hitsfname, weightsfname);
        if (type.equals("count")) {
            processCount(header,hits,params);
        } else if (type.equals("weight")) {
            processWeight(header,hits,params);
        } else if (type.equals("countrange")) {
            processCountRange(header,hits,params);
        } else if (type.equals("weightrange")) {
            processWeightRange(header,hits,params);
        } else if (type.equals("gethits")) {
            processGetHits(header,hits,params);
        } else if (type.equals("gethitsrange")) {
            processGetHitsRange(header,hits,params);
        } else if (type.equals("histogram")) {
            processHistogram(header,hits,params);
        } else if (type.equals("weighthistogram")) {
            processWeightHistogram(header,hits,params);
        } else if (type.equals("getweights")) {
            processGetWeights(header,hits,params);
        } else if (type.equals("getweightsrange")) {
            processGetWeightsRange(header,hits,params);
        } else {
            printInvalid("request type");
        }
        hits = null;
        header = null;
        Lock.readUnLock(this,weightsfname);
        Lock.readUnLock(this,hitsfname);
        Lock.readUnLock(this,headerfname);
     }

    /**
     * performs authentication exchange over the socket and sets the username field
     * if successful.  Returns true if authenticate should continue or is successful.
     * Returns false if authenticate has failed.
    */
    public boolean authenticate() throws IOException {
        /* The SASL client and server give you back bytes to send
           to the other side.  We achieve this by sending a length
           line first (ascii encoded integer followed by '\n')
           and then the raw bytes of the SASL exchange.  Two complexities:
           1) input.read() doesn't necessarily read the expected number
              of bytes all at once, so we have to loop around it until
              it does.
           2) I had problems with isComplete() returning true at different
              times in the client and server.  The server sends an isComplete() byte
              at the end of the loop to tell the client when it's done.
        */
        if (uname == null) {
            uname = readLine();
            if (uname == null) {
                return true;
            }
        }
        if (sasl == null) {
            sasl = Sasl.createSaslServer(Server.SaslMechanisms[0],
                                         "readdb",
                                         socket.getInetAddress().getCanonicalHostName(),
                                         saslprops,
                                         new ServerTaskCallbackHandler(server,uname));
        }
        if (sasl == null || sasl.isComplete()) {
            outstream.write("0\n".getBytes());
            outstream.write((byte)0);
            outstream.flush();
            server.getLogger().log(Level.INFO,"Failed Authentication for " + uname);
            return false;
        }
        while (!sasl.isComplete()) {
            try {
                String l = readLine();
                if (l == null) {
                    return true;
                }
                int length = Integer.parseInt(l);
                byte[] response = new byte[length];
                int read = 0;
               while (read < length) {
                    read += instream.read(response, read, length - read);
                    //                    System.err.println("   read " + read);
                }                
                byte[] challenge = sasl.evaluateResponse(response);
                if (challenge == null) {
                    challenge = new byte[0];
                }
                //                System.err.println("Got challenge of length " +  challenge.length + " in response to length " + read);
                String s = challenge.length + "\n";
                outstream.write(s.getBytes());
                outstream.write(challenge);
                outstream.write(sasl.isComplete() ? (byte)0 : (byte)1);
                outstream.flush();
                //                System.err.println("authenticate complete ? " + sasl.isComplete());
            } catch (Exception e) {
                System.err.println("CAUGHT AN EXCEPTION IN AUTHENTICATE");
                e.printStackTrace();
                outstream.write("0\n".getBytes());
                outstream.write((byte)0);
                outstream.flush();
                break;
            }
        }
        if (sasl.isComplete()) {
            username = sasl.getAuthorizationID();
            sasl.dispose();
            return true;
        } else {
            sasl.dispose();
            server.getLogger().log(Level.INFO,"Failed Authentication for " + uname);
            return false;
        }

    }
    /* returns true iff the user named in the username field is allowed to
       access this file
    */
    public boolean authorizeRead(AlignmentACL acl) {
        return authorize(acl.getReadACL());
    }
    public boolean authorizeWrite(AlignmentACL acl) {
        return authorize(acl.getWriteACL());
    }
    public boolean authorizeAdmin(AlignmentACL acl) {
        return authorize(acl.getAdminACL());
    }
    private boolean authorize(Set<String> acl) {
        if (acl.contains(username)) {
            return true;
        }
        for (String g : acl) {
            if (server.groupContains(username, g)) {
                return true;
            }
        }
        return false;
    }
    /** prints the response header signifying a valid request.  Only happens after
     *  the ServerTask has read enough information from the socket and done
     *  whatever else needs doing to be sure that it can satisfy the request.
     */
    public void printOK() throws IOException {
        printString("OK\n");
    }
    /** prints the response header signifying an invalid request
     */
    public void printInvalid(String reason) throws IOException {
        printString("INVALID " + reason + "\n");
    }
    /** prints the response header signifying lack of permissions */
    public void printAuthError() throws IOException {
        printString("Permission Denied\n");
    }
    /** sends the string s to the client
     */
    public void printString(String s) throws IOException {
//         if (server.debug()) {
//             System.err.println("SEND " + s);
//         }
        outstream.write(s.getBytes());
        outstream.flush();
    }
    public void processAddToGroup(List<String> params) throws IOException {
        String princ = params.get(0);
        String group = params.get(0);
        if (!server.isAdmin(username)) {
            printAuthError();
            return;
        }
        server.addToGroup(this,group,princ);
        printOK();
    }
    /**
     * Reads a line from the socket and returns it
     */
    public String readLine() throws IOException {
        int i = 0;
        boolean done = false;
        while (instream.available() > 0 && 
               (i = instream.read()) != -1) {
            if (i == '\n') {
                done = true;
                break;
            } else {
                buffer[bufferpos++] = (byte)i;
            }
        }
        if (i == -1) {
            if (server.debug()) {
                System.err.println("Setting shouldClose 9 " + socket + " for " + this);
            }
            shouldClose = true;
            return null;
        }

        if (done) {
            String out = new String(buffer,0,bufferpos);
       //     lastbufferpos = 0;
            bufferpos = 0;
//             if (server.debug()) {
//                 System.err.println("READ " +out);
//             }            
            return out;
        } else {
            return null;
        }
    }
    /**
     * input:
     *   alignid
     * output:
     *
     */
    public void processGetACL(List<String> params) throws IOException {
        String alignid = params.get(0);
        String aclFname = server.getACLFileName(alignid);
        Lock.readLock(this,aclFname);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(aclFname);
        } catch (IOException e) {
            printString("No such alignment");
            Lock.readUnLock(this,aclFname);
            return;
        }
        Lock.readUnLock(this,aclFname);
        if (!authorizeAdmin(acl)) {
            printAuthError();
            return;
        }
        printOK();
        StringBuffer sb = new StringBuffer();
        sb.append("READ\n");
        sb.append(acl.getReadACL().size() + "\n");
        for (String s : acl.getReadACL()) {
            sb.append(s + "\n");
        }
        sb.append("WRITE\n");
        sb.append(acl.getWriteACL().size() + "\n");
        for (String s : acl.getWriteACL()) {
            sb.append(s + "\n");
        }
        sb.append("ADMIN\n");
        sb.append(acl.getAdminACL().size() + "\n");
        for (String s : acl.getAdminACL()) {
            sb.append(s + "\n");
        }
        printString(sb.toString());
    }
    public void processSetACL(List<String> params) throws IOException {
        String alignid = params.remove(0);
        String aclFname = server.getACLFileName(alignid);
        Lock.writeLock(this,aclFname);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(aclFname);
        } catch (IOException e) {
            printString("No such alignment");
            Lock.writeUnLock(this,aclFname);
            return;
        }
        if (!authorizeAdmin(acl) && !server.isAdmin(username)) {
            printAuthError();
            Lock.writeUnLock(this,aclFname);
            return;
        }
        for (int i = 0; i < params.size(); i++) {
            String line = params.get(i);
            String pieces[] = line.split(" ");
            /* line format is principal [add|delete] [read|write|admin] */
            Set<String> aclset = pieces[2].equals("admin") ? acl.getAdminACL() : 
                (pieces[2].equals("write") ? acl.getWriteACL() : 
                 (pieces[2].equals("read") ? acl.getReadACL() : null));
            if (aclset == null) {
                printString("Bad ACL Type " + pieces[2] + "\n");
                continue;
            }
            if (pieces[1].equals("add")) {
                aclset.add(pieces[0]);
            } else if (pieces[1].equals("delete")) {
                aclset.remove(pieces[0]);
            } else {
                printString("Bad Operation Type " + pieces[1] + "\n");
                continue;
            }
        }
        acl.writeToFile(server.getACLFileName(alignid));
        server.removeACL(aclFname);
        Lock.writeUnLock(this,aclFname);
        printString("OK\n");                        
    }    
    /** reads two lines from socket: alignment id and chromosome id.
     * returns "exists" or unknown" to indicate whether the 
     * server knows about that pair
    */
    public void processExists(List<String> params) throws IOException {
        String alignid = params.get(0);
        String fname = server.getACLFileName(alignid);
        Lock.readLock(this,fname);
        try {
            AlignmentACL acl = server.getACL(fname);
            if (authorizeRead(acl)) {
                printString("exists\n");
            } else {
                printString("exists but no read permissions\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            printString("unknown\n");
        }
        Lock.readUnLock(this,fname);
    }
    /**
     * Returns the list of chromosomes for an alignment
     */
    public void processGetChroms(List<String> params) throws IOException {
        String alignid = params.get(0);
        String aclfname = server.getACLFileName(alignid);
        String dirname = server.getAlignmentDir(alignid);
        Lock.readLock(this,aclfname);
        Lock.readLock(this,dirname);
        File directory = new File(dirname);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(aclfname);
        } catch (IOException e) {
            printString("No Such Alignment");
            Lock.readUnLock(this,dirname);
            Lock.readUnLock(this,aclfname);
            return;
        }

        if (!authorizeRead(acl)) {
            printAuthError();
            Lock.readUnLock(this,dirname);
            Lock.readUnLock(this,aclfname);
            return;
        }        
        File[] files = directory.listFiles();
        if (files == null) {
            printString("No Such Alignment or couldn't read directory\n");
            Lock.readUnLock(this,dirname);
            Lock.readUnLock(this,aclfname);
            return;
        }

        Set<String> chroms = new HashSet<String>();
        for (int i = 0; i < files.length; i++) {
            String f = files[i].getName();
            if (f.matches("^.*\\.hits$")) {
                chroms.add(f.replaceAll("\\.hits$",""));
            }
        }           
        printOK();
        printString(chroms.size() + "\n");
        for (String f : chroms) {
            printString(f + "\n");
        }            
        Lock.readUnLock(this,dirname);
        Lock.readUnLock(this,aclfname);
    }
    /**
     * Deletes an alignment: the header and hits files, acl file, and the directory are removed
     *
     */
    public void processDeleteAlignment(List<String> params) throws IOException {
        String alignid = params.get(0);
        String dir = server.getAlignmentDir(alignid);        
        String aclFname = server.getACLFileName(alignid);

        Lock.readLock(this,dir);
        Lock.readLock(this,aclFname);

        AlignmentACL acl = server.getACL(aclFname);
        if (!authorizeAdmin(acl)) {
            printAuthError();
            Lock.readUnLock(this,dir);
            Lock.readUnLock(this,aclFname);
            return;
        }
        File directory = new File(server.getAlignmentDir(alignid));
        File[] files = directory.listFiles();
        List<String> toDelete = new ArrayList<String>();
        /* list of files to delete:
           datafiles first, then ACL, then the directory itself
        */
        for (int i = 0; i < files.length; i++) {
            String f = files[i].getName();
            if (f.matches("^.*\\.index$")) {
                String chrom = f.replaceAll("\\.index$","");
                toDelete.add(server.getHeaderFileName(alignid,chrom));
                toDelete.add(server.getHitsFileName(alignid,chrom));
                toDelete.add(server.getWeightsFileName(alignid,chrom));
            }
        }       
        toDelete.add(aclFname);
        toDelete.add(dir);
        boolean allDeleted = true;
        File f;
        for (String fname : toDelete) {
            // lock
            Lock.writeLock(this,fname);
            // remove from cache; could skip this and wait for them to fall out,
            // but this frees up entries sooner
            if (fname.matches(".*index")) {
                server.removeHeader(fname);
            } else if (fname.matches(".*hits")) {
                server.removeHits(fname);
            } else if (fname.matches(".weights")) {
                // skip
            } else if (fname.matches(".txt")) {
                server.removeACL(fname);
            }
            // file system delete
            f = new File(fname);
            allDeleted = allDeleted && f.delete();
            // remove locks.  necessary in case the file is recreated later
            Lock.writeUnLock(this,fname);
        }
        if (allDeleted) {
            printOK();
        } else {
            printString("Partially Deleted\n");
        }
        Lock.readUnLock(this,dir);
        Lock.readUnLock(this,aclFname);
    }
    /** creates or appends to a set of hits.  Reads the alignment, chromosome, and 
     * number of hits that will be sent.  Then expects to read one hit position and one weight per line 
     * (string representation of the numbers) separated by a space
     *
     * If the chromosome file doesn't exist yet, then create a new one and dump in positions and weights.
     * If it does exist, then create a new file and merge the old file with the new
     * set of hits.
     */
    public void processStore(List<String> params) throws IOException {
        String alignid = params.get(0);
        String chromid = params.get(1);
        int numHits = Integer.parseInt(params.get(2));
        printOK();
        if (numHits == 0) {
            printOK();
            return;
        }
        String dir = server.getAlignmentDir(alignid);        
        Lock.writeLock(this,dir);
        File positionsfile, weightsfile, headerfile;
        RandomAccessFile positionsRAF, weightsRAF;
        FileChannel positionsFC, weightsFC;
        ByteBuffer positionsBB, weightsBB;
        headerfile = new File(server.getHeaderFileName(alignid,chromid));
        positionsfile = new File(server.getHitsFileName(alignid,chromid));
        weightsfile = new File(server.getWeightsFileName(alignid,chromid));
        Lock.writeLock(this,headerfile.getPath());
        Lock.writeLock(this,positionsfile.getPath());
        Lock.writeLock(this,weightsfile.getPath());
        try {
            positionsfile.getParentFile().mkdirs();
            weightsfile.getParentFile().mkdirs();
            headerfile.getParentFile().mkdirs();
            if (headerfile.exists()) {
                Header header;
                try {
                    header = server.getHeader(headerfile.getCanonicalPath());
                    AlignmentACL acl = server.getACL(server.getACLFileName(alignid));
                    if (!authorizeRead(acl) || !authorizeWrite(acl)) {
                        printAuthError();

                        Lock.writeUnLock(this,weightsfile.getPath());
                        Lock.writeUnLock(this,positionsfile.getPath());
                        Lock.writeUnLock(this,headerfile.getPath());
                        Lock.writeUnLock(this,dir);
                        return;
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    printAuthError();

                    Lock.writeUnLock(this,weightsfile.getPath());
                    Lock.writeUnLock(this,positionsfile.getPath());
                    Lock.writeUnLock(this,headerfile.getPath());
                    Lock.writeUnLock(this,dir);
                    return;
                }
            } else {
                if (!server.canCreate(username)) {
                    printAuthError();

                    Lock.writeUnLock(this,weightsfile.getPath());
                    Lock.writeUnLock(this,positionsfile.getPath());
                    Lock.writeUnLock(this,headerfile.getPath());
                    Lock.writeUnLock(this,dir);
                    return;
                }
            }
            IntBP positions = new IntBP(numHits);
            FloatBP weights = new FloatBP(numHits);
            ReadableByteChannel rbc = Channels.newChannel(instream);
            Bits.readBytes(positions.bb, rbc);
            Bits.readBytes(weights.bb, rbc);
            if (myorder != clientorder) {
                Bits.flipByteOrder(positions.ib);
                Bits.flipByteOrder(weights.fb);
            }
            long[] hits = new long[numHits];
            for (int j = 0; j < numHits; j++) {
                long val = positions.get(j);
                val = val << 32;
                val += Float.floatToRawIntBits(weights.get(j));
                hits[j] = val;
            }
            positions = null;
            weights = null;
            Arrays.sort(hits);               
            positionsRAF = new RandomAccessFile(positionsfile, "rw");
            weightsRAF = new RandomAccessFile(weightsfile, "rw");
            positionsFC = positionsRAF.getChannel();
            weightsFC = weightsRAF.getChannel();
            if (headerfile.exists() && positionsFC.size() != 0) {
                positionsBB = positionsFC.map(FileChannel.MapMode.READ_WRITE,
                                              0,
                                              positionsFC.size());
                positionsBB.order(myorder);
                weightsBB = weightsFC.map(FileChannel.MapMode.READ_WRITE,
                                          0,
                                          weightsFC.size());
                weightsBB.order(myorder);
                positions = new IntBP(positionsBB);
                weights = new FloatBP(weightsBB);            
                long newsize = numHits * 4 + positionsRAF.getChannel().size();
                File newpositionsfile = new File(server.getHitsFileName(alignid,chromid) + ".tmp");
                File newweightsfile = new File(server.getWeightsFileName(alignid,chromid) + ".tmp");
                newpositionsfile.getParentFile().mkdirs();
                newweightsfile.getParentFile().mkdirs();            
                RandomAccessFile newpositionsRAF, newweightsRAF;
                FileChannel newpositionsFC, newweightsFC;
                newpositionsRAF = new RandomAccessFile(newpositionsfile, "rw");
                newweightsRAF = new RandomAccessFile(newweightsfile, "rw");
                newpositionsFC = newpositionsRAF.getChannel();
                newweightsFC = newweightsRAF.getChannel();
                ByteBuffer newpositionsBB = newpositionsFC.map(FileChannel.MapMode.READ_WRITE,
                                                               0,
                                                               newsize);
                newpositionsBB.order(myorder);
                ByteBuffer newweightsBB = newweightsFC.map(FileChannel.MapMode.READ_WRITE,
                                                           0,
                                                           newsize);
                newweightsBB.order(myorder);
                IntBP newpositions = new IntBP(newpositionsBB);
                FloatBP newweights = new FloatBP(newweightsBB);           
                // now merge the old disk files and the new hits in the array
                // into the new disk file.
                // i is index in hits
                // j is index in positions
                // k is index in newpositions
                int i = 0, j = 0, k = 0;
                int newp = (int)(hits[i] >> 32);
                float neww = Float.intBitsToFloat((int)(hits[i] & 0xFFFFFFFFL));
                int oldp = positions.get(j);
                float oldw = weights.get(j);
                while (i < hits.length && j < positions.size()) {
                    if (newp < oldp || (oldp == newp && neww < oldw)) {
                        newpositions.put(k, newp);
                        newweights.put(k++,neww);
                        i++;
                        if (i < hits.length) {
                            newp = (int)(hits[i] >> 32);
                            neww = Float.intBitsToFloat((int)(hits[i] & 0xFFFFFFFFL));
                        }
                    } else {
                        newpositions.put(k, oldp);
                        newweights.put(k++,oldw);
                        j++;
                        if (j < positions.size()) {
                            oldp = positions.get(j);
                            oldw = weights.get(j);
                        }
                    }               
                }
                while (i < hits.length) {
                    newpositions.put(k, newp);
                    newweights.put(k++,neww);
                    i++;
                    if (i < hits.length) {
                        newp = (int)(hits[i] >> 32);
                        neww = Float.intBitsToFloat((int)(hits[i] & 0xFFFFFFFFL));
                    }
                }
                while (j < positions.size()) {
                    newpositions.put(k, oldp);
                    newweights.put(k++,oldw);
                    j++;
                    if (j < positions.size()) {
                        oldp = positions.get(j);
                        oldw = weights.get(j);
                    }
                }

                // close the old files
                weights = null;
                positions = null;
                weightsBB = null;
                positionsBB = null;
                weightsFC.close();
                positionsFC.close();
                weightsRAF.close();
                positionsRAF.close();
                server.removeHits(positionsfile.getPath());
                server.removeHeader(headerfile.getPath());
                // swap in new data structures
                weights = newweights;
                positions = newpositions;
                weightsBB = newweightsBB;
                positionsBB = newpositionsBB;
                weightsFC = newweightsFC;
                positionsFC = newpositionsFC;
                weightsRAF = newweightsRAF;
                positionsRAF = newpositionsRAF;
                // rename the files into the right spot
                if (!(newpositionsfile.renameTo(positionsfile) &&
                      newweightsfile.renameTo(weightsfile))) {
                    throw new IOException("Couldn't put new files into place");
                }            
            } else {
                positionsBB = positionsFC.map(FileChannel.MapMode.READ_WRITE,
                                              0,
                                              hits.length * 4);
                positionsBB.order(myorder);
                weightsBB = weightsFC.map(FileChannel.MapMode.READ_WRITE,
                                          0,
                                          hits.length * 4);
                weightsBB.order(myorder);
                positions = new IntBP(positionsBB);
                weights = new FloatBP(weightsBB);
                for (int i = 0; i < hits.length; i++) {
                    float w = Float.intBitsToFloat((int)(hits[i] & 0xFFFFFFFFL));
                    positions.put(i, (int)(hits[i] >> 32));
                    weights.put(i, w);
                }            
                AlignmentACL acl = new AlignmentACL();
                try {
                    acl.readFromFile(server.getDefaultACLFileName());
                } catch (IOException e) {
                    // no default acl, so dont' worry.
                }

                acl.getAdminACL().add(username);
                acl.getWriteACL().add(username);
                acl.getReadACL().add(username);
                acl.writeToFile(server.getACLFileName(alignid));        
                server.removeACL(server.getACLFileName(alignid)); // make sure the server doesn't have this ACL cached
            }
            Header header = new Header(positions.ib);
            header.writeIndexFile(headerfile.getCanonicalPath());
            header.close();
            weights = null;
            positions = null;
            weightsBB = null;
            positionsBB = null;
            weightsFC.close();
            positionsFC.close();
            weightsRAF.close();
            positionsRAF.close();

            // remove from cache so they'll be re-read
            server.removeHits(positionsfile.getPath());
            server.removeHeader(headerfile.getPath());

            printOK();
        } catch (Exception e) {
            printString("Failed to store hits");
            e.printStackTrace();
        } finally {
            Lock.writeUnLock(this,weightsfile.getPath());
            Lock.writeUnLock(this,positionsfile.getPath());
            Lock.writeUnLock(this,headerfile.getPath());
            Lock.writeUnLock(this,dir);
            Lock.writeUnLock(this,dir);
        }
    }
    /* reads two lines from socket: alignment id and chromosome id.
       returns number of hits for the specified parameters
    */
    public void processCount(Header header, Hits hits, List<String> params) throws IOException {
        printOK();
        printString(Integer.toString(header.getNumHits()) + "\n");
    }
    /* reads two lines from socket: alignment id and chromosome id.
       returns total weight
    */
    public void processWeight(Header header, Hits hits, List<String> params) throws IOException {
        printOK();
        double total = 0;
        FloatBP f = hits.getWeightsBuffer();
        for (int i = 0; i < f.limit(); i++) {
            total += f.get(i);
        }
        printString(Double.toString(total) + "\n");
    }
    /* returns number ofhits in a range
     */
    public void processCountRange(Header header, Hits hits, List<String> params) throws IOException {        
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        String weightS = params.get(2);
        printOK();
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        int count;
        if (weightS.equals("NaN")) {
            count = hits.getCountBetween(first,last,startpos,stoppos);
        } else {
            count = hits.getCountBetween(first,last,startpos,stoppos,Float.parseFloat(weightS));
        }
        printString(Integer.toString(count) + "\n");        
    }
    /** gets sum of weight in a specified range
     */
    public void processWeightRange(Header header, Hits hits, List<String> params) throws IOException {
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        String weightS = params.get(2);
        printOK();
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        double total = 0;
        FloatBP f = hits.getWeightsBetween(first,last,startpos,stoppos);
        if (weightS.equals("NaN")) {
            for (int i = 0; i < f.limit(); i++) {
                total += f.get(i);
            }
        } else {
            float minweight = Float.parseFloat(weightS);
            for (int i = 0; i < f.limit(); i++) {
                if (f.get(i) > minweight) {
                    total += f.get(i);
                }
            }
        }
        printString(Double.toString(total) + "\n");
    }
    /* returns all hits for alignment/chromosome */
    public void processGetHits(Header header, Hits hits, List<String> params) throws IOException {
        IntBP h = hits.getPositionsBuffer();
        printOK();
        printString(Integer.toString(h.size()) + "\n");
        Bits.sendBytes(h.bb, outchannel);
    }
    /* returns hits in a range for an alignment/chromosome
     */
    public void processGetHitsRange(Header header, Hits hits, List<String> params) throws IOException {
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        String weightS = params.get(2);
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        IntBP h;
        if (weightS.equals("NaN")) {
            h = hits.getHitsBetween(first,last,startpos,stoppos);
        } else {
            h = hits.getHitsBetween(first,last,startpos,stoppos,Float.parseFloat(weightS));
        }
        printOK();
        printString(Integer.toString(h.size()) + "\n");
        Bits.sendBytes(h.bb, outchannel);
    }
    public void processGetWeights(Header header, Hits hits, List<String> params) throws IOException {
        FloatBP w = hits.getWeightsBuffer();
        printOK();
        printString(Integer.toString(w.size()) + "\n");
        Bits.sendBytes(w.bb, outchannel);
    }
    public void processGetWeightsRange(Header header, Hits hits, List<String> params) throws IOException {
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        String weightS = params.get(2);
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        FloatBP w;
        if (weightS.equals("NaN")) {
            w = hits.getWeightsBetween(first,last,startpos,stoppos);
        } else {
            w = hits.getWeightsBetween(first,last,startpos,stoppos,Float.parseFloat(weightS));
        }

        printOK();
        printString(Integer.toString(w.size()) + "\n");
        Bits.sendBytes(w.bb, outchannel);
    }

    
    /* returns a histogram of hits in a region.  Inputs
     * startposition, stopposition, binsize
     *
     * output array of integers is alternating bin-center and count.
     * bins with zero count are not included.
     */
    public void processHistogram(Header header, Hits hits, List<String> params) throws IOException {
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        int binsize = Integer.parseInt(params.get(2));
        String weightS = params.get(3);
        int extension = Integer.parseInt(params.get(4));
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        int[] raw = null;
        if (weightS.equals("NaN")) {
            raw = hits.histogram(first,
                                 last,
                                 startpos,
                                 stoppos,
                                 binsize,
                                 extension);
        } else {
            raw = hits.histogram(first,
                                 last,
                                 startpos,
                                 stoppos,
                                 binsize,
                                 Float.parseFloat(weightS),
                                 extension);
        }
        int n = 0;
        for (int i = 0; i< raw.length; i++) {
            if (raw[i] > 0) {
                n++;
            }
        }
        int[] hist = new int[n*2];
        int pos = 0;
        for (int i = 0; i< raw.length; i++) {
            if (raw[i] > 0) {
                hist[pos*2] = startpos + binsize * i + binsize / 2;
                hist[pos*2+1] = raw[i];
                pos++;
            }
        }
        printOK();
        printString(Integer.toString(hist.length) + "\n");
        Bits.sendInts(hist, outstream, buffer, clientorder);        
    }

    /* returns a histogram of hit weights in a region.  Inputs
     * startposition, stopposition, binsize.  Each bin's
     * value is the total amount of weight in the bin.
     *
     * bins with zero count are not included.
     */
    public void processWeightHistogram(Header header, Hits hits, List<String> params) throws IOException {
        int startpos = Integer.parseInt(params.get(0));
        int stoppos = Integer.parseInt(params.get(1));
        int binsize = Integer.parseInt(params.get(2));
        String weightS = params.get(3);
        int extension = Integer.parseInt(params.get(4));
        int first = header.getFirstIndex(startpos);
        int last = header.getLastIndex(stoppos);
        float[] raw;
        if (weightS.equals("NaN")) {
            raw = hits.weightHistogram(first,
                                       last,
                                       startpos,
                                       stoppos,
                                       binsize,
                                       extension);
        } else {
            raw = hits.weightHistogram(first,
                                       last,
                                       startpos,
                                       stoppos,
                                       binsize,
                                       Float.parseFloat(weightS),
                                       extension);
        }
        int n = 0;
        for (int i = 0; i< raw.length; i++) {
            if (raw[i] > 0) {
                n++;
            }
        }
        int[] parray = new int[n];
        float[] farray = new float[n];
        int pos = 0;
        for (int i = 0; i< raw.length; i++) {
            if (raw[i] > 0) {
                parray[pos] = startpos + binsize * i + binsize / 2;
                farray[pos] = raw[i];
                pos++;
            }
        }
        printOK();
        printString(Integer.toString(parray.length) + "\n");
        Bits.sendInts(parray, outstream, buffer, clientorder);        
        Bits.sendFloats(farray, outstream, buffer, clientorder);
    }

    public void processByteOrder(List<String> params) throws IOException {
        String orderString = params.get(0);
        if (orderString.equals("big")) {
            clientorder = ByteOrder.BIG_ENDIAN;
            printOK();
        } else if (orderString.equals("little")) {
            clientorder = ByteOrder.LITTLE_ENDIAN;
            printOK();
        } else {
            printInvalid("bad byte order " + orderString);
        }
    }
    /**
     * Reads lines up to and including "ENDREQUEST".  Returns a list
     * of all but the last line.
     */
    public List<String> readParameterLines() throws IOException {
        ArrayList<String> output = new ArrayList<String>();
        String line = readLine();
        int i = 0;
        while (line != null && line.length() > 0 && !line.equals("ENDREQUEST") && i++ < MAXPARAMLINES) {
            output.add(line);
            line = readLine();
        }
        return output;
    }
}

/** 
 * SASL callback handler for the authenticate() method.  Provides
 * a username and server name to the CallBack
 */
class ServerTaskCallbackHandler implements CallbackHandler {
    private String uname;
    private Server server;
    public ServerTaskCallbackHandler(Server s, String u) {
        uname = u;
        server = s;
    }

    public void handle(Callback[] callbacks) {
        for (int i = 0; i < callbacks.length; i++) {
            if (callbacks[i] instanceof PasswordCallback) {
                PasswordCallback pc = (PasswordCallback)callbacks[i];
                try {
                    pc.setPassword(server.getPassword(uname).toCharArray());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (callbacks[i] instanceof NameCallback) {
                NameCallback nc = (NameCallback)callbacks[i];
                nc.setName(uname);
            }            
            if (callbacks[i] instanceof AuthorizeCallback) {
                AuthorizeCallback ac = (AuthorizeCallback)callbacks[i];
                ac.setAuthorized(ac.getAuthenticationID().equals(ac.getAuthorizationID()));
            }

        }
    }
}
