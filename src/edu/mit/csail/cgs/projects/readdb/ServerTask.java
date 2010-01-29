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
    private static final int SINGLE = 1, PAIRED = 2;
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
    private Request request;
    private List<String> args;
    private SaslServer sasl;
    private Map<String,String> saslprops;
    private String uname; // temporary, used by authenticate

    public ServerTask(Server serv, Socket s) throws IOException {
        myorder = ByteOrder.nativeOrder();          
        buffer = new byte[8192];
        request = new Request();
        args = new ArrayList<String>();
        saslprops = new HashMap<String,String>();
        saslprops.put("Sasl.POLICY_NOPLAINTEXT","true");
        saslprops.put("Sasl.POLICY_NOANONYMOUS","true");
        server = serv;
        socket = s;
        shouldClose = false;
        username = null;
        uname = null;
        socket.setReceiveBufferSize(Server.BUFFERLEN);
        socket.setSendBufferSize(Server.BUFFERLEN);
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
            if (username == null) {
                if (!authenticate()) {
                    printAuthError();
                    shouldClose = true;
                    return;
                }
                if (username == null) { 
                    return ;
                }
                server.getLogger().log(Level.INFO,"ServerTask " + Thread.currentThread() + " authenticated " + username + " from " + socket.getInetAddress());
                printString("authenticated as " + username + "\n");
            }
            while (true) {
                String p = readLine();
                if (p == null) { 
                    return; 
                } else {
                    if (p.equals("ENDREQUEST")) {
                        String error = request.parse(args);
                        if (error == null) {
                            processRequest();
                        } else {
                            printString("error parsing request: " + error + "\n");
                            System.err.println("Error parsing request from " + args);
                        }
                        args.clear();
                        if (outstream != null) { outstream.flush(); }                            
                        break;
                    } else {
                        args.add(p);
                        if (args.size() > MAXPARAMLINES) {                            
                            shouldClose = true;
                            return;
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
//             if (server.debug()) {
//                 System.err.println("Setting shouldClose 5 " + socket + " for " + this);
//             }
            args.clear();
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
        try {
            if (request.alignid != null) {
                Lock.readLock(this,request.alignid);
            }
            if (request.type.equals("exists")) {
                processExists();            
            } else if (request.type.equals("storesingle")) {
                processSingleStore();
            } else if (request.type.equals("storepaired")) {
                processPairedStore();
            } else if (request.type.equals("bye")) {
                shouldClose = true;
            } else if (request.type.equals("getchroms")) {
                processGetChroms();
            } else if (request.type.equals("getacl")) {
                processGetACL();            
            } else if (request.type.equals("setacl")) {
                processSetACL();
            } else if (request.type.equals("deletealign")) {
                processDeleteAlignment();
            } else if (request.type.equals("byteorder")) {
                processByteOrder();
            } else if (request.type.equals("addtogroup")) {
                processAddToGroup();
            } else if (request.type.equals("shutdown")) {
                server.getLogger().log(Level.INFO,"Received shutdown from " + username);
                if (server.isAdmin(username)) {
                    printOK();
                    server.keepRunning(false);
                } else {
                    printAuthError();
                }
                shouldClose = true;
            } else {
                processFileRequest();
            }
        } catch (IOException e) {
            Lock.readUnLock(this,request.alignid);
            throw e;
        }
        Lock.readUnLock(this,request.alignid);

    }
    /**
     * Handles the subset of requests that deal with a particular
     * file that we expect to exist
     */
    public void processFileRequest()  throws IOException{
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        if (request.alignid == null || request.alignid.length() == 0) {
            printString("null or empty alignment " + request.alignid);
            return;
        }
        File directory = new File(server.getAlignmentDir(request.alignid));
        if (!directory.exists()) {
            printString("No Such Alignment\n");
            return;
        } 
        AlignmentACL acl = null;        
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            // happens if the file doesn't exist or if we can't read it at the OS level
            server.getLogger().log(Level.INFO,String.format("read error on acl for %s : %s",
                                                            request.alignid,
                                                            e.toString()));
            printAuthError();
            return;
        }
        if (!authorizeRead(acl)) {
            server.getLogger().log(Level.INFO,String.format("%s can't read %s",
                                                            username,
                                                            request.alignid));
            printAuthError();
            return;
        }
        Header header;
        Hits hits;
        try {
            if (request.isPaired) {
                hits = server.getPairedHits(request.alignid, request.chromid, request.isLeft);
                header = server.getPairedHeader(request.alignid, request.chromid, request.isLeft);
            } else {
                hits = server.getSingleHits(request.alignid, request.chromid);
                header = server.getSingleHeader(request.alignid, request.chromid);
            }
        } catch (IOException e) {
            // happens if the file doesn't exist or if we can't read it at the OS level
            server.getLogger().log(Level.INFO,String.format("read error on header or hits for %s, %d, %s : %s",
                                                            request.alignid, request.chromid, request.isLeft,
                                                            e.toString()));
            printAuthError();
            return;
        }
        if (request.type.equals("count")) {
            processCount(header,hits);
        } else if (request.type.equals("weight")) {
            processWeight(header,hits);
        } else if (request.type.equals("histogram")) {
            processHistogram(header,hits);
        } else if (request.type.equals("weighthistogram")) {
            processWeightHistogram(header,hits);
        } else if (request.type.equals("gethits")) {
            processGetHits(header,hits);
        } else {
            printInvalid("request type");
        }
        hits = null;
        header = null;
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
    public void processAddToGroup() throws IOException {
        String princ = request.map.get("princ");
        String group = request.map.get("group");
        if (!server.isAdmin(username)) {
            printAuthError();
            return;
        }
        if (princ == null || group == null) {
            printString("Must supply princ and group :" + princ + "," + group + "\n");
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
            shouldClose = true;
            return null;
        }

        if (done) {
            String out = new String(buffer,0,bufferpos);
            bufferpos = 0;
            return out;
        } else {
            return null;
        }
    }
    public void processGetACL() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No such alignment");
            return;
        }
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
    public void processSetACL() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        Lock.writeLock(this,request.alignid);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No such alignment");
            Lock.writeUnLock(this,request.alignid);
            return;
        }
        if (!authorizeAdmin(acl) && !server.isAdmin(username)) {
            printAuthError();
            Lock.writeUnLock(this,request.alignid);
            return;
        }
        for (int i = 0; i < request.list.size(); i++) {
            String line = request.list.get(i);
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
        acl.writeToFile(server.getACLFileName(request.alignid));
        server.removeACL(request.alignid);
        Lock.writeUnLock(this,request.alignid);
        printString("OK\n");                        
    }    
    /** reads two lines from socket: alignment id and chromosome id.
     * returns "exists" or unknown" to indicate whether the 
     * server knows about that pair
    */
    public void processExists() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        try {
            AlignmentACL acl = server.getACL(request.alignid);
            if (authorizeRead(acl)) {
                printString("exists\n");
            } else {
                printString("exists but no read permissions\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            printString("unknown\n");
        }
    }
    /**
     * Returns the list of chromosomes for an alignment
     */
    public void processGetChroms() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = null;
        try {
            acl = server.getACL(request.alignid);
        } catch (IOException e) {
            printString("No Such Alignment");
            return;
        }
        if (!authorizeRead(acl)) {
            printAuthError();
            return;
        }        
        Set<Integer> chroms = server.getChroms(request.alignid,
                                              request.isPaired,
                                              request.isLeft);
        if (chroms == null) {
            printString("No Such Alignment");
            return;
        }

        printOK();
        printString(chroms.size() + "\n");
        for (Integer i : chroms) {
            printString(i + "\n");
        }            
    }
    /**
     * Deletes an alignment: the header and hits files, acl file, and the directory are removed
     *
     */
    public void processDeleteAlignment() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        AlignmentACL acl = server.getACL(request.alignid);
        if (!authorizeAdmin(acl)) {
            printAuthError();
            return;
        }
        File directory = new File(server.getAlignmentDir(request.alignid));
        File[] files = directory.listFiles();
        List<String> toDelete = new ArrayList<String>();
        /* list of files to delete:
           datafiles first, then the directory itself
        */
        for (int i = 0; i < files.length; i++) {
            String name = files[i].getName();
            if (request.isPaired == null) {
                toDelete.add(name);
            } else {
                if (request.isPaired && 
                    (name.indexOf(".prleft.") > 0 ||
                     name.indexOf(".prright.") > 0)) {
                    toDelete.add(name);
                } else if (!request.isPaired && 
                           name.indexOf(".prleft") == -1 &&
                           name.indexOf(".prright") == -1) {
                    toDelete.add(name);
                }
            }
        }       
        if (request.isPaired == null) {
            toDelete.add(server.getACLFileName(request.alignid));
            toDelete.add(directory.getName());
        }
        boolean allDeleted = true;
        File f;
        Lock.writeLock(this,request.alignid);
        for (String fname : toDelete) {
            // file system delete
            f = new File(fname);
            allDeleted = allDeleted && f.delete();
        }
        if (allDeleted) {
            printOK();
        } else {
            printString("Partially Deleted\n");
        }
        Lock.writeUnLock(this,request.alignid);
    }
    /** creates or appends to a set of hits.  
     *
     * If the chromosome file doesn't exist yet, then create a new one and dump in positions and weights.
     * If it does exist, then create a new file and merge the old file with the new
     * set of hits.
     */
    public void processSingleStore() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        int numHits = 0;
        try {
            numHits = Integer.parseInt(request.map.get("numhits"));
        } catch (NumberFormatException e) {
            printString("Invalid numhits value : " + request.map.get("numhits") + "\n");
            return;
        }
        printOK();
        if (numHits == 0) {
            printOK();
            return;
        }
        Lock.writeLock(this,request.alignid);
        List<SingleHit> hits = new ArrayList<SingleHit>();

        /* if the alignment already exists, read in the old hits */
        Set<Integer> chroms = server.getChroms(request.alignid, false,false);
        if (chroms != null) {
            try {        
                /* sure we're allowed to write here */
                AlignmentACL acl = server.getACL(request.alignid);
                if (!authorizeRead(acl) || !authorizeWrite(acl)) {
                    printAuthError();
                    Lock.writeUnLock(this, request.alignid);
                    return;
                }
            } catch (Exception e) {
                e.printStackTrace();
                printAuthError();
                Lock.writeUnLock(this,request.alignid);
                return;
            }
            if (chroms.contains(request.chromid)) {
                try {
                    SingleHits oldhits = server.getSingleHits(request.alignid,
                                                              request.chromid);
                    IntBP positions = oldhits.getPositionsBuffer();
                    FloatBP weights = oldhits.getWeightsBuffer();
                    IntBP las = oldhits.getLASBuffer();
                    for (int i = 0; i < positions.limit(); i++) {
                        hits.add(new SingleHit(request.chromid,
                                               positions.get(i),
                                               weights.get(i),
                                               Hits.getStrandOne(las.get(i)),
                                               Hits.getLengthOne(las.get(i))));
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    printAuthError();
                    Lock.writeUnLock(this,request.alignid);
                    return;
                }
            }
        } else {
            /* this is a new alignment, so set a default ACL */
            AlignmentACL acl = new AlignmentACL();
            try {
                acl.readFromFile(server.getDefaultACLFileName());
            } catch (IOException e) {
                // no default acl, so dont' worry.
            }

            if (!(new File(server.getAlignmentDir(request.alignid)).mkdirs())) {
                System.err.println("Can't create directories for " + request.alignid + ":" + server.getAlignmentDir(request.alignid));
                printAuthError();
                Lock.writeUnLock(this,request.alignid);
                return;
            }
            acl.getAdminACL().add(username);
            acl.getWriteACL().add(username);
            acl.getReadACL().add(username);
            acl.writeToFile(server.getACLFileName(request.alignid));        
            server.removeACL(request.alignid); // make sure the server doesn't have this ACL cached
        }

        IntBP positions = new IntBP(numHits);
        FloatBP weights = new FloatBP(numHits);
        IntBP las = new IntBP(numHits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(positions.bb, rbc);
        Bits.readBytes(weights.bb, rbc);
        Bits.readBytes(las.bb, rbc);

        for (int i = 0; i < numHits; i++) {
            hits.add(new SingleHit(request.chromid,
                                   positions.get(i),
                                   weights.get(i),
                                   Hits.getStrandOne(las.get(i)),
                                   Hits.getLengthOne(las.get(i))));
        }
        positions = null;
        weights = null;
        las = null;
        Collections.sort(hits);
        try {
            SingleHits.writeSingleHits(hits,
                                       server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                       request.chromid);
            SingleHits singlehits = new SingleHits(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                                   request.chromid);
            Header header = new Header(singlehits.getPositionsBuffer().ib);
            header.writeIndexFile(server.getSingleHeaderFileName(request.alignid,
                                                                 request.chromid));
        } catch (IOException e) {
            printString("IOException trying to save files : " + e.toString() + "\n");
            Lock.writeUnLock(this,request.alignid);
            return;
        }
        printOK();
        server.removeSingleHits(request.alignid, request.chromid);
        server.removeSingleHeader(request.alignid, request.chromid);
        Lock.writeUnLock(this,request.alignid);
    }

    public void processPairedStore() throws IOException {
        assert(request != null);
        assert(request.alignid != null);
        assert(request.chromid != null);        
        int numHits = 0;
        try {
            numHits = Integer.parseInt(request.map.get("numhits"));
        } catch (NumberFormatException e) {
            printString("Invalid numhits value : " + request.map.get("numhits") + "\n");
            return;
        }
        try {        
            /* sure we're allowed to write here */
            AlignmentACL acl = server.getACL(request.alignid);
            if (!authorizeRead(acl) || !authorizeWrite(acl)) {
                printAuthError();
                Lock.writeUnLock(this, request.alignid);
                return;
            }
        } catch (IOException e) {
            printAuthError();
            return;
        }

        printOK();
        if (numHits == 0) {
            printOK();
            return;
        }
        Lock.writeLock(this,request.alignid);
        List<PairedHit> hits = new ArrayList<PairedHit>();

        /* if the alignment already exists, read in the old hits */
        Set<Integer> chroms = server.getChroms(request.alignid, true,request.isLeft);
        if (chroms == null) {
            /* this is a new alignment, so set a default ACL */
            AlignmentACL acl = new AlignmentACL();
            try {
                acl.readFromFile(server.getDefaultACLFileName());
            } catch (IOException e) {
                // no default acl, so dont' worry.
            }
            if (!(new File(server.getAlignmentDir(request.alignid))).mkdirs()) {
                System.err.println("Can't create directories for " + request.alignid + ":" + server.getAlignmentDir(request.alignid));
                printAuthError();
                Lock.writeUnLock(this,request.alignid);
                return;
            }
            acl.getAdminACL().add(username);
            acl.getWriteACL().add(username);
            acl.getReadACL().add(username);
            acl.writeToFile(server.getACLFileName(request.alignid));        
            server.removeACL(request.alignid); // make sure the server doesn't have this ACL cached
        }

        List<PairedHit> newhits = new ArrayList<PairedHit>();
        IntBP positions = new IntBP(numHits);
        FloatBP weights = new FloatBP(numHits);
        IntBP las = new IntBP(numHits);
        IntBP otherchrom = new IntBP(numHits);
        IntBP otherpos = new IntBP(numHits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(positions.bb, rbc);
        Bits.readBytes(weights.bb, rbc);
        Bits.readBytes(las.bb, rbc);
        Bits.readBytes(otherchrom.bb,rbc);
        Bits.readBytes(otherpos.bb,rbc);
        for (int i = 0; i < positions.limit(); i++) {
            newhits.add(new PairedHit(request.chromid,
                                      positions.get(i),
                                      Hits.getStrandOne(las.get(i)),
                                      Hits.getLengthOne(las.get(i)),
                                      otherchrom.get(i),
                                      otherpos.get(i),
                                      Hits.getStrandTwo(las.get(i)),
                                      Hits.getLengthTwo(las.get(i)),
                                      weights.get(i)));
        }
        positions = null;
        weights = null;
        las = null;
        otherchrom = null;
        otherpos = null;

        try {
            appendPairedHits(newhits,true);
            for (PairedHit h : newhits) {
                h.flipSides();            
            }
            appendPairedHits(newhits,false);
            PairedHits pairedhits = new PairedHits(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                                   request.chromid, 
                                                   true);
            Header header = new Header(pairedhits.getPositionsBuffer().ib);
            header.writeIndexFile(server.getPairedHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 true));
            pairedhits = new PairedHits(server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"),
                                        request.chromid, 
                                        false);
            header = new Header(pairedhits.getPositionsBuffer().ib);
            header.writeIndexFile(server.getPairedHeaderFileName(request.alignid,
                                                                 request.chromid,
                                                                 false));
        } catch (IOException e) {
            Lock.writeUnLock(this,request.alignid);
            printString("Failed to write hits : " + e.toString());
        }
        server.removePairedHits(request.alignid, request.chromid, true);
        server.removePairedHits(request.alignid, request.chromid, false);
        server.removePairedHeader(request.alignid, request.chromid, true);
        server.removePairedHeader(request.alignid, request.chromid, false);
        printOK();
        Lock.writeUnLock(this,request.alignid);
    }

    private void appendPairedHits(List<PairedHit> newhits,
                                  boolean isLeft) throws IOException {
        Map<Integer,List<PairedHit>> map = new HashMap<Integer,List<PairedHit>>();
        for (PairedHit h : newhits) {
            if (!map.containsKey(h.leftChrom)) {
                map.put(h.leftChrom, new ArrayList<PairedHit>());
            }
            map.get(h.leftChrom).add(h);
        }
        for (int chromid : map.keySet()) {
            List<PairedHit> hits = map.get(chromid);
            try {
                PairedHits oldhits = server.getPairedHits(request.alignid,
                                                          chromid,
                                                          isLeft);
                IntBP positions = oldhits.getPositionsBuffer();
                FloatBP weights = oldhits.getWeightsBuffer();
                IntBP las = oldhits.getLASBuffer();
                IntBP otherchrom = oldhits.getChromsBuffer();
                IntBP otherpos = oldhits.getOtherPosBuffer();
                for (int i = 0; i < positions.limit(); i++) {
                    hits.add(new PairedHit(request.chromid,
                                           positions.get(i),
                                           Hits.getStrandOne(las.get(i)),
                                           Hits.getLengthOne(las.get(i)),
                                           otherchrom.get(i),
                                           otherpos.get(i),
                                           Hits.getStrandTwo(las.get(i)),
                                           Hits.getLengthTwo(las.get(i)),
                                           weights.get(i)));
                }
            } catch (FileNotFoundException e) {
                // this is OK, it just means there were no old hits.  Any other
                // IOException is a problem, so let it propagate out
            } 
            PairedHits.writePairedHits(hits, server.getAlignmentDir(request.alignid) + System.getProperty("file.separator"), chromid, isLeft);
        }
    }

    public void processCount(Header header, Hits hits) throws IOException {
        printOK();
        if (request.start == null && request.end == null && request.minWeight == null && request.isPlusStrand == null) {
            printString(Integer.toString(header.getNumHits()) + "\n");
            return;
        }
        if (request.start == null) {
            request.start = 0;
        }
        if (request.end == null) {
            request.end = Integer.MAX_VALUE;
        }
        int count = 0;
        int first = header.getFirstIndex(request.start == null ? 0 : request.start);
        int last = header.getLastIndex(request.end == null ? Integer.MAX_VALUE : request.end);
        printString(Integer.toString(hits.getCountBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand)) + "\n");
    }
    public void processWeight(Header header, Hits hits) throws IOException {
        printOK();
        if (request.start == null) {
            request.start = 0;
        }
        if (request.end == null) {
            request.end = Integer.MAX_VALUE;
        }
        int first = header.getFirstIndex(request.start);
        int last = header.getLastIndex(request.end);
        printString(Double.toString(hits.getWeightBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand)) + "\n");
    }
    public void processGetHits(Header header, Hits hits) throws IOException {
        int count;
        if (!(hits instanceof PairedHits)) {
            if (request.map.containsKey("wantotherchroms") ||
                request.map.containsKey("wantotherpositions")) {
                printString("invalid columns requested for single-ended data\n");
                return;
            }
        }
        if (request.start == null) {
            request.start = 0;
        }
        if (request.end == null) {
            request.end = Integer.MAX_VALUE;
        }
        int first = header.getFirstIndex(request.start);
        int last = header.getLastIndex(request.end);
        if (request.start == 0 && request.end == Integer.MAX_VALUE && request.minWeight == null && request.isPlusStrand == null) {
            count = header.getNumHits();
        } else {
            count = hits.getCountBetween(first,last,request.start,request.end,request.minWeight, request.isPlusStrand);
        }
        printOK();
        printString(Integer.toString(count) + "\n");
        if (request.map.containsKey("wantpositions")) {
            IntBP p = hits.getHitsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
            Bits.sendBytes(p.bb, outchannel);
        }
        if (request.map.containsKey("wantweights")) {
            FloatBP p = hits.getWeightsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
            Bits.sendBytes(p.bb, outchannel);
        }
        if (request.map.containsKey("wantlengthsandstrands")) {
            IntBP p = hits.getLASBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
            Bits.sendBytes(p.bb, outchannel);
        }
        if (request.map.containsKey("wantotherchroms")) {
            IntBP p = ((PairedHits)hits).getOtherChromsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
            Bits.sendBytes(p.bb, outchannel);
        }
        if (request.map.containsKey("wantotherpositions")) {
            IntBP p = ((PairedHits)hits).getOtherPositionsBetween(first,last,request.start,request.end,request.minWeight,request.isPlusStrand);
            Bits.sendBytes(p.bb, outchannel);
        }        
    }
    public void processHistogram(Header header, Hits hits) throws IOException {
        int binsize = 10;
        if (request.start == null) {
            IntBP ib = hits.getPositionsBuffer();
            request.start = ib.get(0);
        }
        if (request.end == null) {
            IntBP ib = hits.getPositionsBuffer();
            request.end = ib.get(ib.limit()-1);
        }
        try {
            binsize = Integer.parseInt(request.map.get("binsize"));
        } catch (Exception e) {
            printString("missing or invalid bin size : " + request.map.get("binsize") + "\n");
            return;
        }
        boolean extension = request.map.containsKey("extension");
        int first = header.getFirstIndex(request.start);
        int last = header.getLastIndex(request.end);
        int[] raw = hits.histogram(first,
                                   last,
                                   request.start,
                                   request.end,
                                   binsize,
                                   request.minWeight,
                                   request.isPlusStrand,
                                   extension);
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
                hist[pos*2] = request.start + binsize * i + binsize / 2;
                hist[pos*2+1] = raw[i];
                pos++;
            }
        }
        printOK();
        printString(Integer.toString(hist.length) + "\n");
        Bits.sendInts(hist, outstream, buffer, myorder);        
    }

    /* returns a histogram of hit weights in a region.  Inputs
     * startposition, stopposition, binsize.  Each bin's
     * value is the total amount of weight in the bin.
     *
     * bins with zero count are not included.
     */
    public void processWeightHistogram(Header header, Hits hits) throws IOException {
        int binsize = 10;
        if (request.start == null) {
            IntBP ib = hits.getPositionsBuffer();
            request.start = ib.get(0);
        }
        if (request.end == null) {
            IntBP ib = hits.getPositionsBuffer();
            request.end = ib.get(ib.limit()-1);
        }
        try {
            binsize = Integer.parseInt(request.map.get("binsize"));
        } catch (Exception e) {
            printString("missing or invalid bin size : " + request.map.get("binsize") + "\n");
            return;
        }
        boolean extension = request.map.containsKey("extension");
        int first = header.getFirstIndex(request.start);
        int last = header.getLastIndex(request.end);
        float[] raw = hits.weightHistogram(first,
                                           last,
                                           request.start,
                                           request.end,
                                           binsize,
                                           request.minWeight,
                                           request.isPlusStrand,
                                           extension);
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
                parray[pos] = request.start + binsize * i + binsize / 2;
                farray[pos] = raw[i];
                pos++;
            }
        }
        printOK();
        printString(Integer.toString(parray.length) + "\n");
        Bits.sendInts(parray, outstream, buffer, myorder);        
        Bits.sendFloats(farray, outstream, buffer, myorder);
    }

    public void processByteOrder() throws IOException {
        String orderString = request.map.get("order");
        if (orderString == null) {
            printInvalid("must supply order for byte order request\n");
            return;
        }
        if (orderString.equals("big")) {
            clientorder = ByteOrder.BIG_ENDIAN;
            printString(myorder == ByteOrder.BIG_ENDIAN ? "big\n" : "little\n");
        } else if (orderString.equals("little")) {
            clientorder = ByteOrder.LITTLE_ENDIAN;
            printString(myorder == ByteOrder.BIG_ENDIAN ? "big\n" : "little\n");
        } else {
            printInvalid("bad byte order " + orderString);
        }
    }
    public String toString() {         
        return "ServerTask user=" + username + " remoteip=" + socket.getInetAddress() + ":" + socket.getPort();
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
                } catch (NullPointerException e) {
                    System.err.println("No password for " + uname);
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
