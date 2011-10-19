package edu.mit.csail.cgs.projects.readdb;

import java.net.*;
import java.util.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import javax.security.sasl.*;
import javax.security.auth.callback.*;

/**
 * <p>API for remote access to the readdb server.
 *
 *  <p>Calls throw IOException on network errors.
 *  Calls throw ClientException on other errors (authentication, authorization, invalid request ,etc.
 *
 *  <p>Client generally assumes that the hit positions are the 5' end of the
 *  hit. 
 *
 * <p>Client IS NOT REENTRANT.  Do not overlap calls to a single Client object.
 *
 * <p>Most method parameters that are object types (eg Integer, Boolean) are optional.  If a null value
 * is passed then no filtering is done based on that parameter.  

 * <p>Standard parameters shared across methods:
 * <ul>
 * <li> alignid is the name of the alignment.
 * <li> isPaired specifies whether to work on single-ended reads (false) or paired-end reads (true)
 * <li> isLeft is isPaired is true, then isLeft specifies whether to work on the left read (true)(
 *     or right read (false) of the pair
 * <li> plusStrand specifies whether to return only reads on the plus strand (true) or minus strand (false). null
 *     means that reads on both strands should be returned.
 * <li> minWeight specifies the minimum weight of reads (or read pairs) to be returned or included in 
 *     the histogram
 * <li> start, stop specify the lowest (inclusive) and highest (exclusive) coordinates of reads
 *     that should be included in the results
 * </ul>
 *
 * 
 */
public class Client implements ReadOnlyClient {    

    public static final String SaslMechanisms[] = {"CRAM-MD5","DIGEST-MD5"};
    /* socket is the socket to talk to the server.  outstream and instream are from
       the socket
        */       
    Socket socket;
    OutputStream outstream;
    BufferedInputStream instream;
    /* temporary space for receiving data; contents not persistent between method calls */
    byte[] buffer;
    private static final int BUFFERLEN = 8192*20;
    private Request request;
    private boolean printErrors;

    /** Connects to a Readdb server on the specified host and port using the specified 
     * username and password.
     * @throws IOException on network errors
     * @throws ClientException if the client cannot authenticate to the server
     */
    public Client (String hostname,
                   int portnum,
                   String username,
                   String passwd) throws IOException, ClientException {
        init(hostname,portnum,username,passwd);
    }
    /**
     * Creates the default connection
     * as specified by ~/.readdb_passwd or a readdb_passwd found in the classpath
     * Must have keys hostname, port, username, and passwd in a format
     * that java.util.PropertyResourceBundle can read
     *
     * @throws IOException on network errors
     * @throws ClientException if the client cannot authenticate to the server
     */
    public Client() throws IOException, ClientException {
        String homedir = System.getenv("HOME");
        String basename = "readdb_passwd";
        if (System.getenv("READDBROLE") != null) {
            basename = System.getenv("READDBROLE") + basename;
        }
        String fname = homedir + "/." + basename;
        File propfile = new File(fname);
        PropertyResourceBundle bundle = null;
        if (propfile.exists() && propfile.canRead()) {
            bundle = new PropertyResourceBundle(new FileInputStream(propfile));
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening readdb properties from " + propfile);
            }
        } else {
            ClassLoader cl = ClassLoader.getSystemClassLoader();
            URL url = cl.getResource(basename);
            if (System.getenv("DEBUGPW") != null) {
                System.err.println("Opening readdb properties from " + url);
            }            
            if (url != null) {
                bundle = new PropertyResourceBundle(url.openStream());
            } else {
                throw new IOException("Can't read connection properties from " + url);
            }
        }
        String hostname = bundle.getString("hostname");
        String port = bundle.getString("port");
        String username = bundle.getString("username");
        String password = bundle.getString("passwd");
        init(hostname, Integer.parseInt(port), username, password);
    }
    
    private void init(String hostname,
                   int portnum,
                   String username,
                   String passwd) throws IOException, ClientException {
        socket = new Socket(hostname,portnum);
        socket.setTcpNoDelay(true);
        socket.setSendBufferSize(BUFFERLEN);
        socket.setReceiveBufferSize(BUFFERLEN);
        /* linger = true, time = 0 means that the other side gets a reset if the
           socket is closed on our end (eg, java exits).  We turn this off just before
           sending a "bye" to allow for a graceful shutdown.  But the RST in
           other cases lets the server figure out that we've disappeared
        */
        socket.setSoLinger(true,0);
        outstream = socket.getOutputStream();
        instream = new BufferedInputStream(socket.getInputStream());
        buffer = new byte[BUFFERLEN];
        
        if (!authenticate(hostname,username,passwd)) {
            throw new ClientException("Authentication Exception Failed");
        }
        request = new Request();
        printErrors = false;
    }
    /**
     * Determines whether the client will print error messages to STDERR.  Useful for debugging 
     * but may produce unwanted screen output.
     */
    public void printErrors(boolean b) {printErrors = b;}
    /**
     * performs the SASL authentication exchange with the server.  currently called by the constructor
     */
    private boolean authenticate(String hostname, String username, String password) throws IOException {
        SaslClient sasl = null;
        try {
            sendString(username + "\n");
            Map<String,String> props = new HashMap<String,String>();
            props.put("Sasl.POLICY_NOPLAINTEXT","true");
            props.put("Sasl.POLICY_NOANONYMOUS","true");
            sasl = Sasl.createSaslClient(SaslMechanisms,
                                         username,
                                         "readdb",
                                         hostname,
                                         props,
                                         new ClientCallbackHandler(username,password));
            if (sasl == null) { return false; }
            byte[] response = (sasl.hasInitialResponse() ? sasl.evaluateChallenge(new byte[0]) : new byte[0]);
            byte continueByte = 1;
            while (continueByte != 0) {
                outstream.write((response.length + "\n").getBytes());
                outstream.write(response);
                outstream.flush();
                //                System.err.println("Sent " + response.length + " bytes");
                int length = Integer.parseInt(readLine());
                //                System.err.println("Waiting to read " + length + " bytes");
                byte[] challenge = new byte[length];
                int read = 0;
                while (read < length) {
                    read += instream.read(challenge, read, length - read);
                    //                    System.err.println("   read " + read);
                }
                /* the continueByte tells us whether the server
                   expects to do another round.  Necessary because sometimes
                   isComplete() returned true here but the server wasn't done
                */
                continueByte = (byte)instream.read();
                //                System.err.println("Continue byte is " + continueByte);
                if (!sasl.isComplete()) { 
                    response = sasl.evaluateChallenge(challenge);
                } else {
                    response = new byte[0];
                }
            }
            sasl.dispose();
            String status = readLine();
            return (status.equals("authenticated as " + username));
        } catch (SaslException e) {
            e.printStackTrace();
            if (sasl != null) {
                sasl.dispose();
            }
            return false;
        }
    }
    
    /** sends a string to the server and flushes the socket 
     */
    private void sendString(String s) throws IOException {
        //System.err.println("SENDING " + s);
        outstream.write(s.getBytes());
        outstream.flush();
    }
    /** reads one line from the server.  blocking.
     */
    private String readLine() throws IOException {
        int pos = 0;
        int i;
        while ((i = instream.read()) != -1) {
            if (i == '\n') {
                break;
            } else {
                buffer[pos++] = (byte)i;
            }
        }
        String out = new String(buffer,0,pos);
        //System.err.println("READ " + out);
        return out;
    }
    /**
     * Tells the server to shut itself down.  Use this to stop the server process.
     * @throws IOException on network errors
     * @throws ClientException if the user isn't authorized to shut the server down.
     */
    public void shutdown() throws IOException, ClientException{
        request.clear();
        request.type = "shutdown";
        sendString(request.toString());
    }
    /** this was to fix a bug in the server.  You shouldn't need it for general use.
     * Regenerate the index for this alignment and chromosome
     */
    public void reIndex(String align, int chrom, boolean paired) throws IOException, ClientException {
        request.clear();
        request.type = "reindex";
        request.alignid = align;
        request.chromid = chrom;
        request.isPaired = paired;
        sendString(request.toString());
        outstream.flush();
        String response = readLine();
        if (!response.equals("OK")) {
            System.out.println(response);
        }
    }
    /** This was to fix a bug in the server.  You shouldn't need it for general use.
     * Resort the hits for a single-ended alignment and regenerate the index.
     */
    public void checksort(String align, int chrom) throws IOException, ClientException {
        request.clear();
        request.type = "checksort";
        request.alignid = align;
        request.chromid = chrom;
        sendString(request.toString());
        outstream.flush();
        String response = readLine();
        if (!response.equals("OK")) {
            System.out.println(response);
        }
    }
    /**
     * Stores a set of SingleHit objects (representing an un-paired or single-ended read
     * aligned to a genome) in the specified alignment.  The hits are appended
     * to any hits that have already been stored in the alignment.
     
     */
    public void storeSingle(String alignid, List<SingleHit> allhits) throws IOException, ClientException {
        int step = 20000000;
        for (int pos = 0; pos < allhits.size(); pos += step) {
            Map<Integer, List<SingleHit>> map = new HashMap<Integer,List<SingleHit>>();
            for (int i = pos; i < pos + step && i < allhits.size(); i++) {
                SingleHit h = allhits.get(i);
                if (!map.containsKey(h.chrom)) {
                    map.put(h.chrom, new ArrayList<SingleHit>());
                }
                map.get(h.chrom).add(h);
            }
            for (int chromid : map.keySet()) {
                List<SingleHit> hits = map.get(chromid);
                Collections.sort(hits);
                int chunk = step;
                for (int startindex = 0; startindex < hits.size(); startindex += chunk) {
                    int count = ((startindex + chunk) < hits.size()) ? chunk : (hits.size() - startindex);
                    request.clear();
                    request.type="storesingle";
                    request.alignid=alignid;
                    request.chromid = chromid;
                    request.map.put("numhits",Integer.toString(count));
                    sendString(request.toString());
                    String response = readLine();
                    if (!response.equals("OK")) {
                        System.err.println("not-OK response to request: " + response);
                        System.err.println("request was " + request);
                        throw new ClientException(response);
                    }
                    int[] ints = new int[count];
                    for (int i = startindex; i < startindex + count; i++) {
                        ints[i - startindex] = hits.get(i).pos;
                    }        
                    Bits.sendInts(ints, outstream,buffer);
                    float[] floats = new float[count];
                    for (int i = startindex; i < startindex + count; i++) {
                        floats[i - startindex] = hits.get(i).weight;
                        ints[i - startindex] = Hits.makeLAS(hits.get(i).length, hits.get(i).strand);
                    }
                    Bits.sendFloats(floats, outstream,buffer);
                    Bits.sendInts(ints, outstream,buffer);
            
                    System.err.println("Sent " + count + " hits to the server for " + chromid + "," + alignid);
                    outstream.flush();        
                    response = readLine();
                    if (!response.equals("OK")) {
                        throw new ClientException(response);
                    }
                }
            }
        }
    }
    /**
     * Stores a set of PairedHit objects (representing an paired-ended read
     * aligned to a genome) in the specified alignment.  The hits are appended
     * to any hits that have already been stored in the alignment
     */
    public void storePaired(String alignid, List<PairedHit> allhits) throws IOException, ClientException {
        Map<Integer, List<PairedHit>> map = new HashMap<Integer,List<PairedHit>>();
        for (PairedHit h : allhits) {
            if (!map.containsKey(h.leftChrom)) {
                map.put(h.leftChrom, new ArrayList<PairedHit>());
            }
            map.get(h.leftChrom).add(h);
        }
        for (int chromid : map.keySet()) {            
            List<PairedHit> hits = map.get(chromid);
            System.err.println("SENDING PAIRED HITS n="+hits.size() + " for chrom " + chromid);
            int chunk = 10000000;
            for (int startindex = 0; startindex < hits.size(); startindex += chunk) {
                int count = ((startindex + chunk) < hits.size()) ? chunk : (hits.size() - startindex);

                request.clear();
                request.type="storepaired";
                request.alignid=alignid;
                request.chromid=chromid;
                request.isLeft=true;
                request.map.put("numhits",Integer.toString(count));
                sendString(request.toString());
                String response = readLine();
                if (!response.equals("OK")) {
                    System.err.println("not-OK response to request: " + response);
                    System.err.println("request was " + request);
                    throw new ClientException(response);
                }
                int[] ints = new int[count];
                for (int i = startindex; i < startindex + count; i++) {
                    ints[i-startindex] = hits.get(i).leftPos;
                }        
                Bits.sendInts(ints, outstream,buffer);
                float[] floats = new float[count];
                for (int i = startindex; i < startindex + count; i++) {
                    floats[i-startindex] = hits.get(i).weight;
                    ints[i-startindex] = Hits.makeLAS(hits.get(i).leftLength, hits.get(i).leftStrand,
                                                      hits.get(i).rightLength, hits.get(i).rightStrand);

                }
                Bits.sendFloats(floats, outstream,buffer);
                Bits.sendInts(ints, outstream,buffer);
                for (int i = startindex; i < startindex + count; i++) {
                    ints[i-startindex] = hits.get(i).rightChrom;
                }        
                Bits.sendInts(ints, outstream,buffer);
                for (int i = startindex; i < startindex + count; i++) {
                    ints[i-startindex] = hits.get(i).rightPos;
                }        
                Bits.sendInts(ints, outstream,buffer);
                System.err.println("Sent " + count + " hits to the server");
                outstream.flush();        
                response = readLine();
                if (!response.equals("OK")) {
                    if (printErrors) {
                        System.err.println("not-OK response to request: " + response);
                        System.err.println("request was " + request);
                    }
                    throw new ClientException(response);
                }
            }
        }
    }
    
    /** Returns true if the alignment and chromosome exist and are accessible
     * to the user.  Returns false if they don't exist or if they aren't accessible
     */
    public boolean exists(String alignid) throws IOException {
        request.clear();
        request.type="exists";
        request.alignid=alignid;
        sendString(request.toString());
        String response = readLine();
        if (response.equals("exists")) {
            return true;
        } else if (response.equals("unknown")) {
            return false;
        } else {
            return false;
        }
    }
    /**
     * Deletes an alignment (all chromosomes).  isPaired specifies whether to delete
     * the paired or single ended reads.
     */
    public void deleteAlignment(String alignid, Boolean isPaired) throws IOException, ClientException {
        request.clear();
        request.type="deletealign";
        request.isPaired = isPaired;
        request.alignid=alignid;
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
    }
    /**
     * Returns the set of chromosomes that exist for this alignment. 
     */
    public Set<Integer> getChroms(String alignid, boolean isPaired, Boolean isLeft) throws IOException, ClientException {
        request.clear();
        request.type="getchroms";
        request.isLeft = isLeft;
        request.isPaired = isPaired;
        request.alignid=alignid;
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numchroms = Integer.parseInt(readLine());
        Set<Integer> output = new HashSet<Integer>();
        while (numchroms-- > 0) {
            output.add(Integer.parseInt(readLine()));
        }
        return output;
    }
    /**
     * Returns the total number of hits in this alignment.  
     */
    public int getCount(String alignid, boolean isPaired, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
        int count = 0;
        for (int c : getChroms(alignid, isPaired, isLeft)) {
            count += getCount(alignid, c, isPaired, null,null,null,isLeft,plusStrand);
        }
        return count;
    }
    /**
     * Returns the sum of the weights of all hits in this alignment
     */
    public double getWeight(String alignid, boolean isPaired, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
        double total = 0;
        for (int c : getChroms(alignid, isPaired, isLeft)) {
            total += getWeight(alignid, c, isPaired, null, null, null, isLeft, plusStrand);
        }
        return total;
    }

    /** returns the total number of hits on the specified chromosome in the alignment.
     * Any of the object parameters can be set to null to specify "no value"
     */
    public int getCount(String alignid, int chromid, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand)  throws IOException, ClientException {
        request.clear();
        request.type="count";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.isLeft = isLeft == null ? true : isLeft;
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return numhits;
    }
    /** returns the total weight on the specified chromosome in this alignment
     */
    public double getWeight(String alignid, int chromid, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="weight";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.isLeft = isLeft == null ? true : isLeft;
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        return Double.parseDouble(readLine());
    }

    /** 
     * returns the sorted (ascending order) hit positions in the specified range of a chromosome,alignment pair.
     */ 
    public int[] getPositions(String alignid, int chromid, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="gethits";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.isLeft = isLeft;
        request.map.put("wantpositions","1");
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return Bits.readInts(numhits, instream, buffer);        
    }
    /** 
     * returns the hit weights in the specified range of a chromosome,alignment pair.  The weights
     * will be in the same order as the sorted positions returned by getPositions()
     */ 
    public float[] getWeightsRange(String alignid, int chromid, boolean paired, Integer start, Integer stop, Float minWeight, Boolean isLeft, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="gethits";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.isLeft = isLeft;
        request.map.put("wantweights","1");
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return Bits.readFloats(numhits, instream, buffer);        
    }
    public List<SingleHit> getSingleHits(String alignid, int chromid, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="gethits";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = false;
        request.map.put("wantpositions","1");
        request.map.put("wantweights","1");
        request.map.put("wantlengthsandstrands","1");
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        List<SingleHit> output = new ArrayList<SingleHit>();
        int numhits = Integer.parseInt(readLine());
        for (int i = 0; i < numhits; i++) {
            output.add(new SingleHit(chromid,0,(float)0.0,false,(short)0));
        }
        IntBP ints = new IntBP(numhits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(ints.bb, rbc);
        for (int i = 0; i < numhits; i++) {
            output.get(i).pos = ints.get(i);
        }
        FloatBP floats = new FloatBP(numhits);
        Bits.readBytes(floats.bb, rbc);
        for (int i = 0; i < numhits; i++) {
            output.get(i).weight = floats.get(i);
        }
        Bits.readBytes(ints.bb, rbc);
        for (int i = 0; i < numhits; i++) {
            int j = ints.get(i);
            SingleHit h = output.get(i);
            h.length = Hits.getLengthOne(j);
            h.strand = Hits.getStrandOne(j);
        }
        return output;
    }
    public List<PairedHit> getPairedHits(String alignid, int chromid, boolean isLeft, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="gethits";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isLeft = isLeft;
        request.isPaired = true;
        request.map.put("wantpositions","1");
        request.map.put("wantweights","1");
        request.map.put("wantlengthsandstrands","1");
        request.map.put("wantotherchroms","1");
        request.map.put("wantotherpositions","1");
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(String.format("align %s chrom %d: %s", alignid, chromid, response));
        }
        List<PairedHit> output = new ArrayList<PairedHit>();
        int numhits = Integer.parseInt(readLine());
        for (int i = 0; i < numhits; i++) {
            output.add(new PairedHit(chromid,0,false,(short)0,
                                     chromid,0,false,(short)0,(float)0));
        }
        IntBP ints = new IntBP(numhits);
        ReadableByteChannel rbc = Channels.newChannel(instream);
        Bits.readBytes(ints.bb, rbc);
        if (isLeft) {
            for (int i = 0; i < numhits; i++) {
                output.get(i).leftPos = ints.get(i);
            }
        } else {
            for (int i = 0; i < numhits; i++) {
                output.get(i).rightPos = ints.get(i);
            }
        }
        FloatBP floats = new FloatBP(numhits);
        Bits.readBytes(floats.bb, rbc);
        for (int i = 0; i < numhits; i++) {
            output.get(i).weight = floats.get(i);
        }

        Bits.readBytes(ints.bb, rbc);
        if (isLeft) {
            for (int i = 0; i < numhits; i++) {
                int j = ints.get(i);
                PairedHit h = output.get(i);                
                h.leftLength = Hits.getLengthOne(j);
                h.leftStrand = Hits.getStrandOne(j);
                h.rightLength = Hits.getLengthTwo(j);
                h.rightStrand = Hits.getStrandTwo(j);
            }
        } else {
            for (int i = 0; i < numhits; i++) {
                int j = ints.get(i);
                PairedHit h = output.get(i);                
                h.leftLength = Hits.getLengthTwo(j);
                h.leftStrand = Hits.getStrandTwo(j);
                h.rightLength = Hits.getLengthOne(j);
                h.rightStrand = Hits.getStrandOne(j);
            }
        }
        Bits.readBytes(ints.bb, rbc);
        if (isLeft) {
            for (int i = 0; i < numhits; i++) {
                output.get(i).rightChrom = ints.get(i);
            }
        } else {
            for (int i = 0; i < numhits; i++) {
                output.get(i).leftChrom = ints.get(i);
            }
        }

        Bits.readBytes(ints.bb, rbc);
        if (isLeft) {
            for (int i = 0; i < numhits; i++) {
                output.get(i).rightPos = ints.get(i);
            }
        } else {
            for (int i = 0; i < numhits; i++) {
                output.get(i).leftPos = ints.get(i);
            }
        }


        return output;
    }


    /**
     * returns a TreeMap from positions to counts representing a histogram
     * of the hits in a range with the specified binsize.  Bins with a count
     * of zero are not included in the output.
     * 
     * Ex: getHistgram("foo","1+",1,100,10)
     *     6: 5
     *    16: 0
     *    36: 30
     *
     * minweight is the minimum weight for reads to be included in the histogram.
     *
     * dedup is the limit on how many times reads with any given 5' position will be counted.
     *  A value of zero means no limit.  A limit of, eg, 2, means that at most two reads at
     *  any 5' position will be included in the output.  For weighted histograms, the choice of reads
     *  included is unspecified.  For methods that operate on a set of alignments, this many
     *  reads from each alignment will be included.
     *
     * Normally, a read is only counted in a single bin as defined by its position (generally
     * the 5' end of the read).  A non-zero read-Extension counts the read in any
     * bin that you cross within that many bases of the read's position.  A negative value
     * goes backwards (smaller coordinates) and a larger value goes forward.  You need
     * to get the sign right depending on the strandedness of the chromosome that you're
     * working with.
     */
    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getHistogram(alignid, chromid, paired, doReadExtension,binsize,0,start,stop,minWeight,plusStrand,true);
    }
    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand, boolean isLeft) throws IOException, ClientException {
        request.clear();
        request.type="histogram";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.isLeft = isLeft;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.map.put("binsize",Integer.toString(binsize));
        if (dedup > 0) {
            request.map.put("dedup",Integer.toString(dedup));
        }
        if (doReadExtension) {
            request.map.put("extension","1");
        }
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numints = Integer.parseInt(readLine());
        int out[] = Bits.readInts(numints, instream, buffer);
        TreeMap<Integer,Integer> output = new TreeMap<Integer,Integer>();
        for (int i = 0; i < out.length; i += 2) {
            output.put(out[i], out[i+1]);
        }
        return output;
    }
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getWeightHistogram(alignid, chromid, paired, doReadExtension, binsize, 0, start,stop,minWeight,plusStrand);
    }
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        request.clear();
        request.type="weighthistogram";
        request.alignid=alignid;
        request.chromid=chromid;
        request.start = start;
        request.end = stop;
        request.minWeight = minWeight;
        request.isPlusStrand = plusStrand;
        request.isPaired = paired;
        request.map.put("binsize",Integer.toString(binsize));
        if (dedup > 0) {
            request.map.put("dedup",Integer.toString(dedup));
        }

        if (doReadExtension) {
            request.map.put("extension","1");
        }
        sendString(request.toString());        
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        int numints = Integer.parseInt(readLine());
        int out[] = Bits.readInts(numints, instream, buffer);
        float weight[] = Bits.readFloats(numints, instream,buffer);
        TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
        for (int i = 0; i < out.length; i++) {
            output.put(out[i], weight[i]);
        }
        return output;
    }

    public TreeMap<Integer,Integer> getHistogram(Collection<String> alignids, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getHistogram(alignids,chromid,paired,doReadExtension,binsize,0,start,stop,minWeight,plusStrand);
    }
    public TreeMap<Integer,Integer> getHistogram(Collection<String> alignids, int chromid, boolean paired, boolean doReadExtension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        TreeMap<Integer,Integer> output = null;
        for (String alignid : alignids) {
            TreeMap<Integer,Integer> o = getHistogram(alignid,chromid,paired,doReadExtension,binsize,dedup,start,stop,minWeight,plusStrand,true);
            for (int k : o.keySet()) {
                if ((k - start - binsize / 2) % binsize != 0 ) {
                    System.err.println(String.format("Bad key %d for binsize %d and start %d in %s,%d",
                                                     k, binsize, start, alignid,chromid));
                }
            }

            if (output == null) {
                output = o;
            } else {
                for (int k : o.keySet()) {
                    if (output.containsKey(k)) {
                        output.put(k, output.get(k) + o.get(k));
                    } else {
                        output.put(k,o.get(k));
                    }
                }
            }            
        }
        return output;
    }
    public TreeMap<Integer,Float> getWeightHistogram(Collection<String> alignids, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        return getWeightHistogram(alignids,chromid,paired,doReadExtension,binsize,0,start,stop,minWeight,plusStrand);
    }
    public TreeMap<Integer,Float> getWeightHistogram(Collection<String> alignids, int chromid, boolean paired, boolean doReadExtension, int binsize, int dedup, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        TreeMap<Integer,Float> output = null;
        for (String alignid : alignids) {
            TreeMap<Integer,Float> o = getWeightHistogram(alignid,chromid,paired,doReadExtension,binsize,dedup,start,stop,minWeight,plusStrand);
            if (output == null) {
                output = o;
            } else {
                for (int k : o.keySet()) {
                    if (output.containsKey(k)) {
                        output.put(k, output.get(k) + o.get(k));
                    } else {
                        output.put(k,o.get(k));
                    }
                }
            }            
        }
        return output;
    }

    /**
     * Returns a Map from READ, WRITE, and ADMIN to lists of principals that have those privileges on the specified alignment.
     */
    public Map<String,Set<String>> getACL(String alignid) throws IOException, ClientException {
        request.clear();
        request.type="getacl";
        request.alignid=alignid;
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
        Map<String,Set<String>> output = new HashMap<String,Set<String>>();
        fillPartACL(output);
        fillPartACL(output);
        fillPartACL(output);
        return output;
    }
    /* fills one section of the acl output data structure.  A section
       is either read, write, or admin. 
    */
    private void fillPartACL(Map<String,Set<String>> output) throws IOException {
        String type = readLine();
        int entries = Integer.parseInt(readLine());
        Set<String> out = new HashSet<String>();
        while (entries-- > 0) {
            out.add(readLine());
        }
        output.put(type,out);
    }
    /**
     * Applies the specified ACLChangeEntry objects to the acl for this experiment/chromosome.
     */
    public void setACL(String alignid, Set<ACLChangeEntry> changes) throws IOException, ClientException {
        request.clear();
        request.type="setacl";
        request.alignid=alignid;
        for (ACLChangeEntry a : changes) {
            request.list.add(a.toString());
        }    
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }   
    }
    /**
     * Adds the specified user (princ) to a group.
     */
    public void addToGroup(String princ, String group) throws IOException, ClientException {
        request.clear();
        request.type="addtogroup";
        request.map.put("princ",princ);
        request.map.put("group",group);
        sendString(request.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            if (printErrors) {
                System.err.println("not-OK response to request: " + response);
                System.err.println("request was " + request);
            }
            throw new ClientException(response);
        }
    }
    /**
     * Closes this connection to the server.
     */
    public void close() {
        if (socket == null) {
            return;
        }
        try {
            socket.setSoLinger(false,0);
            request.clear();
            request.type="bye";
            sendString(request.toString());
            outstream.close();
            outstream = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            instream.close();
            instream = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            socket.close();
            socket = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }    
    
}

/**
 * SASL callback handler for the authentication exchange
 */
class ClientCallbackHandler implements CallbackHandler {
    private String name, pass;
    public ClientCallbackHandler(String n, String p) {name = n; pass = p;}
    public void handle(Callback[] callbacks) {
        for (int i = 0; i < callbacks.length; i++) {
            if (callbacks[i] instanceof NameCallback) {
                NameCallback nc = (NameCallback)callbacks[i];
                nc.setName(name);
            }
            if (callbacks[i] instanceof PasswordCallback) {
                PasswordCallback pc = (PasswordCallback)callbacks[i];
                pc.setPassword(pass.toCharArray());
            }
        }
    }
}