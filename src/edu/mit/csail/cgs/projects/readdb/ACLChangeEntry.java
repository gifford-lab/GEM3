package edu.mit.csail.cgs.projects.readdb;

/**
 * Represents a change to an ACL.  operation (o) is either ADD or REMOVE
 * and acl (a) is either READ, WRITE, or ADMIN.
 */
public class ACLChangeEntry {

    public static final int ADD = 0, REMOVE = 1;
    public static final String[] ops = {"add","delete"};
    public static final int READ = 0, WRITE = 1, ADMIN = 2;
    public static final String[] acls = {"read","write","admin"};

    private int operation;
    private int acl;
    private String princ;

    public ACLChangeEntry(int o, int a, String p) throws IllegalArgumentException {
        operation = o;
        acl = a;
        princ = p;
        if (o != ADD && o != REMOVE) {
            throw new IllegalArgumentException("Bad operation value : " + o);
        }
        if (a != READ && a != WRITE && a != ADMIN) {
            throw new IllegalArgumentException("Bad acl value : " + a);
        }
    }
    public ACLChangeEntry(String o, String a, String p) throws IllegalArgumentException {
        operation = opCode(o);
        acl = aclCode(a);
        princ = p;
    }
    public String toString() {
        return String.format("%s %s %s",
                             princ, ops[operation], acls[acl]);
    }
    public static int opCode(String opName) throws IllegalArgumentException {
        for (int i = 0; i < ops.length; i++) {
            if (opName.equals(ops[i])) {
                return i;
            }
        }
        throw new IllegalArgumentException("Unknown operation " + opName);
    }
    public static int aclCode(String aclName) throws IllegalArgumentException {
        for (int i = 0; i < acls.length; i++) {
            if (aclName.equals(acls[i])) {
                return i;
            }
        }
        throw new IllegalArgumentException("Unknown acl " + aclName);
    }

}