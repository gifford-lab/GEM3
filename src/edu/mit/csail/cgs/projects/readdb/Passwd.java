package edu.mit.csail.cgs.projects.readdb;

import java.security.MessageDigest;


/**
 * usage:
 *   java Passwd username passwd
 * prints a line for inclusion in a password file
 */

public class Passwd {
    public static void main(String args[]) throws Exception {
        String username = args[0];
        String password = args[1];
        MessageDigest digest = MessageDigest.getInstance("SHA-1");
        password = new String(digest.digest((username + password).getBytes()));
        System.out.println(username + ":" + password);
    }


}