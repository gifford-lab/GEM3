package edu.mit.csail.cgs.projects.readdb;

public class ClientException extends Exception {

    public ClientException (String reason) {
        super(reason);
    }
}