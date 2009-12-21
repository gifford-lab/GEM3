package edu.mit.csail.cgs.projects.readdb;

import java.io.IOException;

public interface Closeable {

    public void close() throws IOException;

}