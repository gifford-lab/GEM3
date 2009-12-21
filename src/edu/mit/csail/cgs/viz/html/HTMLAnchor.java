/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.io.PrintStream;

/**
 * @author Timothy Danford
 *
 */
public class HTMLAnchor implements HTMLElmt {
	
	private String url;
	private HTMLElmt body;
	
	public HTMLAnchor(String url, HTMLElmt body) { 
		this.url = url;
		this.body = body;
	}
	
	public HTMLAnchor(String url, String body) { 
		this.url = url;
		this.body = new HTMLText(body);
	}
	
	public String getURL() { return url; }
	public void setURL(String url) { this.url = url; }

	public void print(PrintStream ps) {
		ps.print("<a href=\"" + url + "\">");
		body.print(ps);
		ps.print("</a>");
	}

}
