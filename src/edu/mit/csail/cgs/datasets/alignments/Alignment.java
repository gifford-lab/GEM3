/*
 * Author: tdanford
 * Date: Aug 25, 2008
 */
package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.DBID;
import edu.mit.csail.cgs.utils.database.Sequence;

/*
create sequence alignment_id;
create table alignment ( 
	id number(10) constraint cst_align_id unique not null,
	params varchar2(1000),
	version number(10) constraint cst_align_version references alignment_version(id),
	score number
);
*/

public interface Alignment { 
	
	public static final String dbRole = "alignment"; // needed for some backwards compatibility issues.  Grr.
	
	public DBID getDBID();
	public String getParams();
	public AlignmentVersion getVersion();
	public Double getScore();
	
	public Collection<AlignBlock> getAlignBlocks();
}

