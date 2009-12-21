package edu.mit.csail.cgs.datasets.function;

import java.sql.SQLException;
import java.util.Collection;

public interface FunctionLoader {

	public Collection<FunctionVersion> getAllVersions()
			throws SQLException;

    public FunctionVersion getVersion(String versionName) 
            throws SQLException;

	public Collection<Category> getCategories(FunctionVersion fv)
			throws SQLException;

    public Collection<Assignment> getAssignments(FunctionVersion version)
        throws SQLException;

    public Collection<Assignment> getAssignments(Category c)
        throws SQLException;

	public Collection<Assignment> getAssignments(String obj, FunctionVersion fv) 
		throws SQLException;
    
    public Collection<Assignment> getAllAssignments(Category c)
        throws SQLException;
    
    public Collection<Category> getChildCategories(Category c) 
        throws SQLException;

    public void close();
}