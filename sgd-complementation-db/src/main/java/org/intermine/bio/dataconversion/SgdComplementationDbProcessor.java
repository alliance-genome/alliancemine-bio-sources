package org.intermine.bio.dataconversion;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import org.apache.log4j.Logger;


public class SgdComplementationDbProcessor {

    private static final Logger LOG = Logger.getLogger(SgdComplementationDbProcessor.class);

    /**
     * Return the results of getting complements
     * @param connection the connection
     * @return the results
     * @throws SQLException if there is a database problem
     */
    protected ResultSet getComplements(Connection connection)
            throws SQLException {

        String query = "select db.format_name as yeast_gene, rdb.pmid, "
                + "direction, dbxref_id, curator_comment, s.format_name "
                + "from nex.dbentity db "
                + "inner join nex.functionalcomplementannotation fca on fca.dbentity_id = db.dbentity_id "
                + "inner join nex.referencedbentity rdb on rdb.dbentity_id = fca.reference_id "
                + "inner join nex.source s on s.source_id = fca.source_id ";
        LOG.info("executing: " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }


}

