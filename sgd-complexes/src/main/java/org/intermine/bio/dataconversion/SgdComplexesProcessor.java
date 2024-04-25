package org.intermine.bio.dataconversion;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import org.apache.log4j.Logger;


public class SgdComplexesProcessor {

    private static final Logger LOG = Logger.getLogger(SgdComplexesProcessor.class);

    /**
     * Return the results of running a query for genes
     * @param connection the connection
     * @return the results
     * @throws SQLException if there is a database problem
     */
    protected ResultSet getComplexes(Connection connection) throws SQLException {

        String query = "select cd.dbentity_id, complex_accession, intact_id, cd.systematic_name, cd.description, properties, " +
                " db.display_name as db_display_name,  ea.format_name as eco_id, array_agg(distinct pmid) as pmid, " +
                " array_agg(distinct ca.display_name) as ca_display_name, array_agg(distinct go.format_name) as goids" +
                " from nex.complexdbentity cd " +
                " inner join nex.dbentity db on db.dbentity_id = cd.dbentity_id" +
                " left join nex.eco ea on cd.eco_id = ea.eco_id" +
                " left join nex.complex_alias ca on ca.complex_id = cd.dbentity_id" +
                " left join nex.complex_reference cr on cr.complex_id = cd.dbentity_id" +
                " left join nex.referencedbentity rdb on cr.reference_id = rdb.dbentity_id" +
                " left join nex.complex_go cg on cg.complex_id = cd.dbentity_id" +
                " left join nex.go go on go.go_id = cg.go_id" +
                //" and ca.alias_type = 'Synonym'" +
                //" where cd.dbentity_id = 1982773"+
                " group by cd.dbentity_id, complex_accession, intact_id, cd.systematic_name, cd.description, properties,  db.display_name,  ea.eco_id "+
                " order by cd.dbentity_id";

        LOG.info("executing: " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }


    /**
     * Return the results of running a query for genes
     * @param connection the connection
     * @return the results
     * @throws SQLException if there is a database problem
     */
    protected ResultSet getComplexInteractions(Connection connection) throws SQLException {

        String query = "select cd.dbentity_id, complex_accession," +
                " range_start, range_end, stoichiometry, psimi2.display_name as role, psimi.display_name as type," +
                " ldb.sgdid as sgdid_1,  array_agg(ldb2.sgdid) as sgdid_2, i.display_name as interactordisplay, i.format_name as interactorid, ib.display_name" +
                " from nex.complexdbentity cd " +
                " inner join nex.complexbindingannotation cba on cd.dbentity_id = cba.complex_id" +
                " inner join nex.taxonomy t on t.taxonomy_id = cba.taxonomy_id" +
                " inner join nex.interactor i on  cba.interactor_id = i.interactor_id" +
                " left join nex.interactor ib on cba.binding_interactor_id = ib.interactor_id" +
                " left join nex.dbentity ldb on ldb.dbentity_id = i.locus_id" +
                " left join nex.dbentity ldb2 on ldb2.dbentity_id = ib.locus_id" +
                " left join nex.psimi psi on psi.psimi_id = cba.binding_type_id" +
                " left join nex.psimi psimi on psimi.psimi_id = i.type_id" +
                " left join nex.psimi psimi2 on psimi2.psimi_id = ib.role_id" +
                //" where cd.dbentity_id = 1982749"+
                " group by cd.dbentity_id, complex_accession," +
                " range_start, range_end, stoichiometry, psimi2.display_name, psimi.display_name,ldb.sgdid, i.display_name, i.format_name, ib.display_name" +
                " order by cd.dbentity_id, ldb.sgdid";

        LOG.info("executing: " + query);
        Statement stmt = connection.createStatement();
        ResultSet res = stmt.executeQuery(query);
        return res;
    }

}
