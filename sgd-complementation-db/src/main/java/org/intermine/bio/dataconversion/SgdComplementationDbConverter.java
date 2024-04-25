package org.intermine.bio.dataconversion;
/*
 * Copyright (C) 2002-2016 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.HashMap;
import java.util.Map;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;


import org.intermine.objectstore.ObjectStoreException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.sql.Database;
import org.intermine.xml.full.Item;
import org.apache.commons.lang.StringUtils;


/**
 *
 * @author
 */
public class SgdComplementationDbConverter extends BioDBConverter {
    //
    private static final String DATASET_TITLE = "Yeast Complementation";
    private static final String DATA_SOURCE_NAME = "SGD-BioGRID curated complementation";
    private final Map<String, Item> genes = new HashMap<String, Item>();
    private final Map<String, Item> publications = new HashMap<String, Item>();
    private static final String TAXON_ID = "4932";
    private static final String H_TAXON_ID = "9606";
    private Item yorganism;
    private Item horganism;
    private String licence;
    private static final SgdComplementationDbProcessor PROCESSOR = new SgdComplementationDbProcessor();


    /**
     * Construct a new SgdComplementationDbConverter.
     * @param database the database to read from
     * @param model the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle Items created
     */
    public SgdComplementationDbConverter(Database database, Model model, ItemWriter writer) throws ObjectStoreException{
        super(database, model, writer, DATA_SOURCE_NAME, DATASET_TITLE);
        yorganism = createItem("Organism");
        yorganism.setAttribute("taxonId", TAXON_ID);
        store(yorganism);
        horganism = createItem("Organism");
        horganism.setAttribute("taxonId", H_TAXON_ID);
        store(horganism);
    }


    /**
     * {@inheritDoc}
     */
    public void process() throws Exception {
        Connection connection = getDatabase().getConnection();
        processComplements(connection);
        storeGenes();
    }


    /**
     * {@inheritDoc}
     */
    public void processComplements(Connection connection) throws SQLException, ObjectStoreException {

        ResultSet res = PROCESSOR.getComplements(connection);

        while (res.next()) {

            String yeastGene = res.getString("yeast_gene");
            String pmid = res.getString("pmid");
            String direction = res.getString("direction");
            String dbxref_id = res.getString("dbxref_id");
            String notes = res.getString("curator_comment");
            String source= res.getString("format_name");

            System.out.println("Processing line..." + yeastGene + "   "+ dbxref_id);

            Item ygene = getGeneItem(yeastGene, "secondaryIdentifier", yorganism);
            Item hgene = getGeneItem(dbxref_id, "secondaryIdentifier", horganism);

            getComplement(notes, direction, source, pmid, ygene, hgene);
            getComplement(notes, direction, source, pmid, hgene, ygene);

        }
    }
    /**
     *
     * @param c
     * @param n
     * @param source
     * @param pmid
     * @return
     * @throws ObjectStoreException
     */
    private void getComplement(String n, String d, String s, String pmid,
                               Item yg, Item hg) throws ObjectStoreException {

		Item pub = publications.get(pmid);
		if (pub == null) {
			pub = createItem("Publication");
			pub.setAttribute("pubMedId", pmid);
			publications.put(pmid, pub);
			store(pub);
		}

        Item comp = createItem("Complement");

        if(StringUtils.isNotEmpty(n)) comp.setAttribute("notes", n);
        comp.setAttribute("source", s);
        comp.setAttribute("direction", d);
        comp.setReference("publication", pub);
        comp.setReference("gene", yg);
        comp.setReference("complement", hg);
        store(comp);

    }

    /**
     *
     * @param geneId
     * @return
     * @throws ObjectStoreException
     */
    private Item getGeneItem(String geneId, String identifier, Item org)
            throws ObjectStoreException {

        Item gene = genes.get(geneId);

        if (gene == null) {
            gene = createItem("Gene");
            genes.put(geneId, gene);
            gene.setAttribute(identifier, geneId);
            gene.setReference("organism", org);
        }

        return gene;

    }

    /**
     *
     * @throws ObjectStoreException
     */

    private void storeGenes() throws ObjectStoreException {
        for (Item gene : genes.values()) {
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }
    /**
     * {@inheritDoc}
     */
    @Override
    public String getDataSetTitle(String taxonId) {
        return DATASET_TITLE;
    }
    /**
     * Set the data licence for these data.
     *
     * @param licence should be URI to data licence.
     */
    public void setLicence(String licence) {
        this.licence = licence;
    }
    /**
     * Get the data licence for these data.
     *
     * @return URI to data licence.
     */
    public String getLicence() {
        return licence;
    }

}


