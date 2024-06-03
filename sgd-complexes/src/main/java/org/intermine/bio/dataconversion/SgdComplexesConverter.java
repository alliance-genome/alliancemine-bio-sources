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

import java.sql.*;
import java.sql.Array;
import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.Iterator;

import org.intermine.objectstore.ObjectStoreException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.sql.Database;
import org.intermine.xml.full.Item;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;


/**
 * 
 * @author
 */
public class SgdComplexesConverter extends BioDBConverter
{
    // 
    private static final String DATASET_TITLE = "Add DataSet.title here";
    private static final String DATA_SOURCE_NAME = "Add DataSource.name here";
    private static final SgdComplexesProcessor PROCESSOR = new SgdComplexesProcessor();
    private String licence;

    private static final Logger LOG = Logger.getLogger(SgdComplexesConverter.class);

    // TODO types (protein and small molecules are processed now) are hardcoded.
    // maybe put this in config file? Or check model to see if type is legal?
    private static final Map<String, String> INTERACTOR_TYPES = new HashMap<String, String>();

    // accession to stored object ID
    private Map<String, String> interactors = new HashMap<String, String>();
    private Map<String, String> synonyms = new HashMap();
    private final Map<String, Item> ecoMap = new HashMap<String, Item>();
    private Map<String, Item> publications = new HashMap();
    private Map<String, String> terms = new HashMap<String, String>();
    private Map<String, Item> complexes = new HashMap();
    private final Map<String, Item> proteins = new HashMap<String, Item>();
    private final Map<String, Item> otherinteractors = new HashMap<String, Item>();

    private static final String TAXON_ID = "559292"; //4932
    private Item yorganism;

    /**
     * Construct a new SgdComplexesConverter.
     * @param database the database to read from
     * @param model the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle Items created
     */
    public SgdComplexesConverter(Database database, Model model, ItemWriter writer) throws ObjectStoreException{
        super(database, model, writer, DATA_SOURCE_NAME, DATASET_TITLE);
        yorganism = createItem("Organism");
        yorganism.setAttribute("taxonId", TAXON_ID);
        store(yorganism);
    }


    /**
     * {@inheritDoc}
     */
    public void process() throws Exception {
        Connection connection = getDatabase().getConnection();
        processComplexes(connection);
        processComplexInteractions(connection);
        storeProteins();
        storeOtherBioentities();
        storeComplexes();
    }


    /**
     * {@inheritDoc}
     */
    public void processComplexes(Connection connection) throws SQLException, ObjectStoreException {

        ResultSet res = PROCESSOR.getComplexes(connection);

        while (res.next()) {
            String productId = res.getString("dbentity_id");
            String complex_accession = res.getString("complex_accession");
            String intact_id = res.getString("intact_id");
            String systematic_name = res.getString("systematic_name");
            String description = res.getString("description");
            String properties = res.getString("properties");
            String display_name = res.getString("db_display_name");
            String eco_id = res.getString("eco_id");
            Array aliases = res.getArray("ca_display_name");
            String[] synonyms = (String[])aliases.getArray();
            Array pubs = res.getArray("pmid");
            Long[] pmids = (Long[])pubs.getArray();
            Array gos = res.getArray("goids");
            String[] goids = (String[])gos.getArray();
            System.out.println("productId is "+ productId);
            processComplex(complex_accession, intact_id, systematic_name, description, properties, display_name, eco_id, synonyms, pmids, goids);
        }
    }

    /**
     *
     * @param accession
     * @param identifier
     * @param systematicName
     * @param description
     * @param properties
     * @param aliases
     * @return
     */
    private void processComplex(String accession, String identifier, String systematicName, String description, String properties,
                                String display_name, String eco_id, String[] aliases, Long[] pmids, String[] goids) throws ObjectStoreException {


        Item complex = createItem("Complex");

        if (StringUtils.isNotEmpty(accession)) {
            complex.setAttribute("accession", accession);
        }
        if (StringUtils.isNotEmpty(identifier)) {
            complex.setAttribute("identifier", identifier);
        }
        if (StringUtils.isNotEmpty(systematicName)) {
            complex.setAttribute("systematicName", systematicName);
        }
        if (StringUtils.isNotEmpty(description)) {
            complex.setAttribute("function", description);
        }
        if (StringUtils.isNotEmpty(properties)) {
            complex.setAttribute("properties", properties);
        }
        if (StringUtils.isNotEmpty(display_name)) {
            complex.setAttribute("name", display_name);
        }
        //add aliases
        for (String alias : aliases) {
            if (StringUtils.isNotEmpty(alias)) {
                String syn = getSynonym(alias);
                complex.addToCollection("synonyms", syn);
            }
        }

        //add pubs
        for ( Long pmid : pmids) {
            if(pmid != null) {
                String pubmedId = Long.toString(pmid);
                if (StringUtils.isNotEmpty(pubmedId)) {
                    String pub = getPublication(pubmedId);
                    complex.addToCollection("publications", pub);
                }
            }
        }
        //add goAnnot
        for ( String goid : goids) {
            if (StringUtils.isNotEmpty(goid)) {
                Item goAnnot = getGoAnnot(goid);
                complex.addToCollection("goAnnotation", goAnnot);
            }
        }
        //evidence
        String ecoId = getEvidence(eco_id);
        complex.setReference("complexEvidence", ecoId);

        complexes.put(accession, complex);

    }

    /**
     *
     * @param goId
     * @return
     * @throws ObjectStoreException
     */
    private Item getGoAnnot(String goId) throws ObjectStoreException {

        String goterm = getTerm(goId);
        Item goAnnotation = createItem("GOAnnotation");
        goAnnotation.setReference("ontologyTerm", goterm);
            try {
                store(goAnnotation);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            return goAnnotation;
    }

    /**
     *
     * @param identifier
     * @return
     * @throws ObjectStoreException
     */

    private String getTerm(String identifier) throws ObjectStoreException {
        String refId = terms.get(identifier);
        if (refId == null) {
            Item ontologyTerm = createItem("GOTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            store(ontologyTerm);
            refId = ontologyTerm.getIdentifier();
            terms.put(identifier, refId);
        }
        return refId;
    }


    /**
     *
     * @param value
     * @return
     * @throws ObjectStoreException
     */
    private String getEvidence(String value) throws ObjectStoreException{

        Item evidence = createItem("ComplexEvidence");

        if (StringUtils.isNotEmpty(value)) {

            Item eco = ecoMap.get(value);
            if(eco == null) {
                eco = createItem("ECOTerm");
                ecoMap.put(value, eco);
                eco.setAttribute("identifier", value);
                try {
                    store(eco);
                } catch (ObjectStoreException e) {
                    throw new ObjectStoreException(e);
                }
            }
            evidence.setReference("ontologyTerm", eco.getIdentifier());
        }

        try {
            store(evidence);
        } catch (ObjectStoreException e) {
            throw new ObjectStoreException(e);
        }

        return evidence.getIdentifier();
    }
    /**
     *
     * @param value
     * @return
     * @throws ObjectStoreException
     */
    private String getSynonym(String value) throws ObjectStoreException {

        if (StringUtils.isEmpty(value)) {
            return null;
        }
        String refId = synonyms.get(value);
        if (refId == null) {
            Item syn = createItem("Synonym");
            syn.setAttribute("value", value);
            refId = syn.getIdentifier();
            try {
                store(syn);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
        return refId;
    }
    private String getPublication(String pubMedId)
            throws ObjectStoreException {
        Item item = publications.get(pubMedId);
        if (item == null) {
             item = createItem("Publication");
            if (StringUtils.isNotEmpty(pubMedId)) {
                item.setAttribute("pubMedId", pubMedId);
            }
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            publications.put(pubMedId, item);
        }
        return item.getIdentifier();
    }

    private void processComplexInteractions(Connection connection) throws ObjectStoreException, SQLException {

        ResultSet res = PROCESSOR.getComplexInteractions(connection);

        while (res.next()) {

            String dbentityId = res.getString("dbentity_id");
            String complex_accession = res.getString("complex_accession");
            String dbentity1 = res.getString("sgdid_1");

            Array interactions = res.getArray("sgdid_2");
            String[] str_interactions = (String[])interactions.getArray();

            String range_start = res.getString("range_start");
            String range_end = res.getString("range_end");
            String stochiometry = res.getString("stoichiometry");
            String role = res.getString("role");
            String type = res.getString("type");
            String int_display_name = res.getString("interactordisplay");
            String int_format_name = res.getString("interactorid");

            Item gene1 = null;

            if(type.equalsIgnoreCase("protein") && dbentity1 != null){
                System.out.println("in type equals protein" + "   "+dbentity1 + "   "+ complex_accession);
                gene1 = getProteinItem(dbentity1);
            }else{
                gene1 = getOtherItem(type, int_display_name, int_format_name);
            }
            processInteractions(complex_accession, gene1, range_start, range_end, stochiometry, role, type, str_interactions);

        }
    }


    /**
     *
     * @param complexId
     * @param ref
     * @param binderRef
     * @param range_start
     * @param range_end
     * @param interactor
     * @throws ObjectStoreException
     */
    private void processInteractions(String complexacc, Item ref, String range_start, String range_end,
                String stochiometry, String role, String type, String[] ints) throws ObjectStoreException {

        Item complex = complexes.get(complexacc);

        //Interactor
        Item interactor = createItem("Interactor");
        if (StringUtils.isNotEmpty(role)) { interactor.setAttribute("biologicalRole", role);}
        if (stochiometry != null) { interactor.setAttribute("stoichiometry", stochiometry); }
        interactor.setAttribute("type", type);
        interactor.setReference("participant", ref);

        //Interaction
        if(ints != null) {

            for (int i = 0; i < ints.length; i++) {
                String binder = ints[i];
                Item binderRef = null;
                if(binder != null) {
                    binderRef = getProteinItem(binder);
                }
                Item interaction = createItem("Interaction");
                interaction.setReference("participant1", ref);
                if (binderRef != null) {interaction.setReference("participant2", binderRef);}
                interaction.setReference("complex", complex);
                store(interaction);
                interactor.addToCollection("interactions", interaction);

                Item detailItem = createItem("InteractionDetail");
                detailItem.setAttribute("type", "physical");
                detailItem.setReference("interaction", interaction);
                //detailItem.setCollection("allInteractors", detail.getAllInteractors());

                processRegions(range_start, range_end, detailItem, ref, binderRef);
                store(detailItem);
            }

        }
        store(interactor);
        if (complex != null) {
            complex.addToCollection("allInteractors", interactor);
        }
    }

    private void processRegions(String startPosition, String endPosition, Item detail, Item ref,
                                Item binderRef)
            throws ObjectStoreException {

        Item location = createItem("Location");

        if (StringUtils.isNotEmpty(startPosition)) { location.setAttribute("start", startPosition);}
        if (StringUtils.isNotEmpty(endPosition))  { location.setAttribute("end", endPosition); }

        location.setReference("locatedOn", ref);
        if(binderRef != null) { location.setReference("feature", binderRef);}

        store(location);

        Item region = createItem("InteractionRegion");
        region.addToCollection("locations", location);
        region.setReference("interaction", detail);
        store(region);

    }

    /**
     *
     * @throws ObjectStoreException
     */

    private void storeComplexes() throws ObjectStoreException {
        for (Item complex : complexes.values()) {
            try {
                store(complex);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

    /**
     *
     * @throws ObjectStoreException
     */

    private void storeProteins() throws ObjectStoreException {
        for (Item protein : proteins.values()) {
            try {
                store(protein);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }


    /**
     *
     * @throws ObjectStoreException
     */

    private void storeOtherBioentities() throws ObjectStoreException {
        for (Item item : otherinteractors.values()) {
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }


    /**
     *
     * @param proteinId
     * @return
     * @throws ObjectStoreException
     */
    private Item getProteinItem(String proteinId)
            throws ObjectStoreException {

        Item protein = proteins.get(proteinId);

        if (protein == null) {
            protein = createItem("Protein");
            proteins.put(proteinId, protein);
            protein.setAttribute("primaryIdentifier", proteinId);
            protein.setReference("organism", yorganism);
        }

        return protein;

    }

    /**
     *
     * @param proteinId
     * @return
     * @throws ObjectStoreException
     */
    private Item getOtherItem(String type, String interactor_dname, String interactor_fname)
            throws ObjectStoreException {

        Item item = otherinteractors.get(interactor_dname);

        if (item == null) {
            item = createItem("BioEntity");
            otherinteractors.put(interactor_dname, item);
            item.setAttribute("primaryIdentifier", interactor_fname);
            item.setAttribute("name", interactor_dname);
            item.setReference("organism", yorganism);
        }

        return item;

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

