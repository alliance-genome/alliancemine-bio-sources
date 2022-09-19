package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2018 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class AllianceGeneticInteractionsConverter extends BioFileConverter
{
    private static final String DATASET_TITLE = "Alliance Genetic Interactions data set";
    private static final String DATA_SOURCE_NAME = "Alliance Genetic Interactions";
    private String licence;
    private Map<String, Item> genes = new HashMap();
    private Map<String, Item> publications = new HashMap();
    private Map<String, Item> interactions = new HashMap();
    private Map<String, String> interactionitems = new HashMap();
    private Map<String, String> interactionterms = new HashMap<String, String>();
    private Map<String, Item> interactiontype = new HashMap();
    private Map<String, Item> interactiondetail = new HashMap();
    private Map<String, Item> experimenttype = new HashMap();
    private Map<String, Item> interactiondetectionmethods = new HashMap();
    private Map<String, Item> psiTerms = new HashMap();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public AllianceGeneticInteractionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception, ObjectStoreException {

        /*ID(s) interactor A     ID(s) interactor B      Alt. ID(s) interactor A Alt. ID(s) interactor B Alias(es)
        interactor A  Alias(es) interactor B  Interaction detection method(s) Publication 1st author(s)
        Publication Identifier(s)       Taxid interactor A      Taxid interactor B      Interaction type(s)
        Source database(s)      Interaction identifier(s)       Confidence value(s)     Expansion method(s)
        Biological role(s) interactor A Biological role(s) interactor B Experimental role(s) interactor A
        Experimental role(s) interactor B       Type(s) interactor A    Type(s) interactor B    Xref(s) interactor A
        Xref(s) interactor B    Interaction Xref(s)     Annotation(s) interactor A      Annotation(s) interactor B
        Interaction annotation(s)       Host organism(s)        Interaction parameter(s)
        Creation date   Update date     Checksum(s) interactor A        Checksum(s) interactor B
        Interaction Checksum(s) Negative        Feature(s) interactor A Feature(s)
         interactor B Stoichiometry(s) interactor A   Stoichiometry(s) interactor B
         Identification method participant A     Identification method participant B


        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int count = 0;
        System.out.println("Processing Interactions...");
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            count++;
            String geneFeatureName = res.getString("dbentity1_id");
            Item gene = genes.get(geneFeatureName);

            String interactionNo = res.getString("annotation_id");
            String referenceNo = res.getString("reference_id");
            String interactionType = "physical interactions";
            String experimentType = res.getString("biogrid_experimental_system");
            String annotationType = res.getString("annotation_type");
            String modification = res.getString("modification");

            String interactingGeneFeatureName = res.getString("dbentity2_id");
            Item interactingGene = genes.get(interactingGeneFeatureName);

            String action = res.getString("bait_hit");
            String[] a = action.split("-");
            String role1 = a[0];
            String role2 = a[1];
            String source = res.getString("display_name");
            String phenotype = ""; //res.getString("phenotype");
            String citation = res.getString("citation");
            String pubmed = res.getString("pmid");

            String interactionRefId = getInteraction(interactionNo,
                    referenceNo, interactionType, experimentType,
                    annotationType, modification, interactingGene, role1, source,
                    phenotype, citation, gene, pubmed);

            store the reverse relationship so that template changes do not have to be made
            if(!geneFeatureName.equals(interactingGeneFeatureName)) {
                String interactionRefId2 = getInteraction(interactionNo,
                        referenceNo, interactionType, experimentType,
                        annotationType, modification, gene, role2, source,
                        phenotype, citation, interactingGene, pubmed, title, volume, page,
                        year, issue, abbreviation, dsId, firstAuthor, dbxrefid, note);
            }
        }
        System.out.println("size of genetic interaction genes:  " + genes.size());
        storeInteractionTypes();
        storeInteractionExperiments();
        //storeInteractionDetails(); <--keep commented
        storeInteractions();
        storeGenes();
        */

    }


    /*
    private String getInteraction(String interactionNo, String referenceNo,
                                  String interactionType, String experimentType,
                                  String annotationType, String modification, Item interactingGene, String action,
                                  String source, String phenotype, String citation, Item gene,
                                  String pubMedId)
            throws ObjectStoreException {

        Item item = getInteractionItem(gene.getIdentifier(), interactingGene.getIdentifier());
        Item detail = createItem("InteractionDetail");

        detail.setAttribute("type", interactionType);
        detail.setAttribute("annotationType", annotationType);
        if (StringUtils.isNotEmpty(modification)) detail.setAttribute("modification", modification);
        if (StringUtils.isNotEmpty(phenotype)) detail.setAttribute("phenotype", phenotype);
        detail.setAttribute("role1", action);
        detail.addToCollection("allInteractors", interactingGene.getIdentifier());
        detail.addToCollection("dataSets", dsetIdentifier);

        String shortType = interactionType.substring(0, interactionType.indexOf(' '));
        detail.setAttribute("relationshipType", shortType); //interactionType
        String unqName = firstAuthor+"-"+pubMedId+"-"+experimentType;
        if (StringUtils.isNotEmpty(note)) detail.setAttribute("note", note);

        //add publication as experiment type
        Item storedExperimentType = experimenttype.get(unqName);
        if(storedExperimentType == null) {

            storedExperimentType = createItem("InteractionExperiment");
            storedExperimentType.setAttribute("name", unqName);
            experimenttype.put(unqName, storedExperimentType);

            String storedTermId = interactionterms.get(experimentType);
            if (storedTermId != null) {
                storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
            } else {
                storedTermId = getInteractionTerm(experimentType);
                storedExperimentType.addToCollection("interactionDetectionMethods", storedTermId );
            }
        }

        //add publication as reference on experiment
        Item storedRef = publications.get(referenceNo);

        if (storedRef != null) {
            storedExperimentType.setReference("publication", storedRef.getIdentifier());
        } else {

            Item pub = createItem("Publication");

            if (StringUtils.isNotEmpty(pubMedId)) {
                pub.setAttribute("pubMedId", pubMedId);
            }
            publications.put(referenceNo, pub);
            storedExperimentType.setReference("publication", pub.getIdentifier());
        }

        detail.setReference("experiment", storedExperimentType.getIdentifier());
        detail.setReference("interaction", item);

        try {
            store(detail);
        } catch (ObjectStoreException e) {
            throw new ObjectStoreException(e);
        }

        interactions.put(item.getIdentifier(), item);
        //interactionitems.put(interactionNo, item);
        String refId = item.getIdentifier();
        return refId;

    }



    private Item getInteractionItem(String refId, String gene2RefId) throws ObjectStoreException {
        MultiKey key = new MultiKey(refId, gene2RefId);
        Item interaction = interactionsnew.get(key);
        if (interaction == null) {
            interaction = createItem("Interaction");
            interaction.setReference("participant1", refId); //gene1
            interaction.setReference("participant2", gene2RefId); //gene2
        }
        return interaction;
    }


    private Item getGene(String g, String org) throws ObjectStoreException {

        Item gene  = genes.get(g);
        if(gene == null) {
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", g);
            gene.setReference("organism", org);
        }
        String geneId = gene.getIdentifier();
        genes.put(g, gene);
        return gene;
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


    private void storeInteractionTypes() throws ObjectStoreException {
        for (Item type : interactiontype.values()) {
            try {
                store(type);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }


    private void storeInteractionExperiments() throws ObjectStoreException {
        for (Item exp : experimenttype.values()) {
            try {
                store(exp);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }


    private void storeInteractionDetails() throws ObjectStoreException {
        for (Item det : interactiondetail.values()) {
            try {
                store(det);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }

    private void storeInteractions() throws ObjectStoreException {
        for (Item intact : interactions.values()) {
            try {
                store(intact);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }


    private void storeGenes() throws ObjectStoreException {
        for (Item gene : genes.values()) {
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
        }
    }
*/

}
