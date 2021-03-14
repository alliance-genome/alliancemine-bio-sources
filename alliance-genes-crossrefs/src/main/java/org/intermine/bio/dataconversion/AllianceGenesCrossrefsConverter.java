package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2009 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */
import org.apache.commons.collections.keyvalue.MultiKey;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.model.bio.BioEntity;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;
import org.intermine.util.FormattedTextParser;

import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.io.Reader;

/**
 *
 * @author
 */
public class AllianceGenesCrossrefsConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(AllianceGenesCrossrefsConverter.class);
    private static final String DATASET_TITLE = "Alliance CrossRef data";
    private static final String DATA_SOURCE_NAME = "Alliance CrossRef data";
    private Map<String, Item> genes = new HashMap();
    private Map<String, Item> organisms = new LinkedHashMap<String, Item>();

    /**
     * Construct a new AllianceGenesCrossrefsConverter.
     * @param database the database to read from
     * @param model the Model used by the object store we will write to with the ItemWriter
     * @param writer an ItemWriter used to handle Items created
     */
    public AllianceGenesCrossrefsConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

        /*
       GeneID  GlobalCrossReferenceID  CrossReferenceCompleteURL       ResourceDescriptorPage  TaxonID
         */
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int count = 0;
        System.out.println("Processing GeneCrossRefCombined.tsv...");
        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();

            if(count <= 14) { count++; continue;}

            String geneId = line[0].trim();
            String xrefId = line[1].trim();

            if(geneId.equals(xrefId)) { continue; }

            String origspecies = line[4].trim();
            String taxonId = origspecies.replace("NCBITaxon:", "");
            if(taxonId.equalsIgnoreCase("Taxon")){ continue; }
            Item organism = newOrganism(taxonId);

           String refId = getGene(geneId, organism);

            if (refId != null) {
                getCrossReference(refId, xrefId);
            }

        }
        storeGenes();

    }

    private String getGene(String geneId, Item org) throws ObjectStoreException {

        Item gene = genes.get(geneId);
        if (gene == null) {
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", geneId);
            gene.setReference("organism", org);
            //store(gene);
            genes.put(geneId, gene);
        }
        return gene.getIdentifier();
    }

    /**
     *
     * @param subjectId
     * @param id
     * @throws ObjectStoreException
     */
    private String getCrossReference(String subjectId, String id)
            throws ObjectStoreException {

        String refId = "";
        String type = "";
        if(id.contains(":")) {
            String[] t = id.split(":");
            type = t[1];
        }else{
            type = id;
        }
        Item crf = createItem("CrossReference");
        crf.setReference("subject", subjectId);
        crf.setAttribute("identifier", id);
        if(!type.equalsIgnoreCase(id)) { crf.setAttribute("dbxreftype", type);}

        try {
            store(crf);
        } catch (ObjectStoreException e) {
            throw new ObjectStoreException(e);
        }

        refId = crf.getIdentifier();
        return refId;

    }

    private Item newOrganism(String taxonId) throws ObjectStoreException {
        Item item = organisms.get(taxonId);
        if (item == null) {
            item = createItem("Organism");
            item.setAttribute("taxonId", taxonId);
            organisms.put(taxonId, item);
            store(item);
        }
        return item;
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


}
