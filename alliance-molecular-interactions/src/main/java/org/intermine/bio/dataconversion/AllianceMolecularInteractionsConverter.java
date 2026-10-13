package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2026 AllianceMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads the TSV emitted by scripts/fetch_interactions.py (molecular-interactions.tsv)
 * and produces Interaction + InteractionDetail items.
 *
 * TSV columns (from fetch_interactions.MOLECULAR_COLUMNS):
 *   0  geneId                      13  interactorBRoleName
 *   1  geneSymbol                   14  interactorARole       <-- dup check
 *   2  geneTaxon                    ...
 *  and so on, see COL_* constants.
 *
 * @author
 */
public class AllianceMolecularInteractionsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceMolecularInteractionsConverter.class);

    private static final String DATASET_TITLE = "Alliance Molecular Interactions";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // TSV column positions - must match scripts/fetch_interactions.py MOLECULAR_COLUMNS
    private static final int COL_GENE_ID                    = 0;
    private static final int COL_GENE_SYMBOL                = 1;
    private static final int COL_GENE_TAXON                 = 2;
    private static final int COL_PARTNER_GENE_ID            = 3;
    private static final int COL_PARTNER_SYMBOL             = 4;
    private static final int COL_PARTNER_TAXON              = 5;
    private static final int COL_INTERACTION_TYPE           = 6;
    private static final int COL_INTERACTION_TYPE_NAME      = 7;
    private static final int COL_INTERACTOR_A_TYPE          = 8;
    private static final int COL_INTERACTOR_A_TYPE_NAME     = 9;
    private static final int COL_INTERACTOR_B_TYPE          = 10;
    private static final int COL_INTERACTOR_B_TYPE_NAME     = 11;
    private static final int COL_INTERACTOR_A_ROLE          = 12;
    private static final int COL_INTERACTOR_A_ROLE_NAME     = 13;
    private static final int COL_INTERACTOR_B_ROLE          = 14;
    private static final int COL_INTERACTOR_B_ROLE_NAME     = 15;
    private static final int COL_INTERACTION_SOURCE         = 16;
    private static final int COL_INTERACTION_SOURCE_NAME    = 17;
    private static final int COL_RELATION                   = 18;
    private static final int COL_INTERACTION_ID             = 19;
    private static final int COL_UNIQUE_ID                  = 20;
    private static final int COL_PUBMED_ID                  = 21;
    private static final int COL_REFERENCE_ID               = 22;
    private static final int COL_SHORT_CITATION             = 23;
    private static final int COL_CROSSREFS                  = 24;
    private static final int COL_DETECTION_METHOD           = 25;
    private static final int COL_DETECTION_METHOD_NAME      = 26;
    private static final int COL_AGGREGATION_DATABASE       = 27;
    private static final int COL_AGGREGATION_DATABASE_NAME  = 28;
    private static final int MIN_COLUMNS                    = 29;

    private final Map<String, String> genes = new HashMap<String, String>();
    private final Map<String, String> publications = new HashMap<String, String>();
    private final Map<String, String> terms = new HashMap<String, String>();
    private final Set<String> storedInteractionUids = new HashSet<String>();

    public AllianceMolecularInteractionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        LOG.info("Processing molecular interactions...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        int skippedHeader = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            String first = line[COL_GENE_ID];
            // Skip the TSV header line itself.
            if ("geneId".equals(first)) {
                skippedHeader++;
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("Molecular: processed " + rows + " rows, skipped " + skippedHeader + " header row(s), unique interactions=" + storedInteractionUids.size());
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String uniqueId = line[COL_UNIQUE_ID].trim();
        if (uniqueId.isEmpty() || storedInteractionUids.contains(uniqueId)) {
            return;
        }
        storedInteractionUids.add(uniqueId);

        String geneId = line[COL_GENE_ID].trim();
        String partnerGeneId = line[COL_PARTNER_GENE_ID].trim();
        if (geneId.isEmpty() || partnerGeneId.isEmpty()) {
            return;
        }

        String geneRef = getGene(geneId, line[COL_GENE_TAXON].trim());
        String partnerRef = getGene(partnerGeneId, line[COL_PARTNER_TAXON].trim());

        Item interaction = createItem("Interaction");
        interaction.setReference("participant1", geneRef);
        interaction.setReference("participant2", partnerRef);
        setIfPresent(interaction, "relation", line[COL_RELATION]);
        setIfPresent(interaction, "aggregationDatabase", line[COL_AGGREGATION_DATABASE]);
        setReferenceIfPresent(interaction, "detectionMethod",
                getInteractionTerm(line[COL_DETECTION_METHOD], line[COL_DETECTION_METHOD_NAME]));
        setReferenceIfPresent(interaction, "interactionSource",
                getInteractionTerm(line[COL_INTERACTION_SOURCE], line[COL_INTERACTION_SOURCE_NAME]));
        store(interaction);

        Item detail = createItem("InteractionDetail");
        detail.setReference("interaction", interaction);
        setIfPresent(detail, "type", line[COL_INTERACTION_TYPE]);
        setIfPresent(detail, "shortName", line[COL_INTERACTION_TYPE_NAME]);
        setIfPresent(detail, "name", uniqueId);
        setIfPresent(detail, "relationshipType", line[COL_RELATION]);
        setIfPresent(detail, "role1", line[COL_INTERACTOR_A_ROLE]);
        setIfPresent(detail, "role1Name", line[COL_INTERACTOR_A_ROLE_NAME]);
        setIfPresent(detail, "role2", line[COL_INTERACTOR_B_ROLE]);
        setIfPresent(detail, "role2Name", line[COL_INTERACTOR_B_ROLE_NAME]);
        setReferenceIfPresent(detail, "participant1Type",
                getInteractionTerm(line[COL_INTERACTOR_A_TYPE], line[COL_INTERACTOR_A_TYPE_NAME]));
        setReferenceIfPresent(detail, "participant2Type",
                getInteractionTerm(line[COL_INTERACTOR_B_TYPE], line[COL_INTERACTOR_B_TYPE_NAME]));
        String pubRef = getPublication(line[COL_PUBMED_ID], line[COL_REFERENCE_ID], line[COL_SHORT_CITATION]);
        if (pubRef != null) {
            Item experiment = createItem("InteractionExperiment");
            experiment.setReference("publication", pubRef);
            if (StringUtils.isNotEmpty(line[COL_DETECTION_METHOD].trim())) {
                String termRef = getInteractionTerm(line[COL_DETECTION_METHOD], line[COL_DETECTION_METHOD_NAME]);
                if (termRef != null) {
                    experiment.addToCollection("interactionDetectionMethods", termRef);
                }
            }
            experiment.setAttribute("name", line[COL_SHORT_CITATION].trim());
            store(experiment);
            detail.setReference("experiment", experiment);
        }
        store(detail);
    }

    private String getGene(String primaryId, String taxon) throws ObjectStoreException {
        String ref = genes.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item gene = createItem("Gene");
        gene.setAttribute("primaryIdentifier", primaryId);
        if (StringUtils.isNotEmpty(taxon)) {
            // taxonCurie is e.g. "NCBITaxon:559292"; parent class expects the numeric taxon id.
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            gene.setReference("organism", getOrganism(taxonId));
        }
        store(gene);
        genes.put(primaryId, gene.getIdentifier());
        return gene.getIdentifier();
    }

    private String getPublication(String pubmedId, String referenceId, String shortCitation) throws ObjectStoreException {
        String p = StringUtils.isNotEmpty(pubmedId) ? pubmedId.trim() : referenceId.trim();
        if (p.isEmpty()) {
            return null;
        }
        String ref = publications.get(p);
        if (ref != null) {
            return ref;
        }
        Item pub = createItem("Publication");
        if (p.startsWith("PMID:")) {
            pub.setAttribute("pubMedId", p.substring("PMID:".length()));
        } else {
            pub.setAttribute("pubXrefId", p);
        }
        if (StringUtils.isNotEmpty(shortCitation)) {
            pub.setAttribute("citation", shortCitation.trim());
        }
        store(pub);
        publications.put(p, pub.getIdentifier());
        return pub.getIdentifier();
    }

    private String getInteractionTerm(String curie, String name) throws ObjectStoreException {
        String c = curie == null ? "" : curie.trim();
        if (c.isEmpty()) {
            return null;
        }
        String ref = terms.get(c);
        if (ref != null) {
            return ref;
        }
        Item term = createItem("InteractionTerm");
        term.setAttribute("identifier", c);
        if (StringUtils.isNotEmpty(name)) {
            term.setAttribute("name", name.trim());
        }
        store(term);
        terms.put(c, term.getIdentifier());
        return term.getIdentifier();
    }

    private static void setIfPresent(Item item, String attrName, String value) {
        if (value != null) {
            String v = value.trim();
            if (!v.isEmpty() && !"-".equals(v)) {
                item.setAttribute(attrName, v);
            }
        }
    }

    private static void setReferenceIfPresent(Item item, String refName, String refId) {
        if (refId != null && !refId.isEmpty()) {
            item.setReference(refName, refId);
        }
    }
}
