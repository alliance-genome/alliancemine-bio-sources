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
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads worm-rnai.tsv (emitted by scripts/fetch_wormmine_rnai.py from
 * WormMine via PathQuery REST) and produces RNAi items linked to a Gene
 * (the inhibited gene), Strain, and Publication.
 *
 * @author
 */
public class AllianceWormRnaiConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceWormRnaiConverter.class);

    private static final String DATASET_TITLE = "WormMine RNAi";
    private static final String DATA_SOURCE_NAME = "WormMine";
    private static final String WORM_TAXON = "6239";

    // Must match scripts/fetch_wormmine_rnai.py COLUMNS
    private static final int COL_RNAI_ID          = 0;
    private static final int COL_METHOD           = 1;
    private static final int COL_TREATMENT        = 2;
    private static final int COL_TEMPERATURE      = 3;
    private static final int COL_GENOTYPE         = 4;
    private static final int COL_DELIVERED_BY     = 5;
    private static final int COL_PHENOTYPE_REMARK = 6;
    private static final int COL_REMARK           = 7;
    private static final int COL_REFERENCE_ID     = 8;
    private static final int COL_STRAIN_ID        = 9;
    private static final int COL_GENE_ID          = 10;
    private static final int MIN_COLUMNS          = 11;

    private final Map<String, String> genes = new HashMap<String, String>();
    private final Map<String, String> strains = new HashMap<String, String>();
    private final Map<String, String> publications = new HashMap<String, String>();
    // RNAi items are accumulated in this map keyed by primaryIdentifier so that
    // the (RNAi, gene) row pairs returned by the WormMine query collapse into a
    // single RNAi Item with a multi-gene `genes` collection.
    private final Map<String, Item> rnaiItems = new HashMap<String, Item>();
    private String organismRef = null;

    public AllianceWormRnaiConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        LOG.info("Processing WormMine RNAi data...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("rnaiId".equals(line[COL_RNAI_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        // Flush all deduplicated RNAi items at the end. Once stored, an Item
        // can't be modified further, so we accumulate first and store last.
        for (Item r : rnaiItems.values()) {
            store(r);
        }
        LOG.info("WormRNAi: processed " + rows + " (RNAi, gene) rows -> "
                + rnaiItems.size() + " distinct RNAi items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String rnaiId = line[COL_RNAI_ID].trim();
        if (rnaiId.isEmpty()) {
            return;
        }
        Item r = rnaiItems.get(rnaiId);
        if (r == null) {
            r = createItem("RNAi");
            r.setAttribute("primaryIdentifier", rnaiId);
            setIfPresent(r, "method", line[COL_METHOD]);
            setIfPresent(r, "treatment", line[COL_TREATMENT]);
            setIfPresent(r, "temperature", line[COL_TEMPERATURE]);
            setIfPresent(r, "genotype", line[COL_GENOTYPE]);
            setIfPresent(r, "deliveredBy", line[COL_DELIVERED_BY]);
            setIfPresent(r, "phenotypeRemark", line[COL_PHENOTYPE_REMARK]);
            setIfPresent(r, "remark", line[COL_REMARK]);
            r.setReference("organism", getWormOrganism());

            String strainId = line[COL_STRAIN_ID].trim();
            if (!strainId.isEmpty()) {
                r.setReference("strain", getStrain(strainId));
            }
            String pubId = line[COL_REFERENCE_ID].trim();
            if (!pubId.isEmpty()) {
                String pubRef = getPublication(pubId);
                if (pubRef != null) {
                    r.setReference("reference", pubRef);
                }
            }
            rnaiItems.put(rnaiId, r);
        }
        // gene is the per-row varying column; everything else is set on first
        // encounter. Multiple inhibitsGene rows for the same RNAi pile up here.
        String geneId = line[COL_GENE_ID].trim();
        if (!geneId.isEmpty()) {
            r.addToCollection("genes", getGene(geneId));
        }
    }

    private String getWormOrganism() throws ObjectStoreException {
        if (organismRef != null) {
            return organismRef;
        }
        organismRef = getOrganism(WORM_TAXON);
        return organismRef;
    }

    private String getGene(String primaryId) throws ObjectStoreException {
        String ref = genes.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item gene = createItem("Gene");
        gene.setAttribute("primaryIdentifier", primaryId);
        gene.setReference("organism", getWormOrganism());
        store(gene);
        genes.put(primaryId, gene.getIdentifier());
        return gene.getIdentifier();
    }

    private String getStrain(String primaryId) throws ObjectStoreException {
        String ref = strains.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item s = createItem("Strain");
        s.setAttribute("primaryIdentifier", primaryId);
        s.setReference("organism", getWormOrganism());
        store(s);
        strains.put(primaryId, s.getIdentifier());
        return s.getIdentifier();
    }

    private String getPublication(String pubToken) throws ObjectStoreException {
        if (pubToken == null || pubToken.trim().isEmpty()) {
            return null;
        }
        String p = pubToken.trim();
        String ref = publications.get(p);
        if (ref != null) {
            return ref;
        }
        Item pub = createItem("Publication");
        if (p.startsWith("PMID:")) {
            pub.setAttribute("pubMedId", p.substring("PMID:".length()));
        } else {
            // WormMine returns WBPaper IDs not PMIDs - store via pubXrefId.
            pub.setAttribute("pubXrefId", p);
        }
        store(pub);
        publications.put(p, pub.getIdentifier());
        return pub.getIdentifier();
    }

    private static void setIfPresent(Item item, String attrName, String value) {
        if (value != null) {
            String v = value.trim();
            if (!v.isEmpty() && !"-".equals(v)) {
                item.setAttribute(attrName, v);
            }
        }
    }
}
