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
 * Reads gene-descriptions.tsv (emitted by
 * scripts/fetch_gene_descriptions.py from the Alliance FMS
 * GENE-DESCRIPTION-JSON bulk files) and produces partial Gene items
 * carrying the AGR auto-curated description fields. Merged with sibling
 * Gene items via Gene.key_primaryidentifier.
 */
public class AllianceGeneDescriptionsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceGeneDescriptionsConverter.class);

    private static final String DATASET_TITLE = "Alliance Curated Gene Descriptions";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_gene_descriptions.py COLUMNS
    private static final int COL_GENE_ID       = 0;
    private static final int COL_GENE_NAME     = 1;
    private static final int COL_AUTO          = 2;
    private static final int COL_GO            = 3;
    private static final int COL_GO_F          = 4;
    private static final int COL_GO_P          = 5;
    private static final int COL_GO_C          = 6;
    private static final int COL_DO            = 7;
    private static final int COL_EXPRESSION    = 8;
    private static final int COL_ORTHOLOGY     = 9;
    private static final int COL_DATA_PROVIDER = 10;
    private static final int MIN_COLUMNS       = 11;

    private final Map<String, String> genes = new HashMap<String, String>();

    public AllianceGeneDescriptionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance gene descriptions...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("geneId".equals(line[COL_GENE_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("GeneDescriptions: emitted " + rows + " partial Gene items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String geneId = line[COL_GENE_ID].trim();
        if (geneId.isEmpty()) {
            return;
        }
        Item gene = createItem("Gene");
        gene.setAttribute("primaryIdentifier", geneId);
        setIfPresent(gene, "autoDescription",        line[COL_AUTO]);
        setIfPresent(gene, "goDescription",          line[COL_GO]);
        setIfPresent(gene, "goFunctionDescription",  line[COL_GO_F]);
        setIfPresent(gene, "goProcessDescription",   line[COL_GO_P]);
        setIfPresent(gene, "goComponentDescription", line[COL_GO_C]);
        setIfPresent(gene, "doDescription",          line[COL_DO]);
        setIfPresent(gene, "expressionDescription",  line[COL_EXPRESSION]);
        setIfPresent(gene, "orthologyDescription",   line[COL_ORTHOLOGY]);
        store(gene);
        genes.put(geneId, gene.getIdentifier());
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
