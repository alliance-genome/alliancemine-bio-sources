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
 * Reads transgenic-alleles.tsv (emitted by scripts/fetch_transgenic_alleles.py
 * from the Alliance /gene/{id}/transgenic-alleles endpoint) and produces
 * TransgenicAllele items linked to their target genes.
 */
public class AllianceTransgenicAllelesConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(AllianceTransgenicAllelesConverter.class);

    private static final String DATASET_TITLE = "Alliance Transgenic Alleles";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_transgenic_alleles.py COLUMNS
    private static final int COL_GENE_ID            = 0;
    private static final int COL_GENE_SYMBOL        = 1;
    private static final int COL_GENE_TAXON         = 2;
    private static final int COL_ALLELE_ID          = 3;
    private static final int COL_ALLELE_SYMBOL      = 4;
    private static final int COL_CONSTRUCT_IDS      = 5;
    private static final int COL_DATA_PROVIDER      = 6;
    private static final int MIN_COLUMNS            = 7;

    private final Map<String, String> genes = new HashMap<String, String>();
    private final Map<String, String> alleles = new HashMap<String, String>();

    public AllianceTransgenicAllelesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance transgenic alleles...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("geneId".equals(line[COL_GENE_ID])) {
                continue;  // TSV header
            }
            processRow(line);
            rows++;
        }
        LOG.info("Transgenic alleles: processed " + rows + " rows");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String geneId = line[COL_GENE_ID].trim();
        String alleleId = line[COL_ALLELE_ID].trim();
        if (alleleId.isEmpty()) {
            return;
        }

        String organismRef = null;
        String taxon = line[COL_GENE_TAXON].trim();
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            organismRef = getOrganism(taxonId);
        }

        String allele = alleles.get(alleleId);
        if (allele == null) {
            Item item = createItem("TransgenicAllele");
            item.setAttribute("primaryIdentifier", alleleId);
            setIfPresent(item, "alleleSymbol", line[COL_ALLELE_SYMBOL]);
            setIfPresent(item, "constructIdentifier", line[COL_CONSTRUCT_IDS]);
            if (organismRef != null) {
                item.setReference("organism", organismRef);
            }
            if (!geneId.isEmpty()) {
                String geneRef = getGene(geneId, organismRef);
                // TransgenicAllele extends Allele which has reference "gene".
                item.setReference("gene", geneRef);
            }
            store(item);
            alleles.put(alleleId, item.getIdentifier());
        }
    }

    private String getGene(String primaryId, String organismRef) throws ObjectStoreException {
        String ref = genes.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item gene = createItem("Gene");
        gene.setAttribute("primaryIdentifier", primaryId);
        if (organismRef != null) {
            gene.setReference("organism", organismRef);
        }
        store(gene);
        genes.put(primaryId, gene.getIdentifier());
        return gene.getIdentifier();
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
