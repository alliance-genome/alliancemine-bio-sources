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
 * Reads gene-crossrefs.tsv (emitted by scripts/fetch_gene_crossrefs.py
 * from the Alliance FMS GENECROSSREFERENCEJSON bulk file) and produces
 * CrossReference items linked to existing Gene rows. The external DB
 * is encoded in the xref id prefix (UniProtKB:..., Ensembl:..., RefSeq:...,
 * NCBI_Gene:..., etc.) which the converter splits into source + accession.
 */
public class AllianceGeneCrossrefsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceGeneCrossrefsConverter.class);

    private static final String DATASET_TITLE = "Alliance Gene Cross-references";
    private static final String DATA_SOURCE_NAME = "Alliance of Genome Resources";

    // Must match scripts/fetch_gene_crossrefs.py COLUMNS
    private static final int COL_GENE_ID    = 0;
    private static final int COL_GENE_TAXON = 1;
    private static final int COL_XREF_ID    = 2;
    private static final int COL_XREF_URL   = 3;
    private static final int COL_XREF_TYPE  = 4;
    private static final int MIN_COLUMNS    = 5;

    private final Map<String, String> genes = new HashMap<String, String>();
    private final Map<String, String> dataSources = new HashMap<String, String>();

    public AllianceGeneCrossrefsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing Alliance gene cross-references...");
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
        LOG.info("GeneCrossrefs: emitted " + rows + " CrossReference items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String geneId = line[COL_GENE_ID].trim();
        String xrefId = line[COL_XREF_ID].trim();
        if (geneId.isEmpty() || xrefId.isEmpty()) {
            return;
        }
        String xrefSource;
        int colon = xrefId.indexOf(':');
        if (colon > 0) {
            xrefSource = xrefId.substring(0, colon);
        } else {
            xrefSource = line[COL_XREF_TYPE].trim();
        }

        Item xref = createItem("CrossReference");
        xref.setAttribute("identifier", xrefId);
        setIfPresent(xref, "url", line[COL_XREF_URL]);
        setIfPresent(xref, "dbxreftype", line[COL_XREF_TYPE]);
        xref.setReference("subject", getGene(geneId, line[COL_GENE_TAXON].trim()));
        if (!xrefSource.isEmpty()) {
            xref.setReference("source", getDataSourceRef(xrefSource));
        }
        store(xref);
    }

    private String getGene(String primaryId, String taxon) throws ObjectStoreException {
        String ref = genes.get(primaryId);
        if (ref != null) {
            return ref;
        }
        Item item = createItem("Gene");
        item.setAttribute("primaryIdentifier", primaryId);
        if (StringUtils.isNotEmpty(taxon)) {
            String taxonId = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            item.setReference("organism", getOrganism(taxonId));
        }
        store(item);
        genes.put(primaryId, item.getIdentifier());
        return item.getIdentifier();
    }

    private String getDataSourceRef(String name) throws ObjectStoreException {
        String ref = dataSources.get(name);
        if (ref != null) {
            return ref;
        }
        Item ds = createItem("DataSource");
        ds.setAttribute("name", name);
        store(ds);
        dataSources.put(name, ds.getIdentifier());
        return ds.getIdentifier();
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
