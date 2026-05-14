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
 * Reads crispr-screens.tsv (emitted by
 * scripts/fetch_crispr_screens.py from the Alliance FMS BIOGRID-ORCS
 * tarballs, filtered to HIT=YES rows) and produces CRISPRScreenResult
 * items keyed on screenId + geneSymbol. Gene reference uses symbol +
 * organism lookup since BIOGRID-ORCS rows do not carry MOD CURIEs.
 */
public class AllianceCrisprScreensConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceCrisprScreensConverter.class);

    private static final String DATASET_TITLE = "BioGRID ORCS CRISPR/RNAi Screens";
    private static final String DATA_SOURCE_NAME = "BioGRID";

    private static final int COL_SCREEN_ID = 0;
    private static final int COL_GENE_ID   = 1;
    private static final int COL_GENE_SYM  = 2;
    private static final int COL_SCORE1    = 3;
    private static final int COL_SCORE2    = 4;
    private static final int COL_HIT       = 5;
    private static final int COL_TAXON     = 6;
    private static final int COL_SOURCE    = 7;
    private static final int COL_PROVIDER  = 8;
    private static final int MIN_COLUMNS   = 9;

    private final Map<String, String> genes = new HashMap<String, String>();

    public AllianceCrisprScreensConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing BioGRID ORCS screen hits...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("screenId".equals(line[COL_SCREEN_ID])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("CRISPRScreens: emitted " + rows + " CRISPRScreenResult items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String screen = line[COL_SCREEN_ID].trim();
        String sym = line[COL_GENE_SYM].trim();
        if (screen.isEmpty() || sym.isEmpty()) {
            return;
        }
        Item r = createItem("CRISPRScreenResult");
        r.setAttribute("screenId", screen);
        r.setAttribute("geneSymbol", sym);
        setIfPresent(r, "score1", line[COL_SCORE1]);
        setIfPresent(r, "score2", line[COL_SCORE2]);
        String hit = line[COL_HIT].trim().toLowerCase();
        if ("yes".equals(hit)) {
            r.setAttribute("hit", "true");
        } else if ("no".equals(hit)) {
            r.setAttribute("hit", "false");
        }
        setIfPresent(r, "source", line[COL_SOURCE]);

        String taxon = line[COL_TAXON].trim();
        if (StringUtils.isNotEmpty(taxon)) {
            String tid = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            String orgRef = getOrganism(tid);
            r.setReference("organism", orgRef);
            r.setReference("gene", getGene(sym, tid));
        }
        store(r);
    }

    private String getGene(String symbol, String taxonId) throws ObjectStoreException {
        String key = symbol + "|" + taxonId;
        String ref = genes.get(key);
        if (ref != null) {
            return ref;
        }
        Item gene = createItem("Gene");
        gene.setAttribute("symbol", symbol);
        gene.setReference("organism", getOrganism(taxonId));
        store(gene);
        genes.put(key, gene.getIdentifier());
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
