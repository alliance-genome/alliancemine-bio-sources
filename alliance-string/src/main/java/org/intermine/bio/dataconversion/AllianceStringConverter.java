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
 * Reads string-interactions.tsv (emitted by scripts/fetch_string.py
 * from STRING DB v12 per-organism downloads) and produces
 * STRINGInteraction items with protein references resolved via UniProt
 * accession, plus organism reference.
 */
public class AllianceStringConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceStringConverter.class);

    private static final String DATASET_TITLE = "STRING v12 functional protein-protein associations";
    private static final String DATA_SOURCE_NAME = "STRING";

    private static final int COL_STRING_A = 0;
    private static final int COL_STRING_B = 1;
    private static final int COL_UP_A     = 2;
    private static final int COL_UP_B     = 3;
    private static final int COL_SCORE    = 4;
    private static final int COL_TAXON    = 5;
    private static final int MIN_COLUMNS  = 6;

    private final Map<String, String> proteinsByAcc = new HashMap<String, String>();

    public AllianceStringConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing STRING interactions...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("stringIdA".equals(line[COL_STRING_A])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("STRING: emitted " + rows + " STRINGInteraction items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String sa = line[COL_STRING_A].trim();
        String sb = line[COL_STRING_B].trim();
        if (sa.isEmpty() || sb.isEmpty()) {
            return;
        }
        Item it = createItem("STRINGInteraction");
        it.setAttribute("stringIdA", sa);
        it.setAttribute("stringIdB", sb);
        String score = line[COL_SCORE].trim();
        if (!score.isEmpty()) {
            it.setAttribute("combinedScore", score);
        }
        String upA = line[COL_UP_A].trim();
        if (!upA.isEmpty()) {
            it.setReference("proteinA", getProtein(upA));
        }
        String upB = line[COL_UP_B].trim();
        if (!upB.isEmpty()) {
            it.setReference("proteinB", getProtein(upB));
        }
        String taxon = line[COL_TAXON].trim();
        if (StringUtils.isNotEmpty(taxon)) {
            String tid = taxon.contains(":") ? taxon.substring(taxon.indexOf(':') + 1) : taxon;
            it.setReference("organism", getOrganism(tid));
        }
        store(it);
    }

    private String getProtein(String acc) throws ObjectStoreException {
        String ref = proteinsByAcc.get(acc);
        if (ref != null) {
            return ref;
        }
        Item p = createItem("Protein");
        p.setAttribute("primaryAccession", acc);
        store(p);
        proteinsByAcc.put(acc, p.getIdentifier());
        return p.getIdentifier();
    }
}
