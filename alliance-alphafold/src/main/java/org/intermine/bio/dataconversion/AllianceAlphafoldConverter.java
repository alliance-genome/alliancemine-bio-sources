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
import java.util.Iterator;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads alphafold.tsv (emitted by scripts/fetch_alphafold.py from
 * gene-crossrefs UniProtKB rows) and produces partial Protein items
 * carrying alphaFoldId + alphaFoldUrl, keyed on primaryAccession so the
 * integration engine merges into existing UniProt-source Protein rows.
 */
public class AllianceAlphafoldConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceAlphafoldConverter.class);

    private static final String DATASET_TITLE = "AlphaFold structure cross-references";
    private static final String DATA_SOURCE_NAME = "EBI AlphaFold";

    // Must match scripts/fetch_alphafold.py COLUMNS
    private static final int COL_ACC      = 0;
    private static final int COL_GENE_ID  = 1;
    private static final int COL_TAXON    = 2;
    private static final int COL_AF_URL   = 3;
    private static final int MIN_COLUMNS  = 4;

    public AllianceAlphafoldConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        LOG.info("Processing AlphaFold structure cross-references...");
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS) {
                continue;
            }
            if ("uniProtAccession".equals(line[COL_ACC])) {
                continue;
            }
            processRow(line);
            rows++;
        }
        LOG.info("AlphaFold: emitted " + rows + " partial Protein items");
    }

    private void processRow(String[] line) throws ObjectStoreException {
        String acc = line[COL_ACC].trim();
        if (acc.isEmpty()) {
            return;
        }
        Item protein = createItem("Protein");
        protein.setAttribute("primaryAccession", acc);
        protein.setAttribute("alphaFoldId", acc);
        protein.setAttribute("alphaFoldUrl", line[COL_AF_URL].trim());
        store(protein);
    }
}
