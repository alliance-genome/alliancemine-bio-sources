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
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Reads chromosome-lengths.tsv (emitted by scripts/fetch_chromosome_lengths.py from
 * NCBI Datasets) and stores partial Chromosome items with a length. They merge with the
 * name-only Chromosomes from alliance-genes on Chromosome.key_primaryidentifier
 * (primaryIdentifier, organism), so create-intergenic-region-features and
 * create-gene-flanking-features, which skip chromosomes without a length, cover every MOD.
 */
public class AllianceChromosomeLengthsConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(AllianceChromosomeLengthsConverter.class);

    private static final String DATASET_TITLE = "NCBI chromosome lengths";
    private static final String DATA_SOURCE_NAME = "NCBI";

    // Must match scripts/fetch_chromosome_lengths.py COLUMNS
    private static final int COL_TAXON = 0;
    private static final int COL_CHROMOSOME = 1;
    private static final int COL_LENGTH = 2;
    private static final int MIN_COLUMNS = 3;

    public AllianceChromosomeLengthsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(Reader reader) throws Exception {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int rows = 0;
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            if (line.length < MIN_COLUMNS || "taxonId".equals(line[COL_TAXON])) {
                continue;
            }
            Item chromosome = createItem("Chromosome");
            chromosome.setAttribute("primaryIdentifier", line[COL_CHROMOSOME].trim());
            chromosome.setAttribute("length", line[COL_LENGTH].trim());
            chromosome.setReference("organism", getOrganism(line[COL_TAXON].trim()));
            store(chromosome);
            rows++;
        }
        LOG.info("ChromosomeLengths: stored " + rows + " Chromosome lengths");
    }
}
