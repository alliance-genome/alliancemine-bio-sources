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

import java.io.StringReader;
import java.util.HashMap;

import junit.framework.TestCase;

import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

public class AllianceDiseaseDetailConverterTest extends TestCase
{
    private static final String HEADER = "diseaseId\tdiseaseName\tgeneId\tgeneTaxon\tgenRel"
            + "\tqualifiers\teco\tslim\tvia\tpmids\trefs\tuniqueId\n";

    private MockItemWriter writer;
    private AllianceDiseaseDetailConverter converter;

    public void setUp() throws Exception {
        writer = new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        converter = new AllianceDiseaseDetailConverter(writer, Model.getInstanceByName("genomic"));
    }

    private static String row(int n, String taxon) {
        return "DOID:" + n + "\tdisease" + n + "\tSGD:S" + n + "\t" + taxon
                + "\tis_implicated_in\t\t\t\t\t\t\tu" + n + "\n";
    }

    private int countGenes() {
        int genes = 0;
        for (Item item : writer.getItems()) {
            if ("Gene".equals(item.getClassName())) {
                genes++;
            }
        }
        return genes;
    }

    public void testRowWithoutTaxonIsSkippedBelowThreshold() throws Exception {
        StringBuilder tsv = new StringBuilder(HEADER);
        for (int i = 0; i < 200; i++) {
            tsv.append(row(i, "NCBITaxon:559292"));
        }
        tsv.append(row(900, ""));
        converter.process(new StringReader(tsv.toString()));
        assertEquals(200, countGenes());
    }

    public void testMissingTaxonAboveOnePercentThrows() throws Exception {
        StringBuilder tsv = new StringBuilder(HEADER);
        for (int i = 0; i < 98; i++) {
            tsv.append(row(i, "NCBITaxon:559292"));
        }
        tsv.append(row(900, ""));
        tsv.append(row(901, ""));
        try {
            converter.process(new StringReader(tsv.toString()));
            fail("2 of 100 rows without taxon must throw");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().startsWith("2 of 100"));
        }
    }
}
