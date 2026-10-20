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

public class CgobHomologsConverterTest extends TestCase
{
    private int countHomologues(String tsv) throws Exception {
        MockItemWriter writer =
                new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        new CgobHomologsConverter(writer, Model.getInstanceByName("genomic"))
                .process(new StringReader(tsv));
        int homologues = 0;
        for (Item item : writer.getItems()) {
            if ("Homologue".equals(item.getClassName())) {
                homologues++;
            }
        }
        return homologues;
    }

    public void testPartialRowPairsOnlyThePresentGenes() throws Exception {
        // yeast gene with two of four Candida columns filled: 2 pairs, 2 rows each
        assertEquals(4, countHomologues("YAL001C\tCAL0001\t---\tCD36_0001\t---\n"));
    }

    public void testRowWithoutYeastGeneIsSkipped() throws Exception {
        assertEquals(0, countHomologues("---\tCAL0001\tCAL0002\t---\t---\n"));
    }

    public void testFullRowPairsAllFour() throws Exception {
        assertEquals(8, countHomologues("YAL001C\tCAL0001\tWO1_0001\tCD36_0001\tCPAR_0001\n"));
    }
}
