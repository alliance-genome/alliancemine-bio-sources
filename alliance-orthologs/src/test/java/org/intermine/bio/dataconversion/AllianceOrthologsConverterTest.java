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

public class AllianceOrthologsConverterTest extends TestCase
{
    public void testHomologueTypeIsOrthologue() throws Exception {
        StringBuilder tsv = new StringBuilder();
        for (int i = 0; i < 16; i++) {
            tsv.append("header line ").append(i).append('\n');
        }
        tsv.append("FB:FBgn0000001\tsym1\tNCBITaxon:7227\tDrosophila melanogaster"
                + "\tSGD:S000000001\tsym2\tNCBITaxon:559292\tSaccharomyces cerevisiae"
                + "\tDIOPT\t5\t10\tYes\tNo\n");
        MockItemWriter writer =
                new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        new AllianceOrthologsConverter(writer, Model.getInstanceByName("genomic"))
                .process(new StringReader(tsv.toString()));
        int homologues = 0;
        for (Item item : writer.getItems()) {
            if ("Homologue".equals(item.getClassName())) {
                homologues++;
                assertEquals("orthologue", item.getAttribute("type").getValue());
            }
        }
        assertEquals(1, homologues);
    }
}
