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
import java.util.HashSet;
import java.util.Set;

import junit.framework.TestCase;

import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

public class AllianceWormRnaiConverterTest extends TestCase
{
    private AllianceWormRnaiConverter converter;
    private MockItemWriter itemWriter;

    public void setUp() throws Exception {
        Model model = Model.getInstanceByName("genomic");
        itemWriter = new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        converter = new AllianceWormRnaiConverter(itemWriter, model);
    }

    private static String row(String rnai, String reference, String strain, String gene) {
        String[] f = new String[11];
        java.util.Arrays.fill(f, "");
        f[0] = rnai;
        f[8] = reference;
        f[9] = strain;
        f[10] = gene;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < f.length; i++) {
            sb.append(i == 0 ? "" : "\t").append(f[i]);
        }
        return sb.toString() + "\n";
    }

    private Set<String> primaryIds(String className) {
        Set<String> ids = new HashSet<String>();
        for (Item i : itemWriter.getItems()) {
            if (className.equals(i.getClassName()) && i.getAttribute("primaryIdentifier") != null) {
                ids.add(i.getAttribute("primaryIdentifier").getValue());
            }
        }
        return ids;
    }

    public void testGeneAndStrainGetWbPrefix() throws Exception {
        converter.process(new StringReader(row("WBRNAi1", "PMID:1", "WBStrain00007286", "WBGene00001831")));
        assertTrue(primaryIds("Gene").contains("WB:WBGene00001831"));
        assertTrue(primaryIds("Strain").contains("WB:WBStrain00007286"));
        assertEquals(1, primaryIds("Gene").size());
    }

    public void testPrefixedIdsAreKept() throws Exception {
        converter.process(new StringReader(row("WBRNAi1", "PMID:1", "WB:WBStrain00007286", "WB:WBGene00001831")));
        assertTrue(primaryIds("Gene").contains("WB:WBGene00001831"));
        assertTrue(primaryIds("Strain").contains("WB:WBStrain00007286"));
    }

    public void testIdWithoutColonFails() throws Exception {
        try {
            converter.process(new StringReader(row("WBRNAi1", "PMID:1", "", "N2")));
            fail("expected IllegalArgumentException");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("N2"));
        }
    }

    public void testBareWbPaperGetsPrefix() throws Exception {
        converter.process(new StringReader(row("WBRNAi1", "WBPaper00001", "", "WBGene00001831")));
        boolean found = false;
        for (Item i : itemWriter.getItems()) {
            if ("Publication".equals(i.getClassName())) {
                found = "WB:WBPaper00001".equals(i.getAttribute("pubXrefId").getValue());
            }
        }
        assertTrue(found);
    }
}
