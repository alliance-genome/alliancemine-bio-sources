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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;

import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

public class AllianceMolecularInteractionsConverterTest extends TestCase
{
    private AllianceMolecularInteractionsConverter converter;
    private MockItemWriter itemWriter;

    public void setUp() throws Exception {
        Model model = Model.getInstanceByName("genomic");
        itemWriter = new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        converter = new AllianceMolecularInteractionsConverter(itemWriter, model);
    }

    private static String row(String gene, String partner, String roleA, String roleB, String uid) {
        String[] f = new String[29];
        java.util.Arrays.fill(f, "");
        f[0] = gene;
        f[2] = "NCBITaxon:9606";
        f[3] = partner;
        f[5] = "NCBITaxon:9606";
        f[8] = "MI:0326";
        f[10] = "MI:0328";
        f[12] = roleA;
        f[14] = roleB;
        f[18] = "interacts_with";
        f[20] = uid;
        f[21] = "PMID:1";
        f[23] = "Smith 2020";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < f.length; i++) {
            sb.append(i == 0 ? "" : "\t").append(f[i]);
        }
        return sb.toString();
    }

    private List<Item> itemsOf(String className) {
        List<Item> out = new ArrayList<Item>();
        for (Item i : itemWriter.getItems()) {
            if (className.equals(i.getClassName())) {
                out.add(i);
            }
        }
        return out;
    }

    private Map<String, String> geneIds() {
        Map<String, String> ids = new HashMap<String, String>();
        for (Item g : itemsOf("Gene")) {
            ids.put(g.getAttribute("primaryIdentifier").getValue(), g.getIdentifier());
        }
        return ids;
    }

    public void testBothParticipantsGetTheInteraction() throws Exception {
        converter.process(new StringReader(row("HGNC:1", "HGNC:2", "MI:0496", "MI:0498", "u1") + "\n"));
        Map<String, String> genes = geneIds();
        List<Item> interactions = itemsOf("Interaction");
        assertEquals(2, interactions.size());
        int forward = 0;
        int reverse = 0;
        for (Item i : interactions) {
            String p1 = i.getReference("participant1").getRefId();
            String p2 = i.getReference("participant2").getRefId();
            if (p1.equals(genes.get("HGNC:1")) && p2.equals(genes.get("HGNC:2"))) {
                forward++;
            } else if (p1.equals(genes.get("HGNC:2")) && p2.equals(genes.get("HGNC:1"))) {
                reverse++;
            }
        }
        assertEquals(1, forward);
        assertEquals(1, reverse);
    }

    public void testMirrorDetailSwapsRolesAndTypes() throws Exception {
        converter.process(new StringReader(row("HGNC:1", "HGNC:2", "MI:0496", "MI:0498", "u1") + "\n"));
        Map<String, String> genes = geneIds();
        Map<String, String> interactionP1 = new HashMap<String, String>();
        for (Item i : itemsOf("Interaction")) {
            interactionP1.put(i.getIdentifier(), i.getReference("participant1").getRefId());
        }
        List<Item> details = itemsOf("InteractionDetail");
        assertEquals(2, details.size());
        for (Item d : details) {
            boolean forward = interactionP1.get(d.getReference("interaction").getRefId())
                    .equals(genes.get("HGNC:1"));
            assertEquals(forward ? "MI:0496" : "MI:0498", d.getAttribute("role1").getValue());
            assertEquals(forward ? "MI:0498" : "MI:0496", d.getAttribute("role2").getValue());
        }
    }

    public void testExperimentIsSharedByBothDetails() throws Exception {
        converter.process(new StringReader(row("HGNC:1", "HGNC:2", "MI:0496", "MI:0498", "u1") + "\n"));
        assertEquals(1, itemsOf("InteractionExperiment").size());
        for (Item d : itemsOf("InteractionDetail")) {
            assertNotNull(d.getReference("experiment"));
        }
    }

    public void testSelfInteractionIsNotMirrored() throws Exception {
        converter.process(new StringReader(row("HGNC:1", "HGNC:1", "MI:0496", "MI:0498", "u1") + "\n"));
        assertEquals(1, itemsOf("Interaction").size());
        assertEquals(1, itemsOf("InteractionDetail").size());
    }

    public void testDuplicateUniqueIdIsStoredOnce() throws Exception {
        converter.process(new StringReader(row("HGNC:1", "HGNC:2", "MI:0496", "MI:0498", "u1") + "\n"
                + row("HGNC:2", "HGNC:1", "MI:0498", "MI:0496", "u1") + "\n"));
        assertEquals(2, itemsOf("Interaction").size());
    }
}
