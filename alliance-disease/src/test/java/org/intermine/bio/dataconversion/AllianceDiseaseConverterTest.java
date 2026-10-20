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

import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.Map;

import junit.framework.TestCase;

import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

public class AllianceDiseaseConverterTest extends TestCase
{
    private AllianceDiseaseConverter converter;
    private MockItemWriter itemWriter;
    private Method publication;

    public void setUp() throws Exception {
        Model model = Model.getInstanceByName("genomic");
        itemWriter = new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        converter = new AllianceDiseaseConverter(itemWriter, model);
        publication = AllianceDiseaseConverter.class.getDeclaredMethod("newPublication", String.class);
        publication.setAccessible(true);
    }

    private String pub(String token) throws Exception {
        return (String) publication.invoke(converter, token);
    }

    private Item find(String ref) {
        for (Item i : itemWriter.getItems()) {
            if (i.getIdentifier().equals(ref)) {
                return i;
            }
        }
        return null;
    }

    public void testPmidPrefixIsStripped() throws Exception {
        Item item = find(pub(" PMID:12345 "));
        assertEquals("12345", item.getAttribute("pubMedId").getValue());
    }

    public void testBareDigitsGoToPubMedId() throws Exception {
        Item item = find(pub("12345"));
        assertEquals("12345", item.getAttribute("pubMedId").getValue());
    }

    public void testOtherIdGoesToPubXrefId() throws Exception {
        Item item = find(pub("MGI:6194238"));
        assertEquals("MGI:6194238", item.getAttribute("pubXrefId").getValue());
        assertNull(item.getAttribute("pubMedId"));
    }

    public void testBareWbPaperGetsPrefix() throws Exception {
        Item item = find(pub("WBPaper00001"));
        assertEquals("WB:WBPaper00001", item.getAttribute("pubXrefId").getValue());
    }

    public void testNonPmidIdIsCached() throws Exception {
        assertEquals(pub("MGI:6194238"), pub("MGI:6194238"));
        assertEquals(pub("PMID:7"), pub("7"));
    }

    public void testEmptyTokenGivesNoPublication() throws Exception {
        assertNull(pub("  "));
        assertNull(pub("PMID:"));
    }
}
