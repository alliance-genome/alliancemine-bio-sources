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

public class GoAnnotationConverterTest extends TestCase
{
    private static String gaf(String pub) {
        return "SGD\tS000000001\tACT1\t\tGO:0005515\t" + pub + "\tIPI\tSGD:S000000002\tF\t\t\t"
                + "gene\ttaxon:559292\t20200101\tSGD\t\t\n";
    }

    private int countAnnotations(String gafText) throws Exception {
        MockItemWriter writer =
                new MockItemWriter(new HashMap<String, org.intermine.model.fulldata.Item>());
        new GoConverter(writer, Model.getInstanceByName("genomic"))
                .process(new StringReader(gafText));
        int annotations = 0;
        for (Item item : writer.getItems()) {
            if ("GOAnnotation".equals(item.getClassName())) {
                annotations++;
            }
        }
        return annotations;
    }

    public void testRepeatedAnnotationWithWithTextIsMerged() throws Exception {
        // same gene, term, qualifier, with-text and publication on two lines
        assertEquals(1, countAnnotations(gaf("PMID:1") + gaf("PMID:1")));
    }
}
