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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class AllianceAlleleDetailConverterTest
{
    @Test
    public void taxonFollowsCuriePrefix() {
        assertEquals("559292", AllianceAlleleDetailConverter.taxonForAlleleId("SGD:S000277594"));
        assertEquals("10090", AllianceAlleleDetailConverter.taxonForAlleleId("MGI:2173"));
        assertEquals("7955", AllianceAlleleDetailConverter.taxonForAlleleId("ZFIN:ZDB-ALT-131217-13114"));
    }

    @Test
    public void unknownOrAmbiguousPrefixHasNoTaxon() {
        assertNull(AllianceAlleleDetailConverter.taxonForAlleleId("Xenbase:XB-ALT-1"));
        assertNull(AllianceAlleleDetailConverter.taxonForAlleleId("noprefix"));
    }
}
