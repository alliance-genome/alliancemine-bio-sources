package org.intermine.bio.dataconversion;

import junit.framework.TestCase;

public class AllianceGenesCoordinatesTest extends TestCase
{
    public void testNullBecomesEmpty() {
        assertEquals("", AllianceGenesConverter.nullToEmpty("null"));
        assertEquals("chrI", AllianceGenesConverter.nullToEmpty("chrI"));
    }

    public void testStrandIsNormalised() {
        assertEquals("1", AllianceGenesConverter.normaliseStrand("+"));
        assertEquals("-1", AllianceGenesConverter.normaliseStrand("-"));
        assertEquals("1", AllianceGenesConverter.normaliseStrand("1"));
        assertEquals("", AllianceGenesConverter.normaliseStrand("null"));
    }

    public void testLengthIsInclusive() {
        assertEquals("10", AllianceGenesConverter.getLength("1", "10"));
        assertEquals("10", AllianceGenesConverter.getLength("10", "1"));
        assertEquals("1", AllianceGenesConverter.getLength("5", "5"));
    }
}
