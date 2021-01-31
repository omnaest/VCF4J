package org.omnaest.genomics.vcf.components.parser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import org.junit.Test;
import org.omnaest.genomics.vcf.components.parser.VCFParser_4_1;
import org.omnaest.genomics.vcf.components.parser.VCFParser_4_1.CommentFilter;

public class VCFParser_4_1Test
{

    @Test
    public void testTest() throws Exception
    {
        CommentFilter commentFilter = new VCFParser_4_1.CommentFilter();

        {
            assertFalse(commentFilter.test("##fileformat=VCFv4.2"));
            assertEquals("VCFv4.2", commentFilter.getCommentMap()
                                                 .get("fileformat")
                                                 .iterator()
                                                 .next());
        }
        {
            assertFalse(commentFilter.test("##reference=/tmp/local_ngs_data//GRCh37.fa"));
            assertEquals("/tmp/local_ngs_data//GRCh37.fa", commentFilter.getCommentMap()
                                                                        .get("reference")
                                                                        .iterator()
                                                                        .next());
        }
    }

}
