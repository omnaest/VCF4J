/*******************************************************************************
 * Copyright 2021 Danny Kunz
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License.  You may obtain a copy
 * of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
 * License for the specific language governing permissions and limitations under
 * the License.
 ******************************************************************************/
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
