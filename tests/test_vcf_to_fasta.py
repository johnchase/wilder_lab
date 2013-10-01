

# File created on 08 Feb 2012
from __future__ import division

__author__ = ""
__copyright__ = ""
__credits__ = ""
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = ""
__status__ = ""

from unittest import TestCase, main
from vcf_to_fasta import (indiv_snp_variation, process_header_entry, 
process_data_entry_line, parse_data_file, seq_slicer, get_start_stop_pos, 
create_entry_header, create_snp_dict, indiv_id_to_snp_dict, create_base_sub, 
create_seq, create_2_base_sub, create_2_seqs, zip_header_and_entry, 
create_entry_header_haploid, zip_header_and_entry_haploid, check_vcf, check_chrom_file)



fp = open('chr22.fa')
fp2 = open('chr5.fa')
#fp3 = open('chr10.fa')



class ExampleTests(TestCase):
    """ Tests the generation of lists from vcf files. """
    
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.example_file1 = example_file1 
        self.example_file2 = example_file2
        self.example_file3 = example_file3
        self.example_file4 = example_file4
        self.fp3 = open('chr10.fa')
        self.fp2 = open('chr5.fa')
    
##Test that indiv_snp_variation returns a list of the allele variation information. 
    def test_indiv_snp_variation(self): 
        """Does the function return correct output when given correct output"""
        input  = '10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.87,-5.00	0|0:0.000:-0.18,-0.48,-2.19	0|0:0.000:-0.12,-0.62,-3.85	0|0:0.000:-0.14,-0.56,-3.22	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.74,-5.00	0|0:0.000:-0.03,-1.19,-5.00	0|0:0.000:-0.05,-0.99,-5.00	0|0:0.000:-0.08,-0.75,-5.00	0|0:0.000:-0.02,-1.28,-5.00	0|0:0.000:-0.00,-2.07,-5.00	0|0:0.000:-0.01,-1.54,-5.00	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.03,-1.20,-5.00	0|0:0.000:-0.06,-0.91,-5.00	0|0:0.000:-0.00,-1.97,-5.00	0|0:0.000:-0.01,-1.51,-5.00'
        expected = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0],[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        self.assertEqual(indiv_snp_variation(input), expected)
       
    def test_indiv_snp_variation_2(self): 
        """Does the function return correct output when given correct output"""
        input  = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000'
        expected = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        self.assertEqual(indiv_snp_variation(input), expected)
        
        
#Tests for process_header_entry##
    def test_process_header_entry(self):
        """ Does this return a list of the fields when given the correct input?"""
        input = '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101'
        expected = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101']
        self.assertEqual(process_header_entry(input), expected)
        

##Tests for process_data_entry_line
    
    def test_process_data_entry_line(self):
        """Does this work when given correct input?"""
        ids = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101']
        line = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000'
        expected = ('10', '89623323', 'G', 'A', [('HG00096',[0,0]), ('HG00097', [0,0]), ('HG00099', [0,0]), ('HG00100', [0,0]), ('HG00101', [0,0])])
        self.assertEqual(process_data_entry_line(line, ids), expected)


    def test_process_data_entry_line_2(self):
        """Does this work when given correct input?"""
        ids = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110']
        line = '10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00'
        expected = ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])
        self.assertEqual(process_data_entry_line(line, ids), expected)


## Test for seq_slicer?
 
                
    def test_seq_slicer(self):     
        actual = fp
        input = (16050151, 16050251)
        expected = 'TTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAGGGTCGCAGGCCGGTTTCTGCTAATTTCTTTAATTCCAAGACAGTCTCAA'
        self.assertEqual(seq_slicer(actual, input ), expected) 

    def test_seq_slicer2(self):     
        actual = fp2
        input = (90466134, 90466314)
        expected = 'AAAGATGATATTTACTTCTCTCATCTCATTCTTGTTTCAAGACCTGTTAGTCATAAGGGCTTTGTGTGGAAAATTACATATGCATGAGTATTTCATAACCAAATGGTTGTAATCTCTTCTACTCATCTGGTTCCACAGTAGGAAGTCAAAAAACTGTGTATTCGTCATTTCATGGTTGATT'

        self.assertEqual(seq_slicer(actual, input), expected)
        
##Tests for parse_data_file##
    
    def test_parse_data_file(self):
        actual = list(parse_data_file(self.example_file1))
        expected =  [('10', '89623323', 'G', 'A', [('HG00096',[0,0]), ('HG00097', [0,0]), ('HG00099', [0,0]), ('HG00100', [0,0]), ('HG00101', [0,0])])]
        self.assertEqual(actual, expected)
        
    def test_parse_data_file_2(self): 
        actual = list(parse_data_file(self.example_file2))
        expected = [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
        self.assertEqual(actual, expected)
        
    def test_parse_data_file_3(self):
        actual = list(parse_data_file(self.example_file3))
        expected =  [('10', '89623323', 'GT', 'A,TC', [('HG00096',[0,0]), ('HG00097', [0,0]), ('HG00099', [0,0]), ('HG00100', [0,0]), ('HG00101', [0,0])])]
        self.assertEqual(actual, expected)
        
    def test_get_start_stop_pos1(self): 
        '''Test to determine that given input as a list the function correctly returns the stop and start postion.''' 
        input = [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]),
                                               ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]),
                                               ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), 
                 ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]),
                                               ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]),
                                               ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
        expected = [89674917, 89674997]
        self.assertEqual(get_start_stop_pos(input), expected)
        
    def test_create_entry_header(self): 
        '''does the function return the correct header given a the list from parse_data_file'''
        input =  [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
        expected = ['>HG00096 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00097 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00099 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00100 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00101 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00102 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00103 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00104 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00106 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00108 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00109 hg19_dna range=chr10:89674917-89674997', 
                    '>HG00110 hg19_dna range=chr10:89674917-89674997']    
        self.assertEqual(create_entry_header(input), expected)
        
    def test_create_entry_header_haploid(self):
        input = [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0])])]
        expected = [['>HG00096_1 hg19_dna range=chr10:89674917-89674997',
                    '>HG00096_2 hg19_dna range=chr10:89674917-89674997'],  
                   ['>HG00097_1 hg19_dna range=chr10:89674917-89674997',
                    '>HG00097_2 hg19_dna range=chr10:89674917-89674997']]
        self.assertEqual(create_entry_header_haploid(input), expected)
        
    def test_create_snp_dict(self): 
        '''does the function return a valid dictionary with correct input?'''
        input = [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
        expected = {89674917: ('T', 'G'), 89674997: ('A', 'G')}
        self.assertEqual(create_snp_dict(input), expected)
    
    def test_indiv_id_to_snp_dict(self): 
        input = [('10', '89674917', 'T', 'G', [('HG00096', [0, 1]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
        expected = {'HG00110': [(89674917, [0, 0]), (89674997, [0, 0])],
                    'HG00097': [(89674917, [0, 0]), (89674997, [0, 0])],
                    'HG00096': [(89674917, [0, 1]), (89674997, [0, 0])],
                    'HG00099': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00103': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00102': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00101': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00100': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00106': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00104': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00109': [(89674917, [0, 0]), (89674997, [0, 0])], 
                    'HG00108': [(89674917, [0, 0]), (89674997, [0, 0])]}
        self.assertEqual(indiv_id_to_snp_dict(input), expected) 
        
    def test_create_base_sub(self):
        alleles = [0, 1]
        ref = 'A'
        alt = 'T'
        expected = 'W'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
    
    def test_create_base_sub_1(self): 
        alleles = [1, 2] 
        ref = 'C'
        alt = 'A,C,G,T' 
        expected = 'M'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
        
    def test_create_base_sub_2(self):
        alleles = [0, 1] 
        ref  = 'GT'
        alt = 'G, GTC'
        expected = 'GT'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
        
    def test_create_base_sub_3(self):
        alleles = [0, 1]
        ref = 'T'
        alt = '.'
        expected = 'T'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
    
    def test_create_base_sub_4(self):
        alleles = [0, 1]
        ref = 'TC'
        alt = 'T'
        expected = 'TC'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
    
    def test_create_base_sub_5(self):
        alleles = [1, 2]
        ref = 'TCG'
        alt = 'TG,T,TCAG'
        expected = 'TCG'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)

    def test_create_base_sub_6(self):
        alleles = [1, 2]
        ref = 'T'
        alt = 'TAA'
        expected = 'T'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)
    
    def test_create_base_sub_7(self): 
        alleles = [1, 1]
        ref = 'T'
        alt = 'C'
        expected = 'C'
        self.assertEqual(create_base_sub(alleles, ref, alt), expected)

#Test the base creation for diploid sequences. 
    def test_create_2_base_sub1(self): 
        alleles = [0, 1]
        ref = 'C' 
        alt = 'G' 
        expected = ['C', 'G']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected) 
    
    def test_create_2_base_sub2(self): 
        alleles = [0, 1]
        ref = 'TC' 
        alt = 'T' 
        expected = ['TC', 'T-']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)

    def test_create_2_base_sub3(self): 
        alleles = [0, 1]
        ref = 'TC' 
        alt = 'TCA' 
        expected = ['TC-', 'TCA']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)
        
    def test_create_2_base_sub4(self): 
        alleles = [0, 1]
        ref = 'TC' 
        alt = 'TG, T' 
        expected = ['TC', 'TG']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)

    def test_create_2_base_sub5(self): 
        alleles = [0, 1]
        ref = 'TC' 
        alt = 'TG, T' 
        expected = ['TC', 'TG']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)
        
    def test_create_2_base_sub6(self): 
        alleles = [0, 1]
        ref = 'TCG'    
        alt = 'TG, T, TCAG' 
        expected = ['TC-G', 'T--G']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)
        
    def test_create_2_base_sub7(self): 
        alleles = [0, 1]
        ref = 'TCG'    
        alt = '., T, TCAG' 
        expected = ['TC-G', '----']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)
        
    def test_create_2_base_sub8(self): 
        alleles = [0, 1]
        ref = 'CTTTTGACTCCCCAAAACTTAATATTTAGCCTATACTTGACTAGAAGTCTTACTGATGACATAATGTTCGTTAATACATATTTTATATATGTGTCAGATAGCATATTTGTATAATAAAGTAAGCTGCAGGAAAAATATTAAAATCATAAAGAAGAGAAAATATACTTACTATTCATTAAGTGGAAGTGGATCCTCATAAAGGTCTTCATCCTCACTGCCTTCACTTTGAGTAGGCCGAGGAGTAGGAGAGAGAGGAAAGGTCAGACTTGCTGTCTCATGGGTGGCAGAGGTAGAAGAAGGTCCACATACAAGTGGTCCGACACAGCTCAAACCGGTTTTGTTCATTGGCCAACTGTAGTTTGATTGAAAGTAATAATAAATGAAGTTTCTGCCTCAGTTCAGTATTATCAAGTCATAGATAGCAAGGGCTGGAAGAAACCTTAGTAGTAATCTCTTTGAGTCTAATTATCATGTAGAATAGGAAATTGCGGTCTAGAAAGGTTAAGTGACTTGTCCAAATTACACAACTAGTTAGAGACATAGCCAGCTCTTAAATCTGACTTCCAGATTTTCACTGTGTCTTCTTTTTTCTGTAACGTGTTGCCTTTTTTAGCCATGAAAAATTAGAAGTTGAACTCTTGTCTTTTCAGGCAGGTGTCAATTTTGGGGTTTTGTTTTGATTTTTGGTTTTTGACATAAAGTACTTTAGTTCTGTGATGTATAAACCGTGAGTTTCTGTTTTTCTCATATACCTGAATACTGTCCATGTGGAAGTTACCTTTTATCTTTACCAGTATTAACACATAAATGGTTATACATAAATACATTGACCACCTTTTATTACTCCAGCTATAGTGGGGAAAACTTTCTTTTCATAACTAGCTAATGTTTTAAAAAGTA'    
        alt = 'C' 
        expected = ['CTTTTGACTCCCCAAAACTTAATATTTAGCCTATACTTGACTAGAAGTCTTACTGATGACATAATGTTCGTTAATACATATTTTATATATGTGTCAGATAGCATATTTGTATAATAAAGTAAGCTGCAGGAAAAATATTAAAATCATAAAGAAGAGAAAATATACTTACTATTCATTAAGTGGAAGTGGATCCTCATAAAGGTCTTCATCCTCACTGCCTTCACTTTGAGTAGGCCGAGGAGTAGGAGAGAGAGGAAAGGTCAGACTTGCTGTCTCATGGGTGGCAGAGGTAGAAGAAGGTCCACATACAAGTGGTCCGACACAGCTCAAACCGGTTTTGTTCATTGGCCAACTGTAGTTTGATTGAAAGTAATAATAAATGAAGTTTCTGCCTCAGTTCAGTATTATCAAGTCATAGATAGCAAGGGCTGGAAGAAACCTTAGTAGTAATCTCTTTGAGTCTAATTATCATGTAGAATAGGAAATTGCGGTCTAGAAAGGTTAAGTGACTTGTCCAAATTACACAACTAGTTAGAGACATAGCCAGCTCTTAAATCTGACTTCCAGATTTTCACTGTGTCTTCTTTTTTCTGTAACGTGTTGCCTTTTTTAGCCATGAAAAATTAGAAGTTGAACTCTTGTCTTTTCAGGCAGGTGTCAATTTTGGGGTTTTGTTTTGATTTTTGGTTTTTGACATAAAGTACTTTAGTTCTGTGATGTATAAACCGTGAGTTTCTGTTTTTCTCATATACCTGAATACTGTCCATGTGGAAGTTACCTTTTATCTTTACCAGTATTAACACATAAATGGTTATACATAAATACATTGACCACCTTTTATTACTCCAGCTATAGTGGGGAAAACTTTCTTTTCATAACTAGCTAATGTTTTAAAAAGTA', 'C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------']
        self.assertEqual(create_2_base_sub(alleles, ref, alt), expected)
                   
#    def test_get_snp_list(self):
#        input = [('10', '89674917', 'T', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])]
#        expected = [89674917, 89674997]
#        self.assertEqual(get_snp_list(input), expected)
    
    def test_create_seq(self):
        file = self.fp3
        actual = [('10', '89674917', 'T', 'G', [('HG00096', [1, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 1])])]
        expected = ['KAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAR']
        self.assertEqual(create_seq(file, actual), expected)
        
    def test_create_seq1(self): 
        file = self.fp3
        actual = [('10', '89674917', 'T', 'G', [('HG00096', [1, 0]), ('HG00097', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 1]), ('HG00097', [0, 0])])]
        expected = ['KAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAR', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA']
        self.assertEqual(create_seq(file, actual), expected)
        
    def test_create_2_seqs(self): 
        file = self.fp3
        actual = [('10', '89674917', 'T', 'G', [('HG00096', [1, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 1])])]
        expected = [['GAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG']]
        self.assertEqual(create_2_seqs(file, actual), expected)
            
            
    def test_create_2_seqs1(self): 
        file = self.fp3
        actual = [('10', '89674917', 'T', 'G', [('HG00096', [1, 0]), ('HG00097', [0, 0])]), ('10', '89674997', 'A', 'G', [('HG00096', [0, 1]), ('HG00097', [1, 0])])]
        expected = [['GAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG'], ['TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA']]
        self.assertEqual(create_2_seqs(file, actual), expected)

    def test_zip_header_and_entry(self):
        seq_list = ['KAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAR', 
                        'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA'] 
        header_list = ['>HG00096 hg19_dna range=chr10:89674917-89674997', 
                       '>HG00097 hg19_dna range=chr10:89674917-89674997']
        expected = [('>HG00096 hg19_dna range=chr10:89674917-89674997',
                      'KAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAR'),          
                    ('>HG00097 hg19_dna range=chr10:89674917-89674997',
                      'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA')]
        self.assertEqual(zip_header_and_entry(seq_list, header_list), expected) 


    def test_zip_header_and_entry_haploid(self):
        seq_list = [['GAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA', 
                         'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG'], 
                        ['TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG', 
                         'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA']]
        header_list = [['>HG00096_1 hg19_dna range=chr10:89674917-89674997',
                    '>HG00096_2 hg19_dna range=chr10:89674917-89674997'],  
                   ['>HG00097_1 hg19_dna range=chr10:89674917-89674997',
                    '>HG00097_2 hg19_dna range=chr10:89674917-89674997']]
        expected = [[('>HG00096_1 hg19_dna range=chr10:89674917-89674997', 'GAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA'), 
                     ('>HG00096_2 hg19_dna range=chr10:89674917-89674997', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG')], 
                    [('>HG00097_1 hg19_dna range=chr10:89674917-89674997', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAG'), 
                     ('>HG00097_2 hg19_dna range=chr10:89674917-89674997', 'TAGGATTACAGGTGTGAGCCACTGCACCCAGCTCTTAGACAAATTTTTTATTCCAAACTTTTTTTATTTTATCATTTGAAA')]]
        self.assertEqual(zip_header_and_entry_haploid(seq_list, header_list), expected)
        
        
    def test_check_chrom_file(self): 
        self.assertRaises(ValueError, check_chrom_file, self.example_file3, self.fp2)
        
#    def test_check_vcf1(self): 
#        self.assertRaises(ValueError, check_vcf, self.example_file4)

#38tests#
        
        
example_file1 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000""".split('\n')
        
        
example_file2 = """##reference=GRCh37																				
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102	HG00103	HG00104	HG00106	HG00108	HG00109	HG00110
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.87,-5.00	0|0:0.000:-0.18,-0.48,-2.19	0|0:0.000:-0.12,-0.62,-3.85	0|0:0.000:-0.14,-0.56,-3.22	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.74,-5.00	0|0:0.000:-0.03,-1.19,-5.00	0|0:0.000:-0.05,-0.99,-5.00	0|0:0.000:-0.08,-0.75,-5.00	0|0:0.000:-0.02,-1.28,-5.00
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00""".split('\n')

example_file3 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	GT	A,TC	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000""".split('\n')

example_file4 = '''##fileformat=VCFv4.1										
"##INFO=<ID=LDAF,Number=1,Type=Float,Description=""MLE Allele Frequency Accounting for LD"">"										
"##INFO=<ID=AVGPOST,Number=1,Type=Float,Description=""Average posterior probability from MaCH/Thunder"">"										
"##INFO=<ID=RSQ,Number=1,Type=Float,Description=""Genotype imputation quality from MaCH/Thunder"">"										
"##INFO=<ID=ERATE,Number=1,Type=Float,Description=""Per-marker Mutation rate from MaCH/Thunder"">"										
"##INFO=<ID=THETA,Number=1,Type=Float,Description=""Per-marker Transition rate from MaCH/Thunder"">"										
"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=""Confidence interval around END for imprecise variants"">"										
"##INFO=<ID=SNPSOURCE,Number=.,Type=String,Description=""indicates if a snp was called when analysing the low coverage or exome alignment data"">"										
##reference=GRCh37										
##reference=GRCh37										
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	"0|0:0.000:-0.03,-1.22,-5.00"	"0|0:0.000:-0.01,-1.62,-5.00"
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	"0|0:0.000:-0.03,-1.23,-5.00"	"0|0:0.000:-0.02,-1.33,-5.00"
10	89675029	rs151009112	T	C	100	PASS	ERATE=0.0005;AA=T;AN=2184;LDAF=0.0018;THETA=0.0005;VT=SNP;RSQ=0.4609;SNPSOURCE=LOWCOV;AC=2;AVGPOST=0.9974;AF=0.0009;ASN_AF=0.0035	GT:DS:GL	"0|0:0.000:-0.03,-1.23,-5.00"	"0|0:0.000:-0.04,-1.05,-5.00"
10	89675036	rs111627758	C	T	100	PASS	ERATE=0.0004;LDAF=0.0047;AA=T;AN=2184;RSQ=0.8612;THETA=0.0005;VT=SNP;AC=9;SNPSOURCE=LOWCOV;AVGPOST=0.9984;AF=0.0041;AMR_AF=0.01;EUR_AF=0.01	GT:DS:GL	"0|0:0.000:-0.03,-1.16,-5.00"	"0|0:0.000:-0.02,-1.28,-5.00"
10	89675296	rs1234224	A	G	100	PASS	ERATE=0.0004;AC=915;RSQ=0.9931;THETA=0.0004;AA=G;AN=2184;VT=SNP;LDAF=0.4188;SNPSOURCE=LOWCOV;AVGPOST=0.9961;AF=0.42;ASN_AF=0.48;AMR_AF=0.43;AFR_AF=0.42;EUR_AF=0.37	GT:DS:GL	"0|0:0.000:-0.01,-1.50,-5.00"	"1|0:1.000:-5.00,-0.00,-4.40"'''.split('\n')




if __name__ == "__main__":
    main()
    