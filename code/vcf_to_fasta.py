from cogent.parse.fasta import MinimalFastaParser
from cogent.app.muscle import align_unaligned_seqs
from cogent import LoadSeqs, DNA
from sys import argv


#create a list of all of the information related to each SNP and the individual results. This data comes from  
def process_data_entry_line (line, ids):
    fields = line.split('\t') 
    chr_n = fields[0] 
    pos_n = fields[1]
    ref = fields[3] 
    alt = fields[4]
    indiv_ids = zip(ids, indiv_snp_variation(line))
    return (chr_n, pos_n, ref, alt, indiv_ids)

#process the information line. This line contains all of the information about the individuals. 
#The function returns a list of the individuals 
def process_header_entry(line):
    fields = line.split() 
    return fields [9:]

#Creates a list of allelic variation for each SNP. The list itself contains one entry of two numbers for each individual. 
def indiv_snp_variation(line): 
    l = line.split('\t')
    list = l[9:]
    variation_list = [] 
    for i in list: 
        variation_list.append(map(int,i.split(':')[0].split('|')))
    return variation_list


def parse_data_file(f): ##name this something else
    ids = None
    for line in f:
        if line.startswith('##'):
            pass
        elif line.startswith('#CHROM'): 
            ids = process_header_entry(line)
        else:
            if ids == None:
                raise ValueError, "Didn't find '#CHROM' line before data lines. Can't continue."
            result = process_data_entry_line(line, ids)
            yield result
            #return result


#This function takes specific sequences from whole chromosome reads. Currently the fasta file must be in the >chrN format. It wil not 
#work with fasta files from the UCSC human genome browser.
def seq_slicer(f, limits):
    fa = [(name, seq) for name, seq in MinimalFastaParser(f)]
    return fa[0][1][limits[0] - 1:limits[1]]

# This function gets the start and stop position of the sequence. This is based on the first and last entry in the vcf file.
#The input here is the output of parse_data_file
def get_start_stop_pos(line):
    limits = []
    limits.append(int(line[0][1]))
    limits.append(int(line[-1][1]))
    return limits


# create the header for each entry. These will later be zipped together with the sequence for each entry. 
def create_entry_header(line):
    ids = get_ids(line) 
    limits = get_start_stop_pos(line)
    header_list = []
    for id in ids:
        header =  ">%s hg19_dna range=chr%d:%d-%d" %(id, int(line[0][0]), limits[0], limits[1])
        header_list.append(header)
    return header_list
    
# creates a list of headers with two entries for each individual, one for each sequence.     
def create_entry_header_haploid(line, genome='hg19_dna'):
    ids = get_ids(line) 
    limits = get_start_stop_pos(line)
    header_list = []
    for id in ids:
        headers = []
        for i in range(1,3):
            header = ">%s_%s %s range=chr%d:%d-%d" %(id, i, genome, int(line[0][0]), limits[0], limits[1])
            headers.append(header) 
        header_list.append(headers)
    return header_list
        
        
        
#create a dictionary relating the snp location to the possible snps.
def create_snp_dict(line):
    snp_dic = {}
    for i in line: 
        snp_dic[int(i[1])] = (i[2], i[3])
    return snp_dic

# Get a list of ids from the parse_data_file output. 
def get_ids(line): 
    ids = []
    for i in line[0][4]:
        ids.append(i[0])
    return ids
      
#create a dictionary where the key is the individual id and the value is the snp location and the allele combination, 
#reference and alternate alleles are not included in this dictionary but can be found in creat_snp_dict
def indiv_id_to_snp_dict(line):
    ids = get_ids(line)   
    id_master = {}
    for i in range(len(ids)):
        mini_id = []
        for j in range(len(line)): 
            mini_id.append((int(line[j][1]), line[j][4][i][1]))
        id_master[(ids[i])] = mini_id
    return id_master

#return the ambiguous base for each snp in each sequence, if there is one. 
def create_base_sub(alleles, ref, alt):
    bases = alt.split(',')
    bases.insert(0, ref)
    ref1 = 'none'
#for the haploid data, indels are going to be ignored. 
    for base in bases: 
        if len(base) > 1 or base == '.': 
            ref1 = ref
            return ref1
        else: 
            pass
    if ref1 == 'none': 
        allele1 = bases[alleles[0]]
        allele2 = bases[alleles[1]]
        #this prevents the function from trying to look up 'AA' (for instance) in the ambiguous_base_dic.
        if allele1 == allele2: 
            return allele1
        else: 
            amb_base_key = allele1 + allele2
            return ambiguous_base_dic[amb_base_key] 

#the ambiguous bases for the haploid data set. 
ambiguous_base_dic = {'AG': 'R', 'CT': 'Y', 'GT': 'K', 
                      'AC': 'M', 'GC': 'S', 'AT': 'W',
                      'GA': 'R', 'TC': 'Y', 'TG': 'K',
                      'CA': 'M', 'CG': 'S', 'TA': 'W'}

#This function is probably not going to be necessary, though may be useful for diploid data.#
#def get_snp_list(line): 
#    snps = []
#    for i in line: 
#        snps.append(int(i[1]))
#    return snps
    
#There has got to be a better way to do this...
#create a list of all of the diploid sequences. This list will be zipped with the header list. 
def create_seq(chrom, line):
    limits = get_start_stop_pos(line)
    id_list = get_ids(line)
    id_dict = indiv_id_to_snp_dict(line)
    snp_dict = create_snp_dict(line)
    seq = seq_slicer(chrom, limits)
    master_list = []
    for id in id_list:
        seq1 = seq.upper()
        # the id list is a list of each individual in the file. This is created to keep the entries in order.
        for snp in id_dict[id]:
        #the snp here is in the format [(89674917, [0, 0]), (89674997, [0, 0])]
            s_pos = snp[0] - limits[0]
            #returns the position of the snp in the sliced sequence. 
            alleles = snp[1]
            ref = snp_dict[snp[0]][0]
            alt = snp_dict[snp[0]][1]
            if alleles == [0, 0]: 
                pass
            else:
                new_base = create_base_sub(alleles, ref, alt)
                seq1 = seq1[:s_pos] + new_base + seq1[s_pos + len(ref):]
        master_list.append(seq1)
    return master_list


#Create the base substitutions that will be applied to the the haploid sequences. This will always return the longest possible
#sequence based on the possible alleles. 
def create_2_base_sub(alleles, ref, alt): 
    bases_list = alt.split(', ')
    bases_list.insert(0, ref)
    seqs = [bases_list[alleles[0]], bases_list[alleles[1]]]  
    scaffold = max(bases_list, key=len)
    substitution = []
    for seq in seqs:
        if len(seq) == len(scaffold): 
            substitution.append(seq)
        elif seq == '.':
            substitution.append('-'*len(scaffold))
        else:
            seq_and_scaffold = [seq, scaffold]
            aligned_seqs = align_unaligned_seqs(seq_and_scaffold, DNA, params={'-gapopen': -0.001})
            aligned_seqs = str(aligned_seqs)
            bases_sub = aligned_seqs.split('\n')
            substitution.append(bases_sub[3])
    return substitution
    
#Create the haploid sequences, two for each id.     
def create_2_seqs(file, line): 
    limits = get_start_stop_pos(line)
    id_list = get_ids(line)
    id_dict = indiv_id_to_snp_dict(line)
    snp_dict = create_snp_dict(line)
    seq = seq_slicer(file, limits)
    master_list = []
    for id in id_list:
        seq1 = seq.upper()
        seq2 = seq.upper()
        for snp in id_dict[id]:
            s_pos = snp[0] - limits[0]
            alleles = snp[1]
            ref = snp_dict[snp[0]][0]
            alt = snp_dict[snp[0]][1]
            if alleles == [0, 0]: 
                pass
            else:
                new_base = create_2_base_sub(alleles, ref, alt)
                seq1 = seq1[:s_pos] + new_base[0] + seq1[s_pos + len(new_base[0]):]
                seq2 = seq2[:s_pos] + new_base[1] + seq2[s_pos + len(new_base[1]):]
        master_list.append([seq1, seq2])
    return master_list

# zip header and seq for diploid sequences. 
def zip_header_and_entry(seq_list, header_list):
        file_entries = zip(header_list, seq_list)
        return file_entries
        
# zip header and seq for diploid sequences.         
def zip_header_and_entry_haploid(seq_list, header_list):
    file_entries = []
    for i in range(len(header_list)): 
            seq_header = zip (header_list[i], seq_list[i])
            file_entries.append(seq_header)
    return file_entries
        
#check to make sure that the chromosome file specified on the command line is the correct chromosome based on the vcf file. 
def check_chrom_file(file, chrom_file):
    file_info = list(parse_data_file(file))
    chrom = file_info[0][0]
    if list(chrom_file)[0].startswith('>chr%s' %chrom):
        return file, chrom_file
    else: 
       raise ValueError, "Chromosome in fasta file does not match chromosome in vcf file or fasta file is in improper format. Fasta file must be in format '>chrN seq'"
    
###########
#Still working on these functions
###########
#NEED TESTING!!!!
###########



def is_number(s):
    try:
        float(s)
        return True
    except:
        return False



def check_vcf(file): 
    for line in file: 
        if line.startswith('##'):
            pass
        elif line.startswith('#CHROM'):
            pass 
        elif is_number(line[0]):
            pass
        else: 
            raise ValueError, "vcf file is in improper format, or may have been modified using excel or another program. Check vcf file in text editor, quotations are not allowed in any field."
        return file




##############
# End function and variable definitions
##############


##############
# Start main execution block
##############


def main():
    script_name, file_name, chrom, output_fasta_fp, ploidy = argv
    
    
    vcf = open(file_name, 'U')
    chrom = open(chrom)
    check_vcf(vcf) 
#    check_chrom_file(vcf,chrom)
    line = list(parse_data_file(vcf))
    

    if ploidy == 'D':
        header_list = create_entry_header(line)
        seq_list = create_seq(chrom, line)
        output = zip_header_and_entry(seq_list, header_list)
        output_file = open(output_fasta_fp,'w')
        for i in output:
            output_file.write(str(i[0]))
            output_file.write('\n')
            output_file.write(str(i[1]))    
            output_file.write('\n'*2)
        output_file.close()               
                       
    elif ploidy == 'H':
        header_list = create_entry_header_haploid(line)
        seq_list = create_2_seqs(chrom, line)
        output = zip_header_and_entry_haploid(seq_list, header_list)
        output_file = open(output_fasta_fp,'w')
        for i in output:
            for j in i:
                output_file.write(str(j[0]))
                output_file.write('\n')
                output_file.write(str(j[1]))
                output_file.write('\n'*2)
        output_file.close()
       

if __name__ == "__main__":
    main()

##############
# End main execution block
##############
