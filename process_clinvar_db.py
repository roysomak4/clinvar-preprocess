import sys
import gzip
import progressbar
import sh


def process_clinvar_db(args):
    # check if required arguments are present
    if len(args) == 0:
        print('Missing arguments...')
        print('Expected command: python process_clinvar_db.py clinvar_releasedate.vcf.gz')
        sys.exit(1)
    
    # get clinvar file
    clinvar_vcf = args[0]

    # Parse vcf file and extract relevant fields
    out_buffer = []
    genome_ref, release_date = get_vcf_fields(clinvar_vcf, out_buffer)

    # write buffer to output file
    write_output(out_buffer, genome_ref, release_date)


def write_output(out_buffer, genome_ref, release_date):
    print('Preparing to write output file...')
    # output clinvar file
    print genome_ref
    print release_date
    
    clinvar_out_file = 'clinvar_{0}_{1}_processed.txt.gz'.format(genome_ref, release_date)
   
    with gzip.open(clinvar_out_file, 'wb') as outf:
        header = 'chr\tpos\tref\talt\tclinvar_id\tclin_sig\treview_status\t'
        header += 'stars\tdisease\tallele_id\tdbsnp_id\n'
        # write header
        outf.write(header)
        print('File header written.')
        # write variant rows
        outf.write('\n'.join(out_buffer) + '\n')
        print('Output file written to {0}.'.format(clinvar_out_file))


def get_vcf_fields(filename, out_buffer):
    genome_ref = ''
    release_date = ''
    print('Calculating VCF file size...')
    gz_app = get_gzip_app()
    row_count = int(sh.bash("-c", "{app} {vcf_file} | wc -l"
                            .format(app=gz_app,
                                    vcf_file=filename)))

    is_vcf_format = False
    clinvar_star_rating = {
        'no_assertion_criteria_provided': '0',
        'no_assertion_provided': '0',
        'no_interpretation_for_the_single_variant': '0',
        'criteria_provided,_single_submitter': '1',
        'criteria_provided,_conflicting_interpretations': '1',
        'criteria_provided,_multiple_submitters,_no_conflicts': '2',
        'reviewed_by_expert_panel': '3',
        'practice_guideline': '4'
    }

    print('Reading Clinvar VCF file...')
    counter = 0
    with progressbar.ProgressBar(max_value=row_count) as pbar:
        for line in read_vcf_file(filename):
            if '#' in line:
                # check if it is a valid VCF format
                if not is_vcf_format:
                    print('Checking VCF file validity...')
                    if 'fileformat=VCFv4' in line:
                        is_vcf_format = True
                        print('VCF file is valid.')
                    else:
                        print('Invalid VCF format or supplied file is not VCF file.')
                        print('Terminating application')
                        sys.exit(1)
                else:
                    if 'fileDate' in line:
                        release_date = line.split('=')[1].replace('-', '').strip('\n')
                        print('Data release date found: {0}'.format(release_date))
                                            
                    if 'reference' in line:
                        genome_ref = line.split('=')[1].strip('\n')
                        print('Reference genome assembly version: {0}'.format(genome_ref))
            else:
                # parse each variant line and extract the relevant fields
                counter +=1
                temp_var_buffer = []
                temp_arr = line.strip('\n').split('\t')
                # check if the chromosome field has a valid entry
                if not is_valid_chr(temp_arr[0]):
                    continue 
                temp_var_buffer.append(temp_arr[0])  # chr
                temp_var_buffer.append(temp_arr[1])  # pos
                temp_var_buffer.append(temp_arr[3])  # ref
                temp_var_buffer.append(temp_arr[4])  # alt
                temp_var_buffer.append(temp_arr[2])  # clinvar variant id
                info_fields = {
                    'clin_sig': '.',
                    'review_status': '.',
                    'stars': '0',
                    'disease': '.',
                    'allele_id': '.',
                    'dbsnp_id': '.'
                }
                process_info_field(temp_arr[7], info_fields, clinvar_star_rating)  # info fields
                temp_var_buffer.append(info_fields['clin_sig'])
                temp_var_buffer.append(info_fields['review_status'])
                temp_var_buffer.append(info_fields['stars'])
                temp_var_buffer.append(info_fields['disease'])
                temp_var_buffer.append(info_fields['allele_id'])
                temp_var_buffer.append(info_fields['dbsnp_id'])
                # add to output buffer
                out_buffer.append('\t'.join(temp_var_buffer))
                pbar.update(counter)
    print('Completed reading Clinvar VCF file.')
    return genome_ref, release_date


def process_info_field(info_field, retinfo, clinvar_star_rating):
    temp_info_arr = info_field.split(';')
    for field in temp_info_arr:
        if 'ALLELEID=' in field:
            retinfo['allele_id'] = field.split('=')[1]
        elif 'CLNDN=' in field:
            retinfo['disease'] = field.split('=')[1].replace('_', ' ')
        elif 'CLNSIG=' in field:
            retinfo['clin_sig'] = field.split('=')[1]
        elif 'RS=' in field:
            retinfo['dbsnp_id'] = field.split('=')[1]
        elif 'CLNREVSTAT=' in field:
            retinfo['review_status'] = field.split('=')[1]
            retinfo['stars'] = clinvar_star_rating[retinfo['review_status']]


def get_gzip_app():
    gz_app = ''
    if sys.platform.startswith("linux"):
        gz_app = "zcat"
    elif sys.platform.startswith("darwin"):
        gz_app = "gunzip -c"
    else:
        pass
    return gz_app


# helper function to parse the vcf file and stream the contents
def read_vcf_file(filename):
    with gzip.open(filename, 'rb') as vcf_file:
        for line in vcf_file:
            yield line


def is_valid_chr(chrom):
    valid_chrom_list = [str(x) for x in range(1,23)]
    valid_chrom_list.extend(['X', 'Y', 'MT'])
    valid_chrom_list.extend(['chr' + str(x) for x in range(1,23)])
    valid_chrom_list.extend(['chrX', 'chrY', 'chrMT'])
    return chrom in valid_chrom_list
    

if __name__ == '__main__':
    process_clinvar_db(sys.argv[1:])
