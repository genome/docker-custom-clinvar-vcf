import unidecode, sys

'''
Take in ReviewStatus and map it to a numerical score, described in
https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/otto/clinvar/clinVarToBed#L222
'''
def status_to_score(status):
    status = status.strip().lower().replace(" ", "")
    score = None

    if status == "nointerpretationforthesinglevariant":
        score = 0
    if status=="criteriaprovided,conflictinginterpretations":
        score = 1
    if status=="criteriaprovided,multiplesubmitters,noconflicts":
        score = 2
    if status=="criteriaprovided,singlesubmitter":
        score = 1
    if status == "practiceguideline":
        score = 4
    if status == "reviewedbyexpertpanel":
        score = 3
    if "noassertion" in status or status=="-":
        score = 0

    #TODO should handle this more gracefully?
    if score is None:
        assert(False)

    return score

'''
Take in dictionary with keys corresponding to fields that should be added to the vcf's INFO column
and return a string that can be directly added, in the form key=value;key=value...
'''
def info_dict_to_string(info_dict):
    fields = [ field + "=" + unidecode.unidecode(unicode(str(value).replace(";", ",").replace(" ", "_"), 'utf-8')) for field,value in info_dict.items() ]
    return ";".join(fields)

'''
BEGIN main portion of script
'''

summary_tsv = sys.argv[1]
clinvar_vcf = sys.argv[2]

info_additions = {}
with open(summary_tsv) as ncbi_tab_delimited:
    lines = ncbi_tab_delimited.readlines()
    del lines[0] #delete header
    for line in lines:
        line = line.strip()
        line = line.split("\t")

        #TODO handle special character \0 ?
        #TODO handle malformed chromosome names?

        if line[16] != "GRCh38":
            continue
        alleleID = int(line[0]) #to be used as dictionary key
        review_status = line[24] #parsed and used to generate another value for final vcf

        info_fields = {}
        info_fields["SCORE"] = status_to_score(review_status)
        info_fields["RCVACC"] = line[11] #RCVaccession in ncbi data
        info_fields["TESTEDINGTR"] = line[27] #TestedinGTR in ncbi data
        info_fields["PHENOTYPELIST"] = line[13] #PhenotypeList in ncbi data
        info_fields["NUMSUBMIT"] = line[25] #NumberSubmitters in ncbi data
        
        #add in the missing value indicator, as specified in http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
        if line[26] == "":
            info_fields["GUIDELINES"] = "."
        else:
            info_fields["GUIDELINES"] = line[26]

        info_additions[alleleID] = info_dict_to_string(info_fields)

with open(clinvar_vcf) as ncbi_vcf:
    all_lines = ncbi_vcf.readlines()
    header = all_lines[:28]
    data = all_lines[28:]
    for line in data:
        line = line.split("\t")
        raw_info = line[7].split(";") #key/value pairs in INFO column are separated by ;
        info_pairs = [ x.split("=") for x in raw_info ]
        info_list = []
        allele_id = None
        for pair in info_pairs:
            if pair[0] == "ALLELEID":
                allele_id = int(pair[1])
            elif pair[0] == "CLNDISDB":
                info_list.append("PHENOTYPE=" + pair[1])
            elif pair[0] == "CLNSIG":
                info_list.append("CLINSIGN=" + pair[1])

        if allele_id is None:
            assert(False) #TODO handle this better?

        final_info = ";".join(info_list) 
        if allele_id not in info_additions:
            continue
        if final_info != "":
            final_info += ";"
        final_info = final_info + info_additions[allele_id]
        output_list = line[:7]
        output_list.append(final_info)
        output_line = "\t".join(output_list)
        print output_line
