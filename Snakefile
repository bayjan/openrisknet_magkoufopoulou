in_file = "Snakefile"
out_file = "line_counts.txt"

# Followings are inputs
# transcriptomics_data_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28878&format=file"
transcriptomics_data_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28878/matrix/GSE28878_series_matrix.txt.gz"
compound_info_url = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/carcin/33/7/10.1093_carcin_bgs182/2/bgs182_Supplementary_Data.zip?Expires=1526557830&Signature=1fvEwy7yzkeqHPknk~X82ScOyQmXDQonDM1UpOdp~y~rQUBThynLmfKCQymB-4~hplt0S~UGbYcTkCdwm0bpJgNxbG2Pmnk1A6oIq9mw-4QLu-5OEqzHEmOcKFqYXrt7UocogQ0~0aLYeKDt6-o3aN8~xLS5rWK2wrYso6VUZ~vJWd80BE4F7BFn5Wib4ToCbXE4HJFa10PXFJGV-cOnk3zAgSY5kK-YOiO2SuqEg22JE-hizmIczYC1As5u3BPu~5nGogMIJ5genwNOIWsAYO-bdeOO~OM3vv5BfnsXDj8Ec5-MFL5CShupQDX16xQuxXDBoS9RZZv2qSfErz-Czg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"

# Followings are outputs from different steps
transcriptomics_data = "GSE28878_series_matrix_good_format.txt"
compound_info = "Supplementary_Data_1.tsv"
compound_info_zip = "supplementary_data.zip"
compound_info_excel = "Supplementary_Data_1.xls"
training_data_compound_info = "training_data_compound_info.tsv"
validation_data_compound_info = "validation_data_compound_info.tsv"


rule all:
    input:
        compound_info, transcriptomics_data, training_data_compound_info, validation_data_compound_info


#####
# Data collection and preparation part STARTS here
#####

rule collect_data:
    input:
        transcriptomics_data, compound_info
    shell:
        "echo 'Collected transcriptomics and compound info data'"

rule get_transcriptomics_data:
    output:
        transcriptomics_data
    shell:
        """
        curl -X GET "{transcriptomics_data_url}" > transcriptomics_data_tmp1.gz
        gunzip transcriptomics_data_tmp1.gz
        echo 'BEGIN{{start=0;}};$1 ~ /series_matrix_table_end/{{start=0;}}; {{ if(start==1){{print;}} }}; $1 ~ /series_matrix_table_begin/{{start=1;}};' > parse_transcript_data.awk
        awk -f parse_transcript_data.awk  transcriptomics_data_tmp1|tr -d '"' > {transcriptomics_data}
        """

rule get_compound_info:
    """This command retrieves the archive file & then only extracts the Excel file we need"""
    output:
        compound_info_excel
    shell:
        'curl -X GET "{compound_info_url}" > {compound_info_zip} && unzip {compound_info_zip} {compound_info_excel}'

rule compound_info_tsv:
    output:
        compound_info
    input:
        compound_info_excel
    run:
        import pandas as pd
        print("input  is:" + str(input))
        print("output is:" + str(output))
        excel_file = pd.read_excel(io=str(input))
        excel_file.to_csv(path_or_buf=str(output), sep="\t")

rule training_compound_info:
    output:
        training_data_compound_info
    input:
        compound_info
    shell:
        """
        cut -f1,10- {input} |awk 'BEGIN{{print("compound\tgenotoxicity");}};NR>3'|head -35|sed -re 's/\+$/GTX/g; s/\-$/NGTX/g; s/-//g; s/\]//g; s/\[//g' > {output}
        """

rule validation_compound_info:
    output:
        validation_data_compound_info
    input:
        compound_info
    shell:
        """
        cut -f1,10- {input} |awk 'NR>41'|sed -re 's/\+$/GTX/g; s/\-$/NGTX/g;  s/[[:punct:]]//g'|awk 'BEGIN{{print("compound\tgenotoxicity");}}{{print}}' > {output}
        """
        

rule find_corresponding_series:
    """Each series has different solvent, so find correct solvents"""
    shell:
        """
        paste <(grep Sample_title GSE28878_series_matrix.txt|cut -f2-|tr -d '"'|tr '\t' '\n') <(grep Series_sample_id GSE28878_series_matrix.txt|cut -f2-|tr -d '"'|sed -re 's/\s*$//'|tr ' ' '\n')|grep 24h|sed -re 's/^Serie\s*//g; s/, HepG2 exposed to\s*/\t/g; s/for 24h, biological rep\s*/\t/g'|awk 'BEGIN{print("series_id\tcompound\treplicate\tarray_name");}{print}' > solvent_to_exposure.tsv
        """

rule calculate_log2_ratio:
    """Calculate the correct log2ratio using the corresponding solvent for each replicate"""

#####
# Data collection and preparation part ENDS here
#####

#####
# Analysis part STARTS here
#####

# Accuracy is calculated differently
rule calculate_accuracy:

#####
# Analysis part ENDS here
#####
