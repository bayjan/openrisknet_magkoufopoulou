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

rule all:
    input:
        compound_info, transcriptomics_data

rule line_counts:
    output:
        out_file
    input:
        in_file
    shell:
        "wc -l {input} > {output}"

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

