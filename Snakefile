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
log2ratio_results = "log2ratio_GSE28878_series.tsv"
solvent2exposure = "solvent_to_exposure.tsv"
solvent2exposure_mapping = "solvent_to_exposure_with_solvent_column.tsv"
compound_array_genotoxicity = "compound_array_genotoxicity_info.tsv"

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
        a=$(tempfile -d .)
        cut -f1,10- {input} |awk 'BEGIN{{print("compound\tgenotoxicity");}};NR>3'|head -35 > $a;
        sed -re 's/\+$/GTX/g; s/\-$/NGTX/g; s/-//g; s/\]//g; s/\[//g' $a > {output}
        rm $a
        """

rule validation_compound_info:
    output:
        validation_data_compound_info
    input:
        compound_info
    shell:
        """
        a=$(tempfile -d .)
        cut -f1,10- {input} |awk 'NR>41'|sed -re 's/\+$/GTX/g; s/\-$/NGTX/g;  s/[[:punct:]]//g' > $a;
        awk 'BEGIN{{print("compound\tgenotoxicity");}}{{print}}' $a| sed -re 's/ppDDT\t/DDT\t/g; s/γHCH\t/HCH\t/g; s/\s+/\t/g'> {output}
        rm $a
        """
        

rule find_corresponding_series:
    """Each series has different solvent, so find correct solvents"""
    output:
        solvent2exposure    
    shell:
        """
        a=$(tempfile -d .)
        b=$(tempfile -d .)
        c=$(tempfile -d .)
        d=$(tempfile -d .)
        e=$(tempfile -d .)
        paste <(grep Sample_title transcriptomics_data_tmp1|cut -f2-|tr -d '"'|tr '\t' '\n') <(grep Series_sample_id transcriptomics_data_tmp1 > $a;
        cut -f2- $a|tr -d '"'|sed -re 's/\s*$//'|tr ' ' '\n')|grep 24h > $b;
        sed -re 's/^Serie\s*//g; s/, HepG2 exposed to\s*/\t/g; s/for 24h, biological rep\s*/\t/g' $b > $c;
        awk 'BEGIN{{print("series_id\tcompound\treplicate\tarray_name");}}{{print}}' $c|sed -re 's/\s+/\t/g' > $d; 
        sed -re 's/DEPH/DEHP/g; s/Ethyl\t/EtAc\t/g; s/NPD\t/NDP\t/g; s/Paracres\t/pCres\t/g; s/Phenol\t/Ph\t/g; s/Resor/RR/g' $d > $e;
        sed -re 's/2-Cl\t/2Cl\t/g' $e> {output}
        rm $a $b $c $d $e
        """

rule compound_array_mapping:
    """Creates a file that stores a mapping between genotoxicity, compound and array information"""
    output:
        compound_array_genotoxicity
    input:
        validation_data_compound_info, training_data_compound_info, solvent2exposure
    shell:
        """
        a=$(tempfile -d .)
        echo -e "series_id\tcompound\treplicate\tarray_name\tgenotoxicity" > {output}; 
        LANG=en_EN join -i -o 2.1,2.2,2.3,2.4,2.5,1.2 -t $'\t' -1 1 -2 2 <(cat {validation_data_compound_info} {training_data_compound_info}|sort -k1.1i,1.3i  -t $'\t' ) <(LANG=en_EN sort -fbi -t $'\t' -k 2 {solvent2exposure}) > $a;
        grep -iv genotoxicity $a >> {output};
        rm $a;
        """

rule calculate_log2_ratio:
    """Calculate the correct log2ratio using the corresponding solvent for each replicate, values are already in log scale, so just subtract solvent values"""
    output:
        log2ratio_results, solvent2exposure_mapping
    input:
        transcriptomics_data, solvent2exposure    
    run:
        import pandas as pd
        print("transcriptomics file is: " + input[0])
        print("treatment2solvent file is: " + input[1])
        transcr_df = pd.read_table(filepath_or_buffer=input[0])
        solvent2exposure_df = pd.read_table(filepath_or_buffer = input[1])
        # Find corresponding solvent ids for each compound
        results_df = pd.DataFrame()
        for val in solvent2exposure_df.series_id.drop_duplicates().values:
            tmp = solvent2exposure_df.loc[solvent2exposure_df.series_id==val,:].query("compound.str.lower() in ['dmso','etoh','pbs']")
            tmp.index=range(tmp.shape[0])          
            tmp.columns = ['solvent_'+str(i) for i in tmp.columns]
            compounds = solvent2exposure_df.loc[solvent2exposure_df.series_id==val,:].query("compound.str.lower() not in ['dmso','etoh','pbs']")['compound'].drop_duplicates()
            for compound in compounds:
                tmp1 = solvent2exposure_df.loc[solvent2exposure_df.series_id==val,:].query("compound=='"+str(compound)+"'")
                tmp1.index = range(tmp1.shape[0])                
                tmp2 = pd.concat([tmp1,tmp],axis=1, join='inner')
                if results_df.shape[0] == 0:
                    results_df = tmp2
                else:
                    results_df = results_df.append(tmp2)
        results_df.to_csv(path_or_buf=str(output[1]), sep="\t")
        # Calculate log2ratio
        log2ratio_df = pd.DataFrame()
        for compound in results_df['compound'].drop_duplicates().values:
            for replicate in results_df[results_df['compound']==compound].replicate:
                compound_array = results_df[results_df['compound']==compound].query('replicate=='+str(replicate)).array_name.values[0]
                solvent_array  = results_df[results_df['compound']==compound].query('replicate=='+str(replicate)).solvent_array_name.values[0]
                tmp3 = transcr_df.loc[:,compound_array] - transcr_df.loc[:,solvent_array]
                tmp3.index = transcr_df.index
                tmp3.columns = [compound_array]
                if log2ratio_df.shape[0] == 0:
                    log2ratio_df = tmp3
                    log2ratio_df.columns = tmp3.columns
                    log2ratio_df.index = tmp3.index
                else:
                    column_names = list(log2ratio_df.columns)
                    column_names.append(compound_array)
                    log2ratio_df = pd.concat([log2ratio_df,tmp3],axis=1)
                    log2ratio_df.columns = column_names
        column_names = list(log2ratio_df.columns)
        column_names.insert(0,'ID_REF')
        log2ratio_df = pd.concat([transcr_df['ID_REF'],log2ratio_df], axis=1)
        log2ratio_df.columns = column_names
        log2ratio_df.to_csv(path_or_buf=str(output[0]), sep="\t", index=False)


#####
# Data collection and preparation part ENDS here
#####

#####
# Analysis part STARTS here
#####
rule t_test:
    """Leave one out t-test is carried for every compound, by leaving out all replicates of a compound"""

# Accuracy is calculated differently
rule calculate_accuracy:

#####
# Analysis part ENDS here
#####
