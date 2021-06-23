task raiss {

    String docker
    String chrom
    String pop
    String name
    String gwas = "COVID_" + name
    
    Array[File] z_scores

    String ref_panel_chr_pop
    String ref_panel_chr=sub(ref_panel_chr_pop,"POP",pop)
    File ref=sub(ref_panel_chr,"CHROM",chrom)
    
    String ref_panel_path_pop
    String ref_panel_path=sub(ref_panel_path_pop,"POP",pop)
        
    
    String ld_slices_path_pop
    String ld_slices_path=sub(ld_slices_path_pop,"POP",pop)

    String ld_slices_chr_pop
    String ld_slices_chr=sub(ld_slices_chr_pop,"POP",pop)
    File ld_slices_file=sub(ld_slices_chr,"CHROM",chrom)
    Array[Array[File]] ld_slices=read_tsv(ld_slices_file)
    
    Float eigen_threshold
    Float r2_threshold

    command <<<

        echo "COVID-19 HGI meta-analysis - RAISS sumstats imputation"

        mkdir z_scores
        mkdir z_scores_imputed
        
        find . -name "z_${gwas}_${chrom}.txt" -exec cp {} ./z_scores/ \;
              
        raiss --chrom ${chrom} \
        --gwas ${gwas} \
        --ref-folder /cromwell_root/${ref_panel_path} \
        --ld-folder /cromwell_root/${ld_slices_path} \
        --zscore-folder z_scores \
        --output-folder z_scores_imputed \
        --ld-type scipy \
        --ref-panel-suffix .bim \
        --eigen-threshold ${eigen_threshold} \
        --R2-threshold ${r2_threshold}

        echo "`date` done"
    >>>

    output {
        File out_imp_scores = "/cromwell_root/z_scores_imputed/z_" + gwas + "_" + chrom + ".txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}


task combine {

    String docker
    File sumstat_lift
    String out_prefix
    Array[File] results
    
    command <<<
        
        echo "`date` Combining imputed Z-scores"     
        time cat \
        <(head -n 1 ${results[0]}) \
        <(for file in ${sep=" " results}; do tail -n+2 $file; done) \
        | bgzip > ${out_prefix}.sumstats_imputed_zscores.txt.gz

        echo "`date` Combining original and imputed sumstats"
        time Rscript /sumstats_imputation/scripts/combine_zscores.R \
        -s ${sumstat_lift} \
        -i ${out_prefix}.sumstats_imputed_zscores.txt.gz \
        -o ${out_prefix}.sumstats_imputed.txt

        bgzip ${out_prefix}.sumstats_imputed.txt
    
        echo "`date` Compressing.."
        time bgzip ${out_prefix}.sumstats_imputed.txt
    
    >>>

    output {
        File z_scores_imputed = "${out_prefix}.sumstats_imputed_zscores.txt.gz"
        File out = "${out_prefix}.sumstats_imputed.txt.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "16 GB"
        disks: "local-disk 50 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}


task plot {

    File sumstat
    String base = basename(sumstat)
    String docker
    Int loglog_ylim

    command <<<

        mv ${sumstat} /cromwell_root/${base}

        Rscript /sumstats_imputation/scripts/qqplot_sumstats_imp.R --file ${base} --bp_col "POS" --chrcol "#CHR" --pval_col "p.value" --snp_col "rsID" --loglog_ylim ${loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: 10*ceil(size(sumstat, "G")) + " GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}



workflow raiss_combine {

    Array[String] chrom_list
    File sumstat_lift
    Array[String] sumstat_file
    Array[File] z_scores

    scatter (chrom in chrom_list){ 
        call raiss {
            input:
            chrom = chrom,
            name = sumstat_file[3],
            z_scores = z_scores,
            pop=sumstat_file[1]
        } 
    }

    call combine {
       input:
       sumstat_lift = sumstat_lift,
       out_prefix = sumstat_file[3],
       results = raiss.out_imp_scores
    }

    call plot {
       input:
       sumstat = combine.out
    }

    output {
       File out = combine.out
    }
}