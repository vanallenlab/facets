workflow facets_workflow {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    # VCF of common SNPs. Please use the a suitable VCF that is consistent with the genome build
    # Read more in the pileup documentation: https://github.com/mskcc/facets/tree/master/inst/extcode
    File pileup_snp_vcf
    File pileup_snp_vcf_index

    Int? min_map_quality = 15
    Int? min_base_quality = 20
    Int? min_read_counts_normal = 25
    Int? min_read_counts_tumor = 0
    Int? pseudo_snps = 100

    Int? preproc_cval = 150
    Int? proc_cval = 150
    Int? maxiter = 10
    Int? seed_initial = 42
    Int? seed_iterations = 10

    Int? preemptible_attempts = 3

    call Pileup {
        input:
            pair_name=pair_name,
            tumor_bam=tumor_bam,
            tumor_bam_index=tumor_bam_index,
            normal_bam=normal_bam,
            normal_bam_index=normal_bam_index,
            snp_vcf=pileup_snp_vcf,
            snp_vcf_index=pileup_snp_vcf_index,
            min_map_quality=min_map_quality,
            min_base_quality=min_base_quality,
            min_read_counts_normal=min_read_counts_normal,
            min_read_counts_tumor=min_read_counts_tumor,
            pseudo_snps=pseudo_snps,
            preemptible_attempts=preemptible_attempts
    }

    call FACETS {
        input:
            pair_name=pair_name,
            pileup=Pileup.pileup,
            ndepth=min_read_counts_normal,
            preproc_cval=preproc_cval,
            proc_cval=proc_cval,
            maxiter=maxiter,
            seed_initial=seed_initial,
            seed_iterations=seed_iterations,
            preemptible_attempts=preemptible_attempts
    }

    call InferWGD {
        input:
            cncf=FACETS.cncf,
            preemptible_attempts=preemptible_attempts
    }

    call InferPGA {
        input:
            cncf=FACETS.cncf,
            preemptible_attempts=preemptible_attempts
    }

    output {
        File facets_pileup = Pileup.pileup
        File facets_genome_segments = FACETS.genome_segments
        File facets_diagnostics_plot = FACETS.diagnostics_plot
        File facets_iterations_plot = FACETS.iterations_plot
        File facets_cncf = FACETS.cncf
        File facets_output = FACETS.summary
        File facets_flags = FACETS.flags
        File facets_emflags = FACETS.emflags
        File facets_iterations = FACETS.iterations
        File facets_rdata = FACETS.rdata
        String facets_purity = FACETS.purity
        String facets_ploidy = FACETS.ploidy
        String facets_log_likelihood = FACETS.log_likelihood
        String facets_dip_log_r = FACETS.dip_log_r
        String facets_seed_used = FACETS.seed_used
        String facets_n_iterations_na_purity = FACETS.n_iterations_na_purity
        String facets_n_segments = FACETS.n_segments
        String facets_n_segments_na_lcn = FACETS.n_segments_na_lcn
        String facets_fraction_mcn_ge2 = InferWGD.fraction_mcn_ge2
        String facets_wgd_bool = InferWGD.wgd_bool
        String facets_percent_genome_altered = InferPGA.pga
    }
}

task Pileup {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File? snp_vcf
    File? snp_vcf_index

    Int? min_map_quality
    Int? min_base_quality
    Int? min_read_counts_normal
    Int? min_read_counts_tumor
    Int? pseudo_snps

    Int? preemptible_attempts
    Int? memoryGB = 4
    Int? diskGB = ceil(1.1 * (size(normal_bam, "G") + size(tumor_bam, "G"))) + 20

    command <<<
        /./snp-pileup --verbose \
            --gzip \
            --min-map-quality ${min_map_quality} \
            --min-base-quality ${min_base_quality} \
            --min-read-counts ${min_read_counts_normal},${min_read_counts_tumor} \
            --pseudo-snps ${pseudo_snps} \
            ${snp_vcf} \
            ${pair_name}.pileup \
            ${normal_bam} \
            ${tumor_bam}
    >>>

    runtime {
        docker: "vanallenlab/facets:v0.5.14.1"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
        preemptible: preemptible_attempts
    }

    output {
        File pileup="${pair_name}.pileup.gz"
    }
}

task FACETS {
    String pair_name
    File pileup
    Int ndepth
    Int preproc_cval
    Int proc_cval
    Int maxiter
    Int seed_initial
    Int seed_iterations

    Int? preemptible_attempts
    Int? memoryGB = 4
    Int? diskGB = ceil(1.1 * (size(pileup, "G"))) + 20

    command <<<
        Rscript /facets.R ${pair_name} ${pileup} ${ndepth} ${preproc_cval} ${proc_cval} ${maxiter} ${seed_initial} ${seed_iterations}
    >>>

    runtime {
        docker: "vanallenlab/facets:v0.5.14-2"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
        preemptible: preemptible_attempts
    }

    output {
        File genome_segments = "${pair_name}.genome_segments.pdf"
        File diagnostics_plot = "${pair_name}.diagnostic_plot.pdf"
        File iterations_plot = "${pair_name}.facets_iterations.pdf"
        File cncf = "${pair_name}.facets_cncf.txt"
        File summary = "${pair_name}.facets_output.txt"
        File flags = "${pair_name}.facets_flags.txt"
        File emflags = "${pair_name}.facets_emflags.txt"
        File iterations = "${pair_name}.facets_iterations.txt"
        File rdata = "${pair_name}.RData"
        String purity = read_string("purity.txt")
        String ploidy = read_string("ploidy.txt")
        String log_likelihood = read_string("log_likelihood.txt")
        String dip_log_r = read_string("dip_log_r.txt")
        String seed_used = read_string("seed_used.txt")
        String n_iterations_na_purity = read_string("number_iterations_with_na_purity.txt")
        String n_segments = read_string("number_segments.txt")
        String n_segments_na_lcn = read_string("number_segments_NA_LCN.txt")
    }
}

task InferWGD {
    File cncf

    Int? preemptible_attempts
    Int? memoryGB = 4
    Int? diskGB = ceil(1.1 * (size(cncf, "G"))) + 5

    command <<<
        python /infer_wgd.py --cncf ${cncf}
    >>>

    runtime {
        docker: "vanallenlab/facets:infer_wgd"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
        preemptible: preemptible_attempts
    }

    output {
        String fraction_mcn_ge2 = read_string("facets_fraction_mcn_ge2.txt")
        String wgd_bool = read_string("facets_wgd_boolean.txt")
    }
}

task InferPGA {
    File cncf

    Int? preemptible_attempts
    Int? memoryGB = 4
    Int? diskGB = ceil(1.1 * (size(cncf, "G"))) + 5

    command <<<
        python /infer_pga.py --cncf ${cncf}
    >>>

    runtime {
        docker: "vanallenlab/facets:infer_pga"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
        preemptible: preemptible_attempts
    }

    output {
        String pga = read_string("facets_fraction_genome_CNA_altered.txt")
    }
}
