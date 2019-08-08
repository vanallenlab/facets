workflow facets_workflow {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File? pileup_snp_vcf = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/common_sites/00-common_all.vcf.gz"
    File? pileup_snp_vcf_index = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/common_sites/00-common_all.vcf.gz"

    Int? min_map_quality = 15
    Int? min_base_quality = 20
    Int? min_read_counts_normal = 25
    Int? min_read_counts_tumor = 0
    Int? pseudo_snps = 100

    Int? cval = 150
    Int? maxiter = 10
    Int? seed_initial = 42
    Int? seed_iterations = 10

    Int? pileup_memoryGB = 3
    Int? pileup_diskGB = 200

	Int? facets_memoryGB = 3
	Int? facets_diskGB = 200

    Int? inferwgd_memoryGB = 3
	Int? inferwgd_diskGB = 50

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
            memoryGB=pileup_memoryGB,
            diskGB=pileup_diskGB,
            preemptible_attempts=preemptible_attempts
    }

    call FACETS {
        input:
            pair_name=pair_name,
            pileup=Pileup.pileup,
            ndepth=min_read_counts_normal,
            cval=cval,
            maxiter=maxiter,
            seed_initial=seed_initial,
            seed_iterations=seed_iterations,
            memoryGB=facets_memoryGB,
            diskGB=facets_diskGB,
            preemptible_attempts=preemptible_attempts
    }

    call InferWGD {
        input:
            cncf=FACETS.cncf,
            memoryGB=inferwgd_memoryGB,
            diskGB=inferwgd_diskGB,
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
        String facets_purity = FACETS.purity
        String facets_ploidy = FACETS.ploidy
        String facets_log_likelihood = FACETS.log_likelihood
        String facets_dip_log_r = FACETS.dip_log_r
        String facets_seed_used = FACETS.seed_used
        String facets_n_iterations_na_purity = FACETS.n_iterations_na_purity
        String facets_fraction_mcn_ge2 = InferWGD.fraction_mcn_ge2
        String facets_wgd_bool = InferWGD.wgd_bool
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

    Int? memoryGB
    Int? diskGB
    Int? preemptible_attempts

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
    Int cval
    Int maxiter
    Int seed_initial
    Int seed_iterations

    Int? memoryGB
    Int? diskGB
    Int? preemptible_attempts

    command <<<
        Rscript /facets.R ${pair_name} ${pileup} ${ndepth} ${cval} ${maxiter} ${seed_initial} ${seed_iterations}
    >>>

    runtime {
        docker: "vanallenlab/facets:v0.5.14.1"
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
        String purity = read_string("purity.txt")
        String ploidy = read_string("ploidy.txt")
        String log_likelihood = read_string("log_likelihood.txt")
        String dip_log_r = read_string("dip_log_r.txt")
        String seed_used = read_string("seed_used.txt")
        String n_iterations_na_purity = read_string("number_iterations_with_na_purity.txt")
    }
}

task InferWGD {
    File cncf

    Int? memoryGB
    Int? diskGB
    Int? preemptible_attempts

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
