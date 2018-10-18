workflow facets_workflow {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File? pileup_snp_vcf = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/facets/exac_hets.vcf"
    File? pileup_snp_vcf_index = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/facets/exac_hets.vcf.idx"
    Int? pileup_min_map_quality = 15
    Int? pileup_min_base_quality = 20
    Int? pileup_min_read_counts = 0

    Int? cval = 150

    call Pileup {
        input:
            pair_name=pair_name,
            tumor_bam=tumor_bam,
            tumor_bam_index=tumor_bam_index,
            normal_bam=normal_bam,
            normal_bam_index=normal_bam_index,
            snp_vcf=pileup_snp_vcf,
            snp_vcf_index=pileup_snp_vcf_index,
            min_map_quality=pileup_min_map_quality,
            min_base_quality=pileup_min_base_quality,
            min_read_counts=pileup_min_read_counts
    }

    call FACETS {
        input:
            pair_name=pair_name,
            pileup=Pileup.pileup,
            cval=cval
    }

    output {
        File facets_pileup = Pileup.pileup
        File facets_genome_segments = FACETS.genome_segments
        File facets_diagnostics_plot = FACETS.diagnostics_plot
        File facets_cncf = FACETS.cncf
        File facets_output = FACETS.summary
        File facets_flags = FACETS.flags
        File facets_emflags = FACETS.emflags
        String facets_purity = FACETS.purity
        String facets_ploidy = FACETS.ploidy
        String facets_log_likelihood = FACETS.log_likelihood
        String facets_dip_log_r = FACETS.dip_log_r
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
    Int? min_read_counts

    Int memoryGB = 3
    Int diskGB = 200

    command <<<
        /./snp-pileup --verbose \
            --min-map-quality ${min_map_quality} \
            --min-base-quality ${min_base_quality} \
            --min-read-counts ${min_read_counts} \
            ${snp_vcf} \
            ${pair_name}.pileup \
            ${normal_bam} \
            ${tumor_bam}
    >>>

    runtime {
        docker: "vanallenlab/facets:v0.5.14"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
    }

    output {
        File pileup="${pair_name}.pileup"
    }
}

task FACETS {
    String pair_name
    File pileup
    Int cval

    Int memoryGB = 3
    Int diskGB = 200

    command <<<
        Rscript /facets.R ${pair_name} ${pileup} ${cval}
    >>>

    runtime {
        docker: "vanallenlab/facets:v0.5.14"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} SSD"
    }

    output {
        File genome_segments = "${pair_name}.genome_segments.pdf"
        File diagnostics_plot = "${pair_name}.diagnostic_plot.pdf"
        File cncf = "${pair_name}.facets_cncf.txt"
        File summary = "${pair_name}.facets_output.txt"
        File flags = "{pair_name}.facets_flags.txt"
        File emflags = "{pair_name}.facets_emflags.txt"
        String purity = read_string("purity.txt")
        String ploidy = read_string("ploidy.txt")
        String log_likelihood = read_string("log_likelihood.txt")
        String dip_log_r = read_string("dip_log_r.txt")
    }
}
