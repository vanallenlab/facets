workflow facets_workflow {
    call FACETS
}

task FACETS {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File reference_fasta
    File reference_index
    File reference_dict

    Int? niter = 1

    Int? memoryGB = 3
    Int? diskSpaceGB = 100

    command <<<
        # move fasta for perl script
        mv ${reference_fasta} /usr/gitc/Homo_sapiens_assembly19.fasta
        mv ${reference_index} /usr/gitc/Homo_sapiens_assembly19.fasta.fai
        mv ${reference_dict} /usr/gitc/Homo_sapiens_assembly19.dict

        # Remember starting directory b/c FireCloud doesn't let you put outputs of files w/ abs paths
        home_dir=$PWD

        cd /usr/gitc

        # Samtools mpileup wrapper
        perl /usr/gitc/scripts_snp_pileup/snp-pileup.pl ${normal_bam} ${tumor_bam} ${pair_name}

        echo "Finished pileup, launching FACETS"

        # mv for R script
        mv /usr/gitc/${pair_name} /usr/gitc/pileup/${pair_name}

        # Gotta add a little vanilla
        Rscript --vanilla /usr/gitc/facets.R ${niter}

        cd /usr/gitc/pileup

        # move to home directory for output
        mv /usr/gitc/pileup/Facets_output.txt $home_dir/Facets_output.txt
        mv /usr/gitc/pileup/Facets_iterations.txt $home_dir/Facets_iterations.txt
        mv /usr/gitc/pileup/${pair_name} $home_dir/${pair_name}_pileup.txt
        mv /usr/gitc/pileup/genome_segments.pdf $home_dir/genome_segments.pdf
        mv /usr/gitc/pileup/diagnostic_plot.pdf $home_dir/diagnostic_plot.pdf
        mv /usr/gitc/pileup/fit_cncf.txt $home_dir/fit_cncf.txt

        cd $home_dir

        sed -n '2p' $home_dir/Facets_output.txt | awk '{print $3}' > purity.txt
        sed -n '2p' $home_dir/Facets_output.txt | awk '{print $4}' > ploidy.txt
    >>>

    output {
        File facetsOutput = "Facets_output.txt"
        File facetsIterations = "Facets_iterations.txt"
        File pileup = "${pair_name}_pileup.txt"
        File genomeSegmentsPdf = "genome_segments.pdf"
        File diagnosticPlotPdf = "diagnostic_plot.pdf"
        File fitCncfTable = "fit_cncf.txt"
        String purity = read_string("purity.txt")
        String ploidy = read_string("ploidy.txt")
    }

    runtime {
        docker: "jakeconway/facets:with_plots"
        memory: "${memoryGB} GB"
        cpu: "1"
        disks: "local-disk ${diskSpaceGB} HDD"
    }

}