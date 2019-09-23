library('ggplot2')
library('facets')

args = commandArgs(trailingOnly=TRUE)
pair_name = args[1]
pileup = args[2]
min_normal_depth = as.numeric(args[3])
cval = as.numeric(args[4])
maxiter = as.numeric(args[5])
seed_initial = as.numeric(args[6])
seed_iterations = as.numeric(args[7])
genome_build = "hg19"

if (seed_iterations <= 0) {
    seed_iterations = 1
}

run_facets = function(seed, pileup, min_normal_depth, cval, genome_build, maxiter) {
    set.seed(seed)

    rcmat = readSnpMatrix(pileup)
    xx = preProcSample(rcmat, ndepth = min_normal_depth, cval = cval, gbuild = genome_build)
    oo = procSample(xx, cval = cval)
    fit = emcncf(oo, maxiter = maxiter)

    return(list(seed_oo=oo, seed_fit=fit))
}

plot_facets_iterations = function(pair_name, seeds_dataframe, median_purity, median_ploidy) {
    title = paste('Stability of FACETS outputs across ', seed_iterations, ' seeds,\n', pair_name, sep='')
    p = ggplot(seeds_dataframe, aes(x=ploidy, y=purity)) + xlim(0, 8) + ylim(0, 1) +
        geom_hline(yintercept=median_purity, linetype='dashed', size=1, alpha=0.2) +
        geom_vline(xintercept=median_ploidy, linetype='dashed', size=1, alpha=0.2) +
        geom_jitter(size = 6, color='#E69F00', alpha=0.5) +
        ggtitle(title) +
        theme(plot.title = element_text(size=16, hjust=0.5)) +
        theme(axis.title.x = element_text(size=16)) +
        theme(axis.title.y = element_text(size=16))

    return(p)
}

seeds = seed_initial:(seed_initial+seed_iterations-1)
seeds_dataframe = data.frame()
for (seed in seeds) {
    seed_list = run_facets(seed, pileup, min_normal_depth, cval, genome_build, maxiter)
    seed_oo = seed_list$seed_oo
    seed_fit = seed_list$seed_fit

    dip_log_r = seed_oo$dipLogR
    purity = as.numeric(signif(seed_fit$purity, 3))
    ploidy = as.numeric(signif(seed_fit$ploidy, 3))

    seed_dataframe = data.frame(seed, purity, ploidy, dip_log_r)
    names(seed_dataframe) <- c("seed", "purity", "ploidy", "dip_log_r")
    seeds_dataframe = rbind(seeds_dataframe, seed_dataframe)
}

seeds_dataframe_filename = paste(pair_name, '.facets_iterations.txt', sep='')
write.table(seeds_dataframe, seeds_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)

idx_na_purity = is.na(seeds_dataframe$purity)
n_iteratons_na_purity = sum(as.numeric(idx_na_purity))
median_purity = median(seeds_dataframe[!idx_na_purity, ]$purity, na.rm=TRUE)
median_ploidy = median(seeds_dataframe[!idx_na_purity, ]$ploidy, na.rm=TRUE)

plot_iterations_filename = paste(pair_name, '.facets_iterations.pdf', sep='')
pdf(plot_iterations_filename, height=10, width=7.5)
plot_facets_iterations(pair_name, seeds_dataframe, median_purity, median_ploidy)
dev.off()

if (!is.na(median_purity)) {
  delta_median_purity = abs(seeds_dataframe$purity - median_purity)
  min_delta_purity = min(delta_median_purity, na.rm=TRUE)
  median_purity_seed = as.numeric(seeds_dataframe[delta_median_purity %in% min_delta_purity,]$seed)
  if (length(median_purity_seed) > 1) {median_purity_seed = median_purity_seed[1]}
  used_seed = as.numeric(median_purity_seed)
} else {
  used_seed = seed_initial
}

facets_list = run_facets(used_seed, pileup, min_normal_depth, cval, genome_build, maxiter)
oo = facets_list$seed_oo
fit = facets_list$seed_fit

purity = as.numeric(signif(fit$purity, 3))
ploidy = as.numeric(signif(fit$ploidy, 3))
log_likelihood = fit$loglik
cncf = fit$cncf
dip_log_r = oo$dipLogR
flags = oo$flags
emflags = fit$emflags

number_segments = NROW(cncf)
number_segments_NA_LCN = sum(is.na(cncf$lcn.em))

genome_segments_filename = paste(pair_name, '.genome_segments.pdf', sep='')
diagnostic_plot_filename = paste(pair_name, '.diagnostic_plot.pdf', sep='')
cncf_dataframe_filename = paste(pair_name, '.facets_cncf.txt', sep='')
summary_dataframe_filename = paste(pair_name, '.facets_output.txt', sep='')
flags_filename = paste(pair_name, '.facets_flags.txt', sep='')
emflags_filename = paste(pair_name, '.facets_emflags.txt', sep='')

if(!is.na(purity)) {
    pdf(genome_segments_filename, height=10, width=7.5)
    plotSample(x=oo, emfit=fit)
    dev.off()

    pdf(diagnostic_plot_filename, height=10, width=7.5)
    logRlogORspider(oo$out, oo$dipLogR)
    dev.off()
} else {
    write.table(NA, genome_segments_filename)
    write.table(NA, diagnostic_plot_filename)
}

df = data.frame(pair_name=pair_name, purity=purity, ploidy=ploidy, digLogR=dip_log_r)
df_flags = data.frame(flags=flags)
df_emflags = data.frame(emflags=emflags)

write.table(cncf, cncf_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df, summary_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df_flags, flags_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df_emflags, emflags_filename, sep='\t', quote=FALSE, row.names=FALSE)

write(purity, 'purity.txt')
write(ploidy, 'ploidy.txt')
write(log_likelihood, 'log_likelihood.txt')
write(dip_log_r, 'dip_log_r.txt')
write(used_seed, 'seed_used.txt')
write(n_iteratons_na_purity, 'number_iterations_with_na_purity.txt')
write(number_segments, 'number_segments.txt')
write(number_segments_NA_LCN, 'number_segments_NA_LCN.txt')

# save RData object needed for add_ccf_to_maf_config method, part of Phylogic preprocessing
out = oo
save(fit, out, file=paste(pair_name, '.RData', sep=''))
