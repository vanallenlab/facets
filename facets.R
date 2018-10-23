library('ggplot2')
library('facets')

args = commandArgs(trailingOnly=TRUE)
pair_name = args[1]
pileup = args[2]
cval = as.numeric(args[3])
maxiter = as.numeric(args[4])

rcmat = readSnpMatrix(pileup)
xx = preProcSample(rcmat)
oo = procSample(xx, cval = cval)
fit = emcncf(oo, maxiter = maxiter)

purity = as.numeric(signif(fit$purity, 3))
ploidy = as.numeric(signif(fit$ploidy, 3))
log_likelihood = fit$loglik
cncf = fit$cncf
dip_log_r = oo$dipLogR
flags = oo$flags
emflags = fit$emflags

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
