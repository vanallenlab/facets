import pandas as pd
import argparse


TOTAL_COPY_NUMBER = 'tcn.em'
MINOR_COPY_NUMBER = 'lcn.em'
MAJOR_COPY_NUMBER = 'mcn.em'
CHROM = 'chrom'
START_POS = 'start'
END_POS = 'end'
BASES = 'bases'
SEGMENT_PERCENT = 'segment_percent'

MCN_THRESHOLD = 2.0
FRACTION_THRESHOLD = 0.50


def subset_autosomal(df):
    return df[df[CHROM].le(22)]


def calculate_bases_covered(start, end):
    return end.subtract(start)


def calculate_fraction_mcn_greater_than_threshold(df):
    # Filling NAs with 1 to be conservative. TO DO: annotate with https://github.com/mskcc/facets/issues/62
    df[MINOR_COPY_NUMBER] = df[MINOR_COPY_NUMBER].fillna(0.0)
    df[MAJOR_COPY_NUMBER] = df[TOTAL_COPY_NUMBER].subtract(df[MINOR_COPY_NUMBER])
    idx_mcn_ge_2 = df[df[MAJOR_COPY_NUMBER].astype(float).ge(MCN_THRESHOLD)].index

    df[BASES] = calculate_bases_covered(df.loc[:, START_POS], df.loc[:, END_POS])
    total_bases = df[BASES].sum()
    segment_percent = df[BASES].divide(total_bases)

    return segment_percent.loc[idx_mcn_ge_2].sum()


def read_cncf(handle):
    return pd.read_csv(handle, sep='\t', comment='#')


def return_wgd_bool(fraction):
    if fraction > FRACTION_THRESHOLD:
        return True
    else:
        return False


def write_output(value, filename):
    f = open(filename, "w")
    f.write(value)


if __name__ == "__main__":
    # Please cite:
    # DOI: 10.1038/s41588-018-0165-1
    # Bielski CM, Zehir A, Penson AV, et al. Genome doubling shapes the evolution and prognosis of advanced cancers.
    # Nat Genet. 2018;50(8):1189-1195.

    arg_parser = argparse.ArgumentParser(prog='Infer WGD from FACETS',
                                         description='Infer WGD from FACETS based on major allele number.')
    arg_parser.add_argument('--cncf', help='CNCF dataframe from FACETS', required=True)
    args = arg_parser.parse_args()

    cncf = read_cncf(args.cncf)
    autosomal = subset_autosomal(cncf)
    fraction = calculate_fraction_mcn_greater_than_threshold(autosomal)
    wgd = return_wgd_bool(fraction)

    write_output(str(round(fraction, 4)), 'facets_fraction_mcn_ge2.txt')
    write_output(str(wgd), 'facets_wgd_boolean.txt')
