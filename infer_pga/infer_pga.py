import pandas as pd
import argparse


TOTAL_COPY_NUMBER = 'tcn.em'
MINOR_COPY_NUMBER = 'lcn.em'
MAJOR_COPY_NUMBER = 'mcn.em'
CHROM = 'chrom'
START_POS = 'start'
END_POS = 'end'
BASES = 'bases'


def calculate_diploid_CN(chrom):
    return 2 if chrom <= 22 else 1


def calculate_bases_covered(start, end):
    return end.subtract(start)


def calculate_fraction_genome_not_diploid(df):
    # uses total copy number !=2 for autosomes and != 1 for sex chromosomes
    df[TOTAL_COPY_NUMBER] = df[TOTAL_COPY_NUMBER].fillna(0.0)
    diploid_CNs = df[CHROM].apply(calculate_diploid_CN)

    idx_not_diploid = df[df[TOTAL_COPY_NUMBER].astype(float) != diploid_CNs].index

    df[BASES] = calculate_bases_covered(df.loc[:, START_POS], df.loc[:, END_POS])
    total_bases = df[BASES].sum()
    segment_percent = df[BASES].divide(total_bases)

    return segment_percent.loc[idx_not_diploid].sum()


def read_cncf(handle):
    return pd.read_csv(handle, sep='\t', comment='#')


def write_output(value, filename):
    f = open(filename, "w")
    f.write(value)


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser(prog='Infer percent genome altered from FACETS',
                                         description='Infer PGA from FACETS based on total copy number.')
    arg_parser.add_argument('--cncf', help='CNCF dataframe from FACETS', required=True)
    args = arg_parser.parse_args()

    cncf = read_cncf(args.cncf)
    fraction = calculate_fraction_genome_not_diploid(cncf)

    write_output(str(round(fraction, 4)), 'facets_fraction_genome_CNA_altered.txt')
