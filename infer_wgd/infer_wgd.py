import pandas as pd
import argparse


TOTAL_COPY_NUMBER = 'tcn.em'
MINOR_COPY_NUMBER = 'lcn.em'
MAJOR_COPY_NUMBER = 'mcn.em'

MCN_THRESHOLD = 2.0
PERCENT_THRESHOLD = 0.50


def calculate_fraction_mcn_greater_than_threshold(df):
    df[MINOR_COPY_NUMBER] = df[MINOR_COPY_NUMBER].fillna(0.0)
    df[MAJOR_COPY_NUMBER] = df[TOTAL_COPY_NUMBER].subtract(df[MINOR_COPY_NUMBER])
    number_segments_mcn_greater_than_threshold = df[df[MAJOR_COPY_NUMBER].ge(MCN_THRESHOLD)].shape[0]
    return number_segments_mcn_greater_than_threshold / df.shape[0]


def read_cncf(handle):
    return pd.read_csv(handle, sep='\t', comment='#')


def return_wgd_bool(percent):
    if percent > PERCENT_THRESHOLD:
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
    fraction = calculate_fraction_mcn_greater_than_threshold(cncf)
    wgd = return_wgd_bool(fraction)

    write_output(str(round(fraction, 4)), 'facets_fraction_mcn_ge2.txt')
    write_output(str(wgd), 'facets_wgd_boolean.txt')
