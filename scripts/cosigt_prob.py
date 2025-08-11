#!/usr/bin/env python

#
# std import
#
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit

#
# third party packages
#
from sklearn.linear_model import HuberRegressor
import numpy as np
import pandas as pd

N_HEADERS = 4

def get_new_content_stats(s):
    M = np.zeros(s.shape[0]-1)
    M[1:] = s.loc[s.index[1]:s.index[-2]]
    return s.loc[s.index[1]:] - M

def compute_tettelin_params(s):
    model = HuberRegressor()
    logy = s.map(np.log)

    params = [np.nan, np.nan]
    if not logy.isna().any() and not np.isinf(logy).any():
        reg = model.fit(s.index.map(np.log).to_frame(), logy)
        params = [np.exp(reg.intercept_), -reg.coef_[0]]
    return pd.Series(index=['K', 'α'], data=params)

def construct_probability_table(df):
    df = df.copy(deep=True)
    df = df.droplevel([0], axis=1)
    df.drop(0, inplace=True)
    df.index -= 1
    df_new = df.apply(get_new_content_stats)
    df_params = df_new.apply(compute_tettelin_params)
    df_res = pd.concat([df_params.iloc[1, :].map(lambda a: df.index[-1]**(1-a)).to_frame(name='P(X>x)').T, df_params])
    df_res.index.name = ''
    df_res.columns.names = df.columns.names
    return df_res


if __name__ == '__main__':
    description='''
    Computes the probability that a countable (bp, node, edge) is not part of the pangenome.
    '''
    parser = ArgumentParser(formatter_class=ADHF, description=description)
    parser.add_argument('histogram', type=open, help='Histogram table computed by panacus')

    out = stdout
    args = parser.parse_args()

    df = pd.read_csv(args.histogram, sep='\t', header=list(range(N_HEADERS)), index_col=[0], comment='#')
    df = df.loc[:, df.columns.map(lambda x: x[0] == 'growth')]
    if not df.shape[1]:
        print('This script requires growth columns computed by panacus, but none seems to be in the input table.', file=stderr)
        exit(1)

    # ignore warnings
    np.seterr(divide='ignore', invalid='ignore')

    # do the magic!
    df_res = construct_probability_table(df)
    df_res.to_csv(stdout, sep='\t')
