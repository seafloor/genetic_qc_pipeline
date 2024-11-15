import pandas as pd
import argparse

if __name__ == '__main__':
    parser.setup...take infile,n_sd,outfile
    args = parser.parse_args()

    df = pd.read_csv(infile, sep='\t')

    outliers = df.loc[df['F'] < np.abs(df['F'].mean() + (n_sd * df['F'].sd())), ['#FID', 'IID']]

    outliers.to_csv(outfile, sep='\t')

