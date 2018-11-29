from oasis import oasisAR1
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='input to oasis wrapper')
parser.add_argument('y_file')
parser.add_argument('gam', type = float)
parser.add_argument('theta', type = float)
parser.add_argument('type', type = str)
parser.add_argument('out_file')

args = parser.parse_args()

# read in tmp data from y_file

import os; print(os.getcwd())

print("trying to open file %s" % args.y_file)
abs_file_path = args.y_file

df = pd.read_csv(abs_file_path, sep=',',header=None)
y = df.values.flatten()


# wraps the oasis ar1 function for easy calls from R
if args.type == 'penalized':
    fit = oasisAR1(y, args.gam, lam = args.theta)
else:
    fit = oasisAR1(y, args.gam, s_min = args.theta)

df = pd.DataFrame([fit[0], fit[1]]).transpose()
df.to_csv(args.out_file, index = False, header = ['calcium', 'spikes'])