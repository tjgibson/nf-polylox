#!/usr/bin/env python3

import polyloxpgen
import sys


# define input and output files
input_file = sys.argv[1]
sample_name = sys.argv[2]

out_name = sample_name + '_pgen'


# run polyloxpgen
df_pgen = polyloxpgen.polylox_pgen(input_file, "./", out_name)