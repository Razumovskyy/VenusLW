from struct import unpack
import numpy as np

def to_format(filename):
    fmt = '%.10f', '%d'
    value, channel = np.loadtxt(filename, unpack=True)
    np.savetxt(filename+'_formatted', np.transpose((value, channel)), fmt=fmt, delimiter='    ')
    print('goog')
    return True

file1 = 'LENTA.__2'
file2 = 'LENTA_BORIS.__2'

to_format(file1)
to_format(file2)

