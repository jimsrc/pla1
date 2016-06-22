import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', type=int)
parser.add_argument('-f', type=float)
parser.add_argument('--IDs', type=str)
parser.add_argument('--file', type=file)

try:
    pa = parser.parse_args()
    ids = pa.IDs
    mylist = map(int, ids.split(','))
    print pa
except IOError, msg:
    parser.error(str(msg))
