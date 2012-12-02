from math import floor, ceil
from pandas import read_csv

def read_bed(bed_file):
    peaks = read_csv(bed_file, sep="\t", header=None)

    peaks.columns = ['chromosome', 'start', 'end', 'name', 'score']
    peaks = peaks.set_index('name')

    return peaks

def clip_to_fit_resolution(peaks, resolution=1):

    if resolution == 1:
        return peaks

    peaks['start'] = peaks['start'].map(lambda x : resolution * int(floor(float(x) / resolution)) )
    peaks['end']   = peaks['end'].map(lambda x : resolution * int(ceil(float(x) / resolution)) )

    return peaks

