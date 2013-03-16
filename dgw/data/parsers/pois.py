import numpy as np

def map_to_bins(poi, other_region, resolution=50, ignore_non_overlaps=False, account_for_strand_information=False):
    """
    Maps the points of interest to the bins in the other region provided.
    Note that this does not validate the chromosome, this needs to be done by the caller

    :param poi:
    :param other_region:
    :param resolution:
    :param ignore_non_overlaps:
    :param account_for_strand_information:
    :return:
    """

    other_start = other_region['start']
    other_end = other_region['end']

    poi = np.asarray(poi, dtype=int)

    if np.min(poi) < other_start or np.max(poi) >= other_end:
        if ignore_non_overlaps:
            return []
        else:
            raise ValueError('Points of interest do not overlap with regions of interest.')
    resolution = int(resolution)
    n_bins = (other_end - other_start) / resolution
    poi_offsetted = poi - other_start
    poi_bins = poi_offsetted / resolution  # This should be integer division

    # Remove duplicates
    poi_bins = np.unique(poi_bins)

    if not account_for_strand_information or other_region['strand'] == '+':
        return poi_bins
    else:
        return sorted((n_bins - 1) - poi_bins)

def from_simple(poi_file, main_regions, resolution=50, ignore_non_overlaps=False, account_for_strand_information=False):
    """
    Parses points of interest from the specified file in the simple form:
    ```index:poi_1,poi_2,poi_3,...,poi_n```, e.g.
    ```MACS_peak_1:23,45,100```

    :param poi_file: file with pois in simple format
    :param main_regions: the main regions the pois will be parsed from
    :return:
    """

    if isinstance(poi_file, basestring):
        poi_file = open(poi_file)
        need_closing = True
    else:
        need_closing = False

    try:
        poi_dict = {}

        for line in poi_file:
            try:
                ix, pois = line.split(':')
            except ValueError:
                raise ValueError('Not a valid poi file in simple format')

            pois = map(int, pois.replace(' ', '').split(','))

            try:
                main_region = main_regions.ix[ix]
            except KeyError:
                # Some of the regions in POI might not exist in main regions, e.g. due to --random-sample constraint ...
                # that's okay
                continue

            bins = map_to_bins(pois, main_region, resolution=resolution, ignore_non_overlaps=ignore_non_overlaps,
                               account_for_strand_information=account_for_strand_information)

            if len(bins) > 0:
                poi_dict[ix] = bins
        return poi_dict
    finally:
        if need_closing:
            poi_file.close()



