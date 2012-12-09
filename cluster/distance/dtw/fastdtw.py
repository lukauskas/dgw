'''
Created on 13 Nov 2012

@author: saulius
'''
from sys import maxint as MAXINT
import numpy as np

def fast_cityblock(a, b):
    '''
    Similar to scipy.spatial.distance.cityblock but skips the array validations
    @param a:
    @param b:
    @return:
    '''

    return np.abs(a-b).sum()

def traceback_path(x, y, cost_matrix):

    # Trace back the path
    min_cost_path = []

    i = len(x)-1
    j = len(y)-1
    min_cost_path.append((i,j))

    while (i > 0) or (j > 0):
        # Find costs of moving in three possible directions:

        # diagonal
        if i >= 1 and j >=1:
            diag_cost = cost_matrix[i-1,j-1]
        else:
            diag_cost = float('inf')

        # left
        if i >= 1:
            left_cost = cost_matrix[i-1, j]
        else:
            left_cost = float('inf')

        # down
        if j >= 1:
            down_cost = cost_matrix[i, j-1]
        else:
            down_cost = float('inf')

        # determine where to move in. 
        # Prefer moving diagonally or towards i==j axis if there
        # are ties

        if diag_cost <= left_cost and diag_cost <= down_cost:
            i -= 1
            j -= 1
        elif left_cost < diag_cost and left_cost < down_cost:
            i -= 1
        elif down_cost < diag_cost and diag_cost < left_cost:
            j -= 1
        elif i <= j: # left_cost == right_cost < diag_cost
            j -= 1
        else: # left_cost == right_cost < diag_cost
            i -= 1

        min_cost_path.append((i,j))

    min_cost_path.reverse()

    return min_cost_path

def dtw(x, y, distance_function=fast_cityblock, return_aligned_path=False):
    cost_matrix = np.empty([len(x), len(y)])

    for i in range(len(x)):
        for j in range(len(y)):
            local_dist = distance_function(x[i], y[j])
            if i == 0 and j == 0:
                cost_matrix[i,j] = local_dist
            elif i == 0:
                #assert(cost_matrix[i, j-1] is not None)
                cost_matrix[i,j] = cost_matrix[i, j-1] +\
                                   local_dist
            elif j == 0:
                #assert(cost_matrix[i-1, j] is not None)
                cost_matrix[i,j] = cost_matrix[i-1, j] +\
                                   local_dist
            else:

                #assert(cost_matrix[i, j-1] is not None)
                #assert(cost_matrix[i-1, j] is not None)
                #assert(cost_matrix[i-1, j-1] is not None)
                min_global_cost = min(cost_matrix[i-1, j],
                    cost_matrix[i, j-1],
                    cost_matrix[i-1, j-1])

                cost_matrix[i,j] = min_global_cost + local_dist

    # Minimal cost is at the top-right corner of the matrix
    min_cost = cost_matrix[len(x)-1,len(y)-1]
    if not return_aligned_path:
        return min_cost
    else:
        min_cost_path = traceback_path(x, y, cost_matrix)
        return min_cost, min_cost_path

def constrained_dtw(x, y, window, distance_function=fast_cityblock, return_aligned_path=False):
    #assert(isinstance(window, DTWWindow))

    cost_matrix = window.get_cost_matrix()

    for (i, j) in window:
        local_dist = distance_function(x[i], y[j])
        if i == 0 and j == 0:
            cost_matrix[i,j] = local_dist
        elif i == 0:
            #assert(cost_matrix[i, j-1] is not None)
            cost_matrix[i,j] = cost_matrix[i, j-1] +\
                               local_dist
        elif j == 0:
            #assert(cost_matrix[i-1, j] is not None)
            cost_matrix[i,j] = cost_matrix[i-1, j] +\
                               local_dist
        else:

            #assert(cost_matrix[i, j-1] is not None)
            #assert(cost_matrix[i-1, j] is not None)
            #assert(cost_matrix[i-1, j-1] is not None)
            min_global_cost = min(cost_matrix[i-1, j],
                cost_matrix[i, j-1],
                cost_matrix[i-1, j-1])

            cost_matrix[i,j] = min_global_cost + local_dist

    min_cost = cost_matrix[len(x)-1, len(y)-1]
    if not return_aligned_path:
        return min_cost
    else:
        min_cost_path = traceback_path(x, y, cost_matrix)
        return min_cost, min_cost_path

def shrink_time_series(x, shrink_factor):
    '''
    Shrinks time series x by a factor specified as shrink_factor.
    The series is shrunk by averaging components.
    
    Examples:
    >>> shrink_time_series([1,2,3,4,5,6,7,8,9], 2.0)
    [1.5, 4.0, 6.5, 8.5]
    >>> shrink_time_series([1,2,3,4,5,6], 3.0)
    [2.0, 5.0]

    This is direct translation from a similar function in PAA.java in FastDTW implementation.
    
    @param x: time series to be shrunk. should be list like and support direct element access
    @param shrink_factor: ratio to reduce the size by. E.g
    @type shrink_factor: float
    '''

    #assert(shrink_factor > 1)

    original_size = len(x)
    # Get the size of shrunk time series
    reduced_size = int(float(original_size) / shrink_factor)

    # Determine size of a sampled font (might be fractional)
    reduced_point_size = float(original_size) / reduced_size

    # Keep track of points being averaged
    read_from = 0
    read_to   = None

    ans = []
    aggregate_sizes = []
    while read_from < original_size:

        read_to = int(round(reduced_point_size*(len(ans)+1)) -1 )

        points_to_read = read_to - read_from +1

        val_sum = 0

        for i in range(read_from, read_to+1):
            val_sum += x[i]

        val_sum = float(val_sum) / points_to_read
        ans.append(val_sum)
        aggregate_sizes.append(points_to_read)
        read_from = read_to + 1

    return ans, aggregate_sizes

def expanded_res_window(x, y, shrunk_x, shrunk_y,
                        agg_x, agg_y,
                        low_res_path, radius):

    if radius > 0:
        raise Exception, 'Not implemented yet'

    # variables to keep track of current location of higher-res 
    # projected path

    # TODO: that looks fishy: refers to the min values of low_res_path
    # but uses them for the larger seq
    current_i = low_res_path[0][0]
    current_j = low_res_path[0][1]

    last_warped_i = MAXINT
    last_warped_j = MAXINT

    window = DTWWindow(len(x), len(y))

    for warped_i, warped_j in low_res_path:
        block_size_i = agg_x[warped_i]
        block_size_j = agg_y[warped_j]

        # if the path moved up or diagonally, the next cell's values
        # on the x axis will be larger
        if warped_j > last_warped_j:
            current_j += agg_y[last_warped_j]

        if warped_i > last_warped_i:
            current_i += agg_x[last_warped_i]

        # If diagonal move was performed add 2 cells to the edges 
        # to create a continuous path with even with:
        #                        |_|_|x|x|     then mark      |_|_|x|x|
        #    ex: projected path: |_|_|x|x|  --2 more cells->  |_|X|x|x|
        #                        |x|x|_|_|        (X's)       |x|x|X|_|
        #                        |x|x|_|_|                    |x|x|_|_| 
        if (warped_j > last_warped_j) and (warped_i > last_warped_i):
            window.mark_visited(current_i - 1, current_j)
            window.mark_visited(current_i, current_j - 1)

        # fill in cells that are created by projection from cell
        # in a lower resolution to higher resolution
        for ii in xrange(0, block_size_i):
            window.mark_visited(current_i+ii, current_j)
            window.mark_visited(current_i+ii, current_j+block_size_j-1)

        last_warped_i = warped_i
        last_warped_j = warped_j

    # Expanding window goes here
    # TODO: implement. Currently I assume radius = 0
    return window

class DTWWindow(object):

    _rows    = None
    _columns = None

    _min_values = None
    _max_values = None

    _size = None
    def __init__(self, columns, rows):
        self._rows = rows
        self._columns = columns

        self._min_values = [None] * columns
        self._max_values = [None] * columns
        self._size = 0

    @property
    def columns(self):
        return self._columns

    @property
    def rows(self):
        return self._rows

    @property
    def size(self):
        return self._size

    @property
    def boundaries(self):
        '''
           Returns min_values and max values for each column
        @return:
        '''
        return zip(self._min_values, self._max_values)

    def max_value_for(self, i):
        return self._max_values[i]

    def min_value_for(self, i):
        return self._min_values[i]

    @property
    def min_column(self):
        return 0

    @property
    def min_row(self):
        return 0

    def mark_visited(self, column, row):
        #assert(column >= 0)
        #assert(column <= self.columns)
        if self._min_values[column] is None:
            self._min_values[column] = row
            self._max_values[column] = row
            self._size += 1
        elif self._min_values[column] > row:
            self._size += self._min_values[column] - row
            self._min_values[column] = row
        elif self._max_values[column] < row:
            self._size += row - self._max_values[column]
            self._max_values[column] = row

    def __iter__(self):
        for i in xrange(self.min_column, self.columns):
            for j in xrange(self.min_value_for(i), self.max_value_for(i)+1):
                yield (i, j)


    def get_cost_matrix(self):
        return WindowMatrix(self)

class WindowMatrix(object):
    _window = None

    __values = None
    INFINITY = float('inf')
    def __init__(self, window):
        '''
        
        @param window:
        @type window: DTWWindow
        '''
        self._window = window

        self.__init_cell_values()

    def __init_cell_values(self):
        values = {}
        window_boundaries = self._window.boundaries

        for col, (min_row, max_row) in enumerate(window_boundaries):
            for row in xrange(min_row, max_row+1):
                values[(col, row)] = None

        self.__values = values

    def contains(self, col, row):
        return self._window.min_value_for(col) <= row <= self._window.max_value_for(col)

    def __getitem__(self, key):

        try:
            return self.__values[key]
        except KeyError:
            return self.INFINITY

    def __setitem__(self, key, value):
        try:
            self.__values[key] = value
        except KeyError:
            raise IndexError, 'Cannot set {0} as matrix does not contain this entry'.format(key)

def fast_dtw(x, y, distance_function=fast_cityblock, return_aligned_path = False):

    radius = 0

    min_ts_size = radius + 2

    len_x = len(x)
    len_y = len(y)
    if len_x <= min_ts_size or len_y <= min_ts_size:
        return dtw(x, y, distance_function, return_aligned_path=return_aligned_path)
    else:
        SHRINK_FACTOR = 2.0
        shrunk_x, agg_x = shrink_time_series(x, SHRINK_FACTOR)
        #assert(sum(agg_x) == len(x))
        shrunk_y, agg_y = shrink_time_series(y, SHRINK_FACTOR)
        #assert(sum(agg_y) == len(y))

        _, low_res_path = fast_dtw(shrunk_x, shrunk_y, distance_function, return_aligned_path=True)

        window = expanded_res_window(x, y,
            shrunk_x, shrunk_y,
            agg_x, agg_y,
            low_res_path,
            radius)

        return constrained_dtw(x, y, window, distance_function=distance_function,
            return_aligned_path=return_aligned_path)
