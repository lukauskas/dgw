from collections import namedtuple
from itertools import izip
Point = namedtuple('Point', ["location", "type"])

    
class Peak(object):
    
    chromosome = None
    __start      = None
    __end        = None
    __resolution = None
    
    __data     = None
    __points_of_interest = None
    
    def __init__(self, chromosome, start, end, resolution=1):
        self.chromosome = chromosome
        self.__start      = int(start)
        self.__end        = int(end)
        self.__resolution = resolution
        
        self.data = []
        self.points_of_interest = set()

    @property
    def points_of_interest(self):
        return self.__points_of_interest

    @points_of_interest.setter
    def points_of_interest(self, value):
        self.__points_of_interest = value

    def add_point_of_interest(self, point):
        self.points_of_interest.add(point)
    
    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        '''
        End of the peak. Please note that the end is not inclusive.
        
        That is Peak(chr1, 0, 100) will have 100 elements in it that are integers in [0,99], but
        will not contain the bp=100
        '''
        return self.__end

    @start.setter
    def start(self, value):
        self.__start = int(value)

    @end.setter
    def end(self, value):
        self.__end = int(value)

    @property
    def data(self):
        return self.__data
    @data.setter
    def data(self, value):
        self.__data = value
    
    @property
    def resolution(self):
        return self.__resolution
    
    def __data_bin_of_point(self, point):
        adjusted_point = point - self.start
        bin_id = adjusted_point / self.resolution
        return bin_id
        
    def expanded_data_relative_to_point(self, point):
        ans = []
        for pos in xrange(self.start, self.end):
            bin_of_pos = self.__data_bin_of_point(pos)
            n = self.data[bin_of_pos]
            adjusted_pos = pos - point
            ans.append( (adjusted_pos, n) )
        
        return ans
    
    def data_relative_to_point(self, point):
        point_bin = self.__data_bin_of_point(point)
        
        ans = []
        for i, n in enumerate(self.data):
            rel_bin = i - point_bin
            ans.append((rel_bin, n))
        
        return ans
            
    def data_relative_to_start(self):
        return list(enumerate(self.data))
    
    def expanded_data_relative_to_start(self):
        return self.expanded_data_relative_to_point(self.start)
    
    def find_interesting_points_from_set(self, interesting_points_set):
        for pos in xrange(self.start, self.end):
            if (self.chromosome, pos) in interesting_points_set:
                self.add_point_of_interest(pos)
        
    
    def __repr__(self):
        return '<Peak({0!r}, {1!r}, {2!r})>'.format(self.chromosome, 
                                             self.start, 
                                             self.end) 
    
