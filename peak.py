from collections import namedtuple
from itertools import izip
Point = namedtuple('Point', ["location", "type"])

class Peak:
    
    chromosome = None
    __start      = None
    __end        = None
    
    __data     = None
    __points_of_interest = None
    
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.__start      = int(start)
        self.__end        = int(end)
        
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

    def data_relative_to_point(self, point):
        for pos, n in izip(xrange(self.start, self.end), self.data):
            adjusted_pos = pos - point
            yield (adjusted_pos, n)
    
    def find_interesting_points_from_set(self, interesting_points_set):
        for pos in xrange(self.start, self.end):
            if pos in interesting_points_set:
                self.add_point_of_interest(pos)
        
        
    def __repr__(self):
        return '<Peak({0!r}, {1!r}, {2!r})>'.format(self.chromosome, 
                                             self.start, 
                                             self.end) 
    
