from collections import namedtuple
from itertools import izip
Point = namedtuple('Point', ["location", "type"])

class Peak:
    
    chromosome = None
    __start      = None
    __end        = None
    
    __data     = None
    __points_of_interest = None
    
    data = property(get_data, set_data, None, None)
    start = property(get_start, set_start, None, None)
    end = property(get_end, set_end, None, None)
    points_of_interest = property(get_points_of_interest, set_points_of_interest, None, None) 
    
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start      = start
        self.end        = end
        
        self.data = []
        self.points_of_interest = set()

    def get_points_of_interest(self):
        return self.__points_of_interest


    def set_points_of_interest(self, value):
        self.__points_of_interest = value

    def add_point_of_interest(self, point):
        self.points_of_interest.add(point)

    def get_start(self):
        return self.__start


    def get_end(self):
        return self.__end


    def set_start(self, value):
        self.__start = int(value)


    def set_end(self, value):
        self.__end = int(value)


    def get_data(self):
        return self.__data

    def set_data(self, value):
        self.__data = value

    def get_data_relative_to_point(self, point):
        location = point.location
        
        for pos, n in izip(xrange(self.start, self.end), self.data):
            adjusted_pos = pos - location
            yield (adjusted_pos, n)
    
    def find_interesting_points_from_set(self, interesting_points_set, points_type):
        for pos in xrange(self.start, self.end):
            if pos in interesting_points_set:
                self.add_point_of_interest(Point(location=pos, type=points_type))
        
        
    def __repr__(self):
        return '<Peak({0!r}, {1!r}, {2!r})>'.format(self.chromosome, 
                                             self.start, 
                                             self.end)    
    
