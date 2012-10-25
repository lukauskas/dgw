'''
Created on 16 Oct 2012

@author: saulius
'''
import os
import gzip

from functools import wraps
from peak import Peak


def memoised(func):
    
    cache = {}
    
    @wraps(func)
    def memoise_response(*args):
        try:
            return cache[(func.__name__, args)]
        except KeyError:
            ans = func(*args)
            cache[(func.__name__, args)] = ans
            return ans
    
    return memoise_response
        
class FileDescription(object):
    
    basename = None
    
    def __init__(self, basename, **kwargs):
        self.basename = basename
        
        for key, value in kwargs.iteritems():
            setattr(self, key, value)
        

    @classmethod
    def from_description_file_line(cls, line):
        line = line.split('\t')
        
        basename = line[0]
        parameters = line[1]
        
        parameters = parameters.split('; ')
        
        kwargs = {}
        for p in parameters:

            p = p.split('=')
            kwargs[p[0]] = p[1]
        
        return cls(basename, **kwargs)

def open_gzipped_file(filename):
    f = gzip.open(filename, 'rb')
    return f
   
def find_description_file_from_filename(filename):
    return os.path.join(os.path.dirname(filename), 'files.txt')

def process_peak_file(filename):
    dt = generate_description_table_from_file(find_description_file_from_filename(filename))
    desc = dt[os.path.basename(filename)]
    track_name = desc.tableName
    
    if desc.type == 'broadPeak':
        if 'broadPeak' in filename:
            line_parsing_function = Peak.parse_from_broadpeak_data_row
        elif 'narrowPeak' in filename:
            line_parsing_function = Peak.parse_from_narrowpeak_data_row
        else:
            raise Exception, 'Unknown filetype'
    
    f = open_gzipped_file(filename)
    
    for line in f:
        yield line_parsing_function(track_name, line)
    
    f.close()
        
@memoised
def generate_description_table_from_file(filename):
    f = open(filename, 'r')
    descriptions = {}
    
    for line in f:
        file_description = FileDescription.from_description_file_line(line)
        descriptions[file_description.basename] = file_description
    
    return descriptions

    