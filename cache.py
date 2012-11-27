'''
Created on 14 Nov 2012

@author: saulius
'''
import hashlib
import pickle
import os
from functools import wraps

CACHE_DIR = '_cache'
CACHE_DISABLED = os.environ.get('DISABLE_CACHING', False) != False
if CACHE_DISABLED:
    print "Cache is disabled"
    
def __cache_filename(func, args, kwargs, unique_hash=None):
    if unique_hash is None:
        important_info = (func.__name__, args, sorted(kwargs.items()))
        hash_str = pickle.dumps(important_info)
        unique_hash = hashlib.md5(hash_str).hexdigest()
  
    filename = '{0}-{1}'.format(func.__name__, unique_hash)
    return os.path.join(CACHE_DIR, filename)
    
def cached(func, unique_hash=None):
    
    @wraps(func)
    def f(*args, **kwargs):
        if CACHE_DISABLED:
            return func(*args, **kwargs)
        
        cache_filename = __cache_filename(func, args, kwargs)
        
        try:
            f = open(cache_filename, 'rb')
            return pickle.load(f)
        
        except IOError:
            result = func(*args, **kwargs)
            
            f = open(cache_filename, 'wb')
            pickle.dump(result, f)
            f.close()
            
            return result
    
    return f

        
    