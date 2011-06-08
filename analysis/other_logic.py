__author__ = 'GG'

def uniqify(value_list, idfun=None):
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    """Retain only unique elements while preserving original order."""
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in value_list:
        marker = idfun(item[0])
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

def ensure_dir(dir_path):
    """Check that the directory exists; if not, create it."""
    from os import path, makedirs
    abs_path = path.abspath(dir_path)
    if not path.exists(abs_path):
        try: makedirs(abs_path)
        except Exception as message: 
            status = 1
            # TODO: make graceful fail or request input if interactive mode
        else: 
            message = 'created path'
            status = 0
    else: 
        message = 'path exists'
        status = 0
    report = {'message': message, 'status': status}
    return abs_path, report

def create_id(prefix):
    """Create a unique identification number with a given prefix."""
    # If several IDs with the same prefix are created in quick succession or
    # in parallel, some may end up being the same. It would be better to use
    # another method that adds a random element. But it's good to have the
    # date of creation in the ID, for provenance tracking purposes.
    from time import time
    timestamp = str(time()).split('.')
    unique_id = prefix+'_'+timestamp[0]
    return unique_id