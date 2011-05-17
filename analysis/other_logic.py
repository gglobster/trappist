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
