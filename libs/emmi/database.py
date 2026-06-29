import glob
import os
#import numpy as np

def convert_string(s: str):
    if s.isdigit():
        return int(s)
    elif s.replace('.', '', 1).isdigit() and s.count('.') == 1:
        return float(s)
    else:
        return s

def build_database(pattern):
    database = []
    files = glob.glob(pattern)
    for f in sorted(files):
        abspath = os.path.abspath(f)
        basename = os.path.basename(f)
        tag, ext = os.path.splitext(basename)
        if ext != '.tif':
            continue
        print(f' --- new database entry: {tag}')
        dbentry = { 'path':abspath }
        fields = tag.split('_')
        for field in fields:
            if '=' not in field:
                continue
            key, value = field.split('=', 1)
            dbentry[key] = convert_string(value)
        database.append(dbentry)
    return database

def build_coordinates(database):
    x_vals = sorted({d['x'] for d in database if 'x' in d}, reverse=True) # if 'x' in d})
    y_vals = sorted({d['y'] for d in database if 'y' in d}, reverse=True) # if 'y' in d})
    coords = [ [ (y, x) for x in x_vals ] for y in y_vals ]
    return coords

def get_filename(database, conditions):
    for d in database:
        if all(d[k] == v for k, v in conditions.items()):
            return d['path']
    return None
