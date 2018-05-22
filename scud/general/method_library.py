#!/usr/bin/env cctbx.python

from __future__ import division
import numpy as np
import sys

def change_element_in_json(json_file = None,
                           label_name = None,
                           new_element = None,
                           l = None):
    """
    Load a json file and change one of the elements in it
    """
    import json
    if l != None:
        l.process_message('Renaming json element')
    with open(json_file,'r') as infile:
        tmp = json.load(infile)
    tmp[label_name] = new_element
    json.dumps(tmp,sort_keys=True)
    with open(json_file,'w') as outfile:
        json.dump(tmp,outfile)

def datetime_string():
    '''
    Return a string of current date in form of YYYYMMDD
'''
    import datetime as dt
    return dt.date.today().isoformat().replace('-','')

def input_mtz(ls):
    '''
    Returns list of MTZ files
'''
    mtz_list = []
    for f in ls:
        if f[-4:] == '.mtz': mtz_list.append(f)
    return mtz_list

def input_pdb(ls):
    '''
    Returns list of PDB files
'''
    pdb_list = []
    for f in ls:
        if f[-4:] == '.pdb': pdb_list.append(f)
    return pdb_list

def print_title(title):
    '''
    Fancy way of printing Title
'''
    print('#########################\n#\t%s\t#\n#########################'%title)

def vecDistSquare(vec1,vec2):
    '''
    Flex method to calculate distance between cart lists
'''
    from scitbx.array_family import flex
    dvec = flex.sqrt((vec1-vec2).dot())
    return dvec**2
    
