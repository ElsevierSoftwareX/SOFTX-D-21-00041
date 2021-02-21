#!/usr/bin/env python3
from __future__ import print_function
from __future__ import division

import sys,os
from glob import glob

def replace_header(header, fname, write=True,
                   top='/**', bot='*/'):

    print(fname)
    
    # get license
    file1 = open(header, 'r') 
    header_lines = file1.readlines() 
    
    # get file
    file1 = open(fname, 'r') 
    lines = file1.readlines() 

    # clean end of files
    tmp = list()
    for line in lines:
        tmp.append(line.rstrip())
    lines = tmp

    # find start and end of header
    first = None
    last = None
    brief = "TODO"
    for i,line in enumerate(lines):
        if not first and line.startswith(top):
            first = i
        elif not last and line.endswith(bot):
            last = i
        elif '@brief' in line:
            if line.split()[-1] != '@brief':
                brief = line.split()[-1]
        else:
            pass
    #print(first,last)

    # if first is None, there is no header
    if first is not None and last is not None:
        lines = lines[last+1:]
    #print('first new line:',lines[0])

    # add licence
    nlines = header_lines + lines

    # replace info
    tmp = list()
    for line in nlines:
        l = line.rstrip()
        l = l.replace("__FILENAME__",fname.split('/')[-1])
        l = l.replace("__BRIEF__",brief)
        tmp.append(l)
    nlines = tmp

    # writing to file
    if write:
        with open(fname, 'w') as file2:
            for line in nlines:
                file2.write(line + os.linesep)

# -----------------------------------------
def get_all_files(extension):
    files = glob('./**/*'+extension, recursive=True)
    files = [f for f in files if './build' not in f] # remove build files
    return files

# -----------------------------------------
def update():
    fnames = get_all_files('.hh') + get_all_files('.cc') + get_all_files('.hh.in')
    for fname in fnames:
        replace_header('.headers_for_c.txt',fname)

    fnames = get_all_files('CMakeLists.txt') + get_all_files('cmake.in')
    for fname in fnames:
        replace_header('.headers_for_cmake.txt',fname,top='###',bot='####')

        
# -----------------------------------------
if __name__ == '__main__':
    usage = """./replace_header.py"""
    if len(sys.argv) != 1:
        sys.exit(usage)

    update()
