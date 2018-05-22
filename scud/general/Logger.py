#!/usr/bin/env cctbx.python

import os

"""
Class that helps writing log files, this to save outputs
that would normally go to the screen
"""

class Logger(object):

    def __init__(self,filename='log.txt'):
        """
        Initialize Logger class, open file, print header
        """
        self.filename = filename
        self.log = self._openLog()
        self._writeHeader()
        
    def _openLog():
        return open(self.filename,'w')

    def _writeHeader():
        print >> self.log, "HEADERINFO HERE"
