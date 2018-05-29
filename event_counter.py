#!/usr/bin/env python

from I3Tray import *
from icecube import icetray, dataclasses, dataio
    
def main(file):
    tray = I3Tray()

    global count
    count = 0 
    
    tray.AddModule("I3Reader","reader", FilenameList = [file])
    
    def event_counter(frame):
    	global count
#        if "LineFit" not in frame.keys():
#            print count #frame["LineFit"]
    	count += 1
    
    tray.AddModule(event_counter, 'counter', Streams = [icetray.I3Frame.Physics])
    
    tray.AddModule("TrashCan", "cleanup")
    
    tray.Execute()
    
    tray.Finish()	
    
    del tray
    print count
    return count

if __name__ == '__main__':
    # test1.py executed as script
    # do something
    import os
    import sys
    import shutil
    from os.path import expandvars
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-i", dest="data_file", help="Input data file.")
    (options, args) = parser.parse_args()

    file = "%s" % (options.data_file)	

    main(file)
