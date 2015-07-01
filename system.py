from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class System(WESTSystem):
    """
    System for P53 folding and unfolding.
    """

    def initialize(self):
        """
        Initializes system
        """
        self.pcoord_ndim  = 3
        self.pcoord_len   = 6
        self.pcoord_dtype = numpy.float32
        binbounds         = [ 0.00, 2.80, 2.88, 3.00, 3.10, 3.29, 3.79, 3.94,
                              4.12, 4.39, 5.43, 5.90, 6.90, 7.90, 8.90, 9.90,
                             10.90,11.90,12.90,13.90,14.90,15.90,float('inf')]
        self.bin_mapper   = RectilinearBinMapper([binbounds])

        self.bin_target_counts      = numpy.empty((self.bin_mapper.nbins,),
                                        numpy.int)
        self.bin_target_counts[...] = 4

def _2D_pcoord_loader_color_tracker(fieldname, coord_file, segment, single_point=False):
    """
    This function loads a 2-dimensional progress coordinate, performs some logic to track color,
    then returns the 2+1 dimensional progress coordinate to the system to be processed.
    In this tutorial, there are 3 dimensions specified in this file; runseg.sh returns two of them.
    The third is calculated here.
    Note that we are defining our states only based on one progress coordinate dimension, in this example.

    **Arguments:**
        :*fieldname*:      Key at which to store dataset
        :*coord_filename*: Temporary file from which to load coordinates
        :*segment*:        WEST segment
        :*single_point*:   Data to be stored for a single frame
                           (only false half the time)
    """

    # These are the raw coordinates.
    coord_raw = numpy.loadtxt(coord_file, dtype=numpy.float32) 
    # These are the states; they are left inclusive, and right exclusive, which is consistent with the normal
    # binning procedure.
    color_bins = [(0.0,0.20),(11.60,float('inf'))]
    unknown_state = 2
    system = westpa.rc.get_system_driver()

    if single_point == True:
        npts = 1
    else:
        npts = system.pcoord_len

    coords = numpy.empty((npts,2), numpy.float32)
    colors = numpy.empty((npts), numpy.float32)
    if single_point == True:
        colors[:] = unknown_state
        for istate,state_tuple in enumerate(color_bins):
            # Note that here, we are using the first dimension and first dimension alone.
            # The shape of the returned coord_raw is slightly different if single_point evalues
            # to true.
            # We evalulate whether or not we're in a state; if not, we leave it as in the
            # unknown state.
            if coord_raw[0] >= state_tuple[0] and coord_raw[0] < state_tuple[1]:
                colors[:] = istate
        coords[:] = coord_raw[...]
    else:
        # If we're not the first point, we set the state to be what it was in the beginning
        # of the iteration.  We only want to update the state when we update a bin for purposes
        # of state tracking.
        colors[:] = segment.pcoord[0][2]
        coords[:] = coord_raw[...]

    for istate,state_tuple in enumerate(color_bins):
        if coords[-1,0] >= state_tuple[0] and coords[-1,0] < state_tuple[1]:
            colors[-1] = istate
    
    # We require different stacking behavior to return things in the proper order
    # depending on how many points we have.  I could probably clean this up.
    if single_point == True:
        segment.pcoord = numpy.hstack((coords[0,0],coords[0,1],colors[:]))
    else:
        # This could easily be modified to return N dimensions.
        segment.pcoord = numpy.swapaxes(numpy.vstack((coords[:,0],coords[:,1],colors[:])), 0, 1)

def coord_loader(fieldname, coord_filename, segment, single_point=False):
    """
    Loads and stores coordinates

    **Arguments:**
        :*fieldname*:      Key at which to store dataset
        :*coord_filename*: Temporary file from which to load coordinates
        :*segment*:        WEST segment
        :*single_point*:   Data to be stored for a single frame
                           (should always be false)
    """
    # Load coordinates
    n_frames = 6
    n_atoms  = 2
    coord    = numpy.loadtxt(coord_filename, dtype = numpy.float32)
    coord    = numpy.reshape(coord, (n_frames, n_atoms, 3))

    # Save to hdf5
    segment.data[fieldname] = coord

def log_loader(fieldname, log_filename, segment, single_point=False):
    """
    Loads and stores log

    **Arguments:**
        :*fieldname*:    Key at which to store dataset
        :*log_filename*: Temporary file from which to load log
        :*segment*:      WEST segment
        :*single_point*: Data to be stored for a single frame
                         (should always be false)
    """
    # Load log
    with open(log_filename, 'r') as log_file:
        raw_text = [line.strip() for line in log_file.readlines()]

    # Determine number of fields
    n_frames = 6
    n_fields = 0
    line_i   = 0
    starts   = []
    while line_i < len(raw_text):
        line = raw_text[line_i]
        if len(line.split()) > 0:
            start = line.split()[0]
            if start in starts:
                break
            else:
                try:
                    float(start)
                    n_fields += len(line.split())
                except ValueError:
                    starts.append(start)
        line_i += 1
    dataset = numpy.zeros((n_frames, n_fields), numpy.float32)
#    print(dataset.shape, starts)

    # Parse data
    line_i  = 0
    frame_i = 0
    field_i = 0
    while line_i < len(raw_text):
        line = raw_text[line_i]
#        print(line_i, frame_i, field_i, line)
        if len(line.split()) > 0:
            start = line.split()[0]
            try:
                float(start)
            except ValueError:
                starts.append(start)
            if start not in starts:
                for field in line.split():
                    dataset[frame_i, field_i] = float(field)
                    if field_i == n_fields - 1:
                        frame_i += 1
                        field_i  = 0
                    else:
                        field_i += 1
        line_i += 1

    # Save to hdf5
    segment.data[fieldname] = dataset
