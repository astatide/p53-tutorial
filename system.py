from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
import westpa
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
        self.pcoord_ndim  = 2
        self.pcoord_len   = 11
        self.pcoord_dtype = numpy.float32
        # As the RMSD coordinate is taken relative to the coil, aligned on the coil,
        # it will remain sensitive to coil changes.  It's best to assume the maximum is
        # not dissimilar to the maximum for the distance; something around 57 A, as
        # that would take into account the peptide flipping completely around.
        # However, we must bin much finer.
        self.rmsd_binbounds         = [0.0+0.4*i for i in xrange(0,19)] + \
                                      [8.0+0.8*i for i in xrange(0,19)] + \
                                      [24.0+11.0*i for i in xrange(0,3)] + [float('inf')]
        # As this is the end to end length of the P53 peptide, the bins should cover
        # everything from the coil to the fully extended peptide, or 3.8 A * 15 or so.
        # (Including caps probably overestimates length, but that's alright)
        # Total of 57 angstroms.
        #self.dist_binbounds         = [0.0+10.0*i for i in xrange(0,6)] + [float('inf')]

        # It's best not to place these at the integer boundaries, due to 
        # oddities with the way numpy/h5py stores the values inside the west.h5 file.
        # Given that we are starting in the coil conformation, the 'unknown state'
        # (that is, 1.5 to float, or 2) will never be used; our bins will never be more
        # than 66% filled.

        self.color_binbounds = [-0.5,0.5,1.5,float('inf')]

        # A simple rectilinear binmapper, with the third dimension as color, to ensure good sampling.
        self.bin_mapper   = RectilinearBinMapper([self.rmsd_binbounds, self.color_binbounds])

        self.bin_target_counts      = numpy.empty((self.bin_mapper.nbins,),
                                        numpy.int)
        self.bin_target_counts[...] = 4

def pcoord_loader_color_tracker(fieldname, coord_file, segment, single_point=False):
    """
    This function loads a 1-dimensional progress coordinate, performs some logic to track color,
    then returns the 2 dimensional progress coordinate to the system to be processed.
    In this tutorial, there are 2 dimensions specified in this file; runseg.sh returns one of them.
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
    # It's difficult to ascertain what is truly 'folded' and 'unfolded' for these without a prior
    # free energy profile; thankfully, we just need some rough estimates.  In the worst case scenario,
    # we devolve to the original Huber and Kim sampling scheme.
    color_bins = [(0.0,0.20),(11.60,float('inf'))]
    unknown_state = 2
    system = westpa.rc.get_system_driver()

    if single_point == True:
        npts = 1
    else:
        npts = system.pcoord_len

    coords = numpy.empty((npts), numpy.float32)
    colors = numpy.empty((npts), numpy.float32)
    #coords = numpy.empty((npts,system.pcoord_ndim), numpy.float32)
    #colors = numpy.empty((npts), numpy.float32)
    if single_point == True:
        colors[:] = unknown_state
        for istate,state_tuple in enumerate(color_bins):
            # Note that here, we are using the first dimension and first dimension alone.
            # The shape of the returned coord_raw is slightly different if single_point evalues
            # to true.
            # We evalulate whether or not we're in a state; if not, we leave it as in the
            # unknown state.
            # Swap this line to enable an N-dimensional pcoord, using the 1st dimension
            # as the state definition.
            #if coord_raw[0] >= state_tuple[0] and coord_raw[0] < state_tuple[1]:
            if coord_raw >= state_tuple[0] and coord_raw < state_tuple[1]:
                colors[:] = istate
        coords[:] = coord_raw[...]
    else:
        # If we're not the first point, we set the state to be what it was in the beginning
        # of the iteration.  We only want to update the state when we update a bin for purposes
        # of state tracking.
        # Swap lines to enable multiple pcoord dimensions, then change dimensions.
        #colors[:] = segment.pcoord[0][2]
        colors[:] = segment.pcoord[0][1]
        coords[:] = coord_raw[...]

    for istate,state_tuple in enumerate(color_bins):
        #if coords[-1,0] >= state_tuple[0] and coords[-1,0] < state_tuple[1]:
        if coords[-1] >= state_tuple[0] and coords[-1] < state_tuple[1]:
            colors[-1] = istate
    
    # We require different stacking behavior to return things in the proper order
    # depending on how many points we have.  I could probably clean this up.
    if single_point == True:
        # Again, swap lines.
        #segment.pcoord = numpy.hstack((coords[0,0],coords[0,1],colors[:]))
        segment.pcoord = numpy.hstack((coords[:],colors[:]))
    else:
        # This could easily be modified to return N dimensions.
        #segment.pcoord = numpy.swapaxes(numpy.vstack((coords[:,0],coords[:,1],colors[:])), 0, 1)
        segment.pcoord = numpy.swapaxes(numpy.vstack((coords[:],colors[:])), 0, 1)

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
