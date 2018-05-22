#!/usr/bin/env dials.python

from __future__ import division

help_message = '''
ik heb ook geen idee
'''

# create phil_scope - Copied from other DIALS scripts...?
from libtbx.phil import parse
import numpy as np
from dials.array_family import flex
import iotbx.phil
import wx
from dials.algorithms.indexing import indexer

phil_scope = iotbx.phil.parse("""\
image_viewer {
  brightness = 120
    .type = int
  color_scheme = grayscale rainbow *heatmap invert
    .type = choice
  show_beam_center = True
    .type = bool
  show_resolution_rings = False
    .type = bool
  show_ice_rings = False
    .type = bool
  show_ctr_mass = False
    .type = bool
  show_max_pix = False
    .type = bool
  show_all_pix = True
    .type = bool
  show_shoebox = False
    .type = bool
  show_predictions = False
    .type = bool
  show_miller_indices = False
    .type = bool
  show_indexed = False
    .type = bool
  show_integrated = False
    .type = bool
  show_mask = False
    .type = bool
  show_mask2 = False
    .type = bool
  display = *image mean variance dispersion sigma_b \
            sigma_s threshold global_threshold
    .type = choice
  nsigma_b = 6
    .type = float(value_min=0)
  nsigma_s = 3
    .type = float(value_min=0)
  global_threshold = 0
    .type = float(value_min=0)
  kernel_size = 3,3
    .type = ints(size=2, value_min=1)
  min_local = 2
    .type = int
  gain = 1
    .type = float(value_min=0)
  sum_images = 1
    .type = int(value_min=1)
    .expert_level = 2
  untrusted_polygon = None
    .multiple = True
    .type = ints(value_min=0)
  d_min = None
    .type = float(value_min=0)
  mask = None
    .type = str
    .help = path to mask pickle file
}

predict_reflections = False
  .type = bool
  .help = Predict reflections if no reflections provided in input

include scope dials.algorithms.profile_model.factory.phil_scope
include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

""", process_includes=True)

class Script(object):

    def __init__(self):
        '''Initialise the script.'''
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = "usage: %s [options] experiment.json reflections.json" %libtbx.env.dispatcher_name

        # Create the parser and process options
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            read_experiments=True,
            read_reflections=True,
            read_datablocks=True)
        params, options = self.parser.parse_args(show_diff_phil=True)
        print params

    def run(self):
        from dials.util.options import flatten_experiments
        from dials.util.options import flatten_reflections
        from dials.util.options import flatten_datablocks
        from dials.util import log
        from logging import info, debug

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)
        experiments = flatten_experiments(params.input.experiments)
        reflections = flatten_reflections(params.input.reflections)
        datablocks = flatten_datablocks(params.input.datablock)
        if len(experiments) == 0:
            self.parser.print_help()
            return
        assert len(experiments) == 1

        imageset = experiments[0].imageset
        data = imageset.get_raw_data(0)[0]
#        detector = imageset.get_detector()
#        help(detector)
        for i in reflections[0].keys(): print i
        x_lim, y_lim = data.all()[1],data.all()[0]

        from dials.util.spotfinder_wrap import spot_wrapper
        self.wrapper = spot_wrapper(params=params)
        self.wrapper.display(
            datablock=datablocks[0],
            experiments=experiments,
            reflections=reflections)
        
        refl = indexer.indexer_base.map_spots_pixel_to_mm_rad(
            reflections[0],
            imageset.get_detector(), imageset.get_scan())

        indexer.indexer_base.map_centroids_to_reciprocal_space(
            refl, imageset.get_detector(), imageset.get_beam(),
            imageset.get_goniometer())
        reflections[0].extend(refl)

        for i in reflections[0].keys(): print i 

        for rlp in reflections[0]['rlp']:
            print rlp

        x, y = detector[reflection['panel']].millimeter_to_pixel(
            reflection['xyzcal.mm'][:2])
        x, y = map_coords(x+ 0.5, y + 0.5, reflection['panel'])

            
        quit()
        
        # Mark pixels using different reflection attributes:
        bbox = False
        calpx = True
        # Make px in bbox 0
        if bbox:
            for refls in reflections:
                for box in refls['bbox']:
                    box_pix = [(y,x) for x in range(box[2],box[3]) for y in range(box[0],box[1])]
                    for pix in box_pix:
                        if 0 < pix[0] < y_lim and 0 < pix[1] < x_lim:
                            data[pix] = 0

        # make px in call.px-10+10 0
        if calpx:
            for refls in reflections:
                for px in refls['xyzcal.px']:
                    box_px = [(y,x) for x in range(int(px[1]-10),int(px[1]+10)) for y in range(int(px[0]-10),int(px[0]+10))]
                    for pix in box_px:
                        if 0 < pix[0] < y_lim and 0 < pix[1] < x_lim:
                            data[pix] = 0

        tst_view(p=params,
                 refl=reflections,
                 datablock=None,
                 experiments=experiments)
        
    def median_filter_setup():
        # get images and experimental parameters
        imageset = experiments[0].imageset
        beam = imageset.get_beam()
        detector = imageset.get_detector()
        # Find image array (intensity values)
        data = imageset.get_raw_data(0)[0]
        # Limits / size of image
        lim_x, lim_y = data.all()[1],data.all()[0]
        # Box size for filtering
        box_size = 15

        # Generate separate jobs - this can be done smarter!!!
        print 'Creating easy_mp input...'
        from libtbx import easy_mp
        mp_list = []
        for x in range(box_size,lim_x-box_size):
            for y in range(box_size,lim_y-box_size):
                mp_list.append([x,y,box_size,data])

        # Run multiprocessor job
        print 'Starting MP-jobs...'
        result = easy_mp.pool_map(
            func=median_filt,
            args=mp_list)

        # Place results in image
        print 'Creating new image...'
        for i in result:
            data[i[1],i[0]]=i[2]

        # Saving image using 'merge' tool from dials.command_line
        print 'saving image...'
        from dials.command_line.merge_cbf import merge_cbf
        merge_cbf(imageset, n_images=1, out_prefix="tst_")
        
        # Show image in viewer using dials.image_viewer
        print 'showing image...'
        tst_view(p=params,
                 refl=reflections,
                 datablock=None,
                 experiments=experiments)

def median_filt(ar):
    '''
    Small method that is being started millions of times when median filtering the image
    '''
    x = ar[0]
    y = ar[1]
    box_size = ar[2]
    data = ar[3]
    xr = range(x-box_size, x+box_size)
    yr = range(y-box_size, y+box_size)
    i = int(np.median(np.array([data[yt,xt] for xt in xr for yt in yr])))
    return [x,y,i]

def tst_view(p=None,refl=None,datablock=None,experiments=None):
    '''
    Loads the dials.image_viewer with adapted image 
    '''
    import dials.command_line.image_viewer as im_viewer
    import libtbx.load_env
    runner = im_viewer.Script(
        params=p,
        reflections=refl,
        datablock=datablock,
        experiments=experiments)
    # Run the script
    runner()

def run(l):
    '''
    starter method, When halraiser is active error messages will be replaced by link to DIALS website
    '''
    finished = False # Determines type of error
    if finished:
        from dials.util import halraiser
        try:
            script = Script()
            script.run()
        except Exception as e:
            halraiser(e)
    else:
        script = Script()
        script.run()

def print_shit(par,h):
    '''
    help printer, prints type, and when h=True, the help page for an object.
    '''
    print type(par)
    if h: print help(par)
