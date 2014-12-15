# Render aligned models for visualization.
#
# John D. Chodera <choderaj@mskcc.org> - 17 Feb 2013
#
# PREREQUISITES
#
# * PyMOL
# http://www.pymol.org/

# PARAMETERS

import sys
from ast import literal_eval

# Process only these targets, if specified.
# e.g. -targets '["SRC_HUMAN_PK0_P12931", "ABL1_HUMAN_PK0_P00519"]'
try:
    process_only_these_targets = literal_eval( sys.argv[ sys.argv.index('-targets') + 1 ] )
except ValueError:
    process_only_these_targets = False

if process_only_these_targets:
    print 'Processing only these targets:'
    print process_only_these_targets

width = 640
height = 480
dpi = -1
ray = 1

overwrite = False # if True, will overwrite existing files

#
# INITIALIZE PYMOL
#

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

import pymol.cmd

pymol.cmd.bg_color('white')
pymol.cmd.set('ray_trace_mode', '0')
pymol.cmd.set('ambient', '1')
pymol.cmd.set('reflect', '0')
pymol.cmd.set('antialias', '1')
pymol.cmd.set('two_sided_lighting', 'on')
pymol.cmd.set('cartoon_fancy_helices', '1')

#
# GET ABSOLUTE PATHS
#

import os.path

# Input files.
targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

#
# READ TEMPLATE AND TARGET INDICES
#

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
infile = open(targets_index_filename, 'r')
targets = [ line.strip() for line in infile ]
infile.close()
#print "targets:"
#print targets

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
infile = open(templates_index_filename, 'r')
templates = [ line.strip() for line in infile ]
infile.close()
#print "templates:"
#print templates

#
# SUBROUTINES
#

def make_contact_sheet(fnames,(ncols,nrows),(photow,photoh),
                       (marl,mart,marr,marb),
                       padding):
    """\
    Make a contact sheet from a group of filenames:

    fnames       A list of names of the image files
    
    ncols        Number of columns in the contact sheet
    nrows        Number of rows in the contact sheet
    photow       The width of the photo thumbs in pixels
    photoh       The height of the photo thumbs in pixels

    marl         The left margin in pixels
    mart         The top margin in pixels
    marr         The right margin in pixels
    marl         The left margin in pixels

    padding      The padding between images in pixels

    returns a PIL image object.
    """

    import Image
    import ImageFont
    import ImageDraw
    
    # Read in all images and resize appropriately
    #imgs = [Image.open(fn).resize((photow,photoh)) for fn in fnames]
    
    imgs = list()
    for image_filename in fnames:
        image = Image.open(image_filename)
        image = image.resize((photow,photoh), Image.ANTIALIAS)
        fontsize = 10
        #font = ImageFont.truetype("HelveticaNeueLight.ttf", fontsize)
        font = ImageFont.load_default()
        elements = os.path.split(image_filename)
        label = os.path.basename(elements[0])
        draw = ImageDraw.Draw(image)
        print label
        draw.text((fontsize, 0), label, (0,0,0), font=font)
        imgs.append(image)

    # Calculate the size of the output image, based on the
    #  photo thumb sizes, margins, and padding
    marw = marl + marr
    marh = mart + marb

    padw = (ncols-1)*padding
    padh = (nrows-1)*padding
    isize = (ncols*photow+marw+padw,nrows*photoh+marh+padh)

    # Create the new image. The background doesn't have to be white
    white = (255,255,255)
    inew = Image.new('RGB',isize,white)

    # Insert each thumb:
    for irow in range(nrows):
        for icol in range(ncols):
            left = marl + icol*(photow+padding)
            right = left + photow
            upper = mart + irow*(photoh+padding)
            lower = upper + photoh
            bbox = (left,upper,right,lower)
            try:
                img = imgs.pop(0)
            except:
                break
            inew.paste(img,bbox)
    return inew

# 
# RENDER IMAGES
#

original_directory = os.getcwd()

for target in targets:
    
    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    print "-------------------------------------------------------------------------"
    print "Imaging '%s'" % (target)
    print "-------------------------------------------------------------------------"

    # Process all templates.
    for template in templates:

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        model_filename = os.path.join(model_directory, 'aligned.pdb')
        if not os.path.exists(model_filename): continue

        image_filename = os.path.join(model_directory, 'image.png')
        if os.path.exists(image_filename) and not overwrite: continue
        
        print "-------------------------------------------------------------------------"
        print "Imaging '%s' => '%s'" % (template, target)
        print "-------------------------------------------------------------------------"
        
        os.chdir(model_directory)
        
        model_filename = 'aligned.pdb'
        image_filename = 'image.png'

        pymol.cmd.delete('all')
        pymol.cmd.load(model_filename)
        pymol.cmd.hide('lines', 'all')
        pymol.cmd.show('cartoon', 'all')
        pymol.cmd.color('red', 'all')
        pymol.cmd.zoom('all')

        core_filename = os.path.join(target_directory, 'core.txt')
        if os.path.exists(core_filename): 
            infile = open(core_filename, 'r')
            for line in infile:
                elements = line.split()
                segment_start = int(elements[0])
                segment_end = int(elements[1])
                pymol.cmd.color('grey', 'resi %d-%d' % (segment_start, segment_end))
            infile.close()

        pymol.cmd.png(image_filename, width, height, dpi, ray)

        import time
        time.sleep(3.0)

        os.chdir(original_directory)    

    #
    # GATHER SEQUENCE IDENTITIES
    #
        
    import deprecated_commands
    import os.path
    import numpy

    # Compile sequence identities.
    image_filenames = list()
    seqids = list()
    for template in templates:

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory): continue

        seqid_filename = os.path.join(model_directory, 'sequence-identity.txt')
        if not os.path.exists(seqid_filename): continue

        image_filename = os.path.join(model_directory, 'image.png')
        if not os.path.exists(image_filename): continue

        infile = open(seqid_filename, 'r')
        seqid = float(infile.readline().strip())
        infile.close()

        image_filenames.append(image_filename)
        seqids.append(seqid)

    
    indices = numpy.argsort(-numpy.array(seqids))
    image_filenames = [ image_filenames[i] for i in indices ]
    seqids = [ seqids[i] for i in indices ].reverse()
    print seqids

    # 
    # MAKE CONTACT SHEET
    #
    
    ncols = 10
    nrows = 10
    photow = width
    photoh = height
    marl = 50
    marr = marl 
    mart = marl
    marb = marl
    padding = marl

    print "Rendering contact sheet..."
    image = make_contact_sheet(image_filenames, (ncols,nrows), (photow,photoh), (marl,mart,marr,marb), padding)
    image_filename = os.path.join(target_directory, 'models.png')
    image.save(image_filename)

pymol.cmd.quit()

