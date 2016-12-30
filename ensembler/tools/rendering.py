import subprocess
import tempfile
import os
import warnings
import shutil
import ensembler
import mdtraj
from PIL import Image, ImageFont, ImageDraw


class PymolRender(object):
    def __init__(self, width=640, height=480):
        self.width = width
        self.height = height

    def run_pymol(self, input_filepath=None, img_filepath=None):
        input_filepath_abs = os.path.abspath(input_filepath)
        img_filepath_abs = os.path.abspath(img_filepath)
        pml_file_text = """\
load {INPUT_FILEPATH}

bg white
set ray_trace_mode, 0
set ambient, 1
set reflect, 0
set antialias, 1
set two_sided_lighting, on
set cartoon_fancy_helices, 1

set_color dblue, [0.204, 0.298, 0.384]

hide lines, all
show cartoon, all
color grey20, ss h
color grey20, ss s
color grey, ss l+''
zoom all

ray 640,480
png {IMG_FILEPATH}
"""
        pml_file_text = pml_file_text.replace('{INPUT_FILEPATH}', input_filepath_abs)
        pml_file_text = pml_file_text.replace('{IMG_FILEPATH}', img_filepath_abs)

        tmpdir = tempfile.mkdtemp()

        try:
            pml_filepath = os.path.join(tmpdir, 'render.pml')
            with open(pml_filepath, 'w') as pml_file:
                pml_file.write(pml_file_text)

            subprocess.call(['pymol', '-qc', pml_filepath])

        finally:
            shutil.rmtree(tmpdir)


class Rendering(object):
    def __init__(self):
        self.targets, templates_resolved_seq = ensembler.core.get_targets_and_templates()
        self.templates = templates_resolved_seq

    def render_template_structure_resolved(self, templateid):
        input_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, templateid+'.pdb')
        img_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, templateid+'.png')
        pr = PymolRender()
        pr.run_pymol(input_filepath=input_filepath, img_filepath=img_filepath)

    def render_template_structure_modeled_loops(self, templateid):
        input_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, templateid+'.pdb')
        img_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, templateid+'.png')
        pr = PymolRender()
        pr.run_pymol(input_filepath=input_filepath, img_filepath=img_filepath)

    def render_model(self, targetid, templateid):
        input_filepath = os.path.join(ensembler.core.default_project_dirnames.models, targetid, templateid, 'model.pdb.gz')
        img_filepath = os.path.join(ensembler.core.default_project_dirnames.models, targetid, templateid, 'model.png')
        pr = PymolRender()
        pr.run_pymol(input_filepath=input_filepath, img_filepath=img_filepath)


class RenderOnGrid(object):
    def __init__(self, nrows=1, ncols=2, width=640, height=480):
        self.structure_filepaths = []
        self.labels = []
        self.nrows = nrows
        self.ncols = ncols
        self.width = width
        self.height = height
        self._tmpdir = None

    def add_template_structures_resolved(self, templateids, labels=None):
        if type(templateids) is str:
            templateids = [templateids]

        if labels is None:
            labels = [None] * len(templateids)

        for i in range(len(templateids)):
            templateid = templateids[i]
            template_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_resolved, templateid+'.pdb')
            self.structure_filepaths.append(template_filepath)
            self.labels.append(labels[i])

    def add_template_structures_modeled_loops(self, templateids, labels=None):
        if type(templateids) is str:
            templateids = [templateids]

        if labels is None:
            labels = [None] * len(templateids)

        for i in range(len(templateids)):
            templateid = templateids[i]
            template_filepath = os.path.join(ensembler.core.default_project_dirnames.templates_structures_modeled_loops, templateid+'.pdb')
            self.structure_filepaths.append(template_filepath)
            self.labels.append(labels[i])

    def add_model(self, targetid, templateid, label=None):
        structure_filepath = os.path.join(ensembler.core.default_project_dirnames.models, targetid, templateid, 'model.pdb.gz')
        self.structure_filepaths.append(structure_filepath)
        self.labels.append(label)

    def render(self, tiled_img_filename):
        if len(self.structure_filepaths) == 0:
            warnings.warn('No structures selected.')
            return
        self._tmpdir = tempfile.mkdtemp()
        try:
            self._align_structures()
            self._render_aligned_structures()
            self._tile_images(tiled_img_filename)
        finally:
            shutil.rmtree(self._tmpdir)

    def _align_structures(self):
        for irow in range(self.nrows):
            row_first_col = irow * self.ncols
            self._ref_traj = mdtraj.load_pdb(self.structure_filepaths[row_first_col])
            self._heavy_atoms = self._ref_traj.topology.select('not name H')
            for icol in range(self.ncols):
                alignment_index = row_first_col + icol
                structure_filepath = self.structure_filepaths[alignment_index]
                self._align_structure(structure_filepath, alignment_index)

    def _align_structure(self, structure_filepath, alignment_index):
        traj = mdtraj.load_pdb(structure_filepath)
        traj.superpose(self._ref_traj, atom_indices=self._heavy_atoms, parallel=False)
        aligned_structure_filepath = os.path.join(self._tmpdir, 'aligned%d.pdb' % alignment_index)
        traj.save(aligned_structure_filepath)

    def _render_aligned_structures(self):
        pr = PymolRender()
        for structure_index in range(len(self.structure_filepaths)):
            aligned_structure_filepath = os.path.join(self._tmpdir, 'aligned%d.pdb' % structure_index)
            # if structure_index == 0:
            #     shutil.copy(aligned_structure_filepath, '.')
            img_filepath = os.path.join(self._tmpdir, 'img%d.png' % structure_index)   # TODO
            pr.run_pymol(input_filepath=aligned_structure_filepath, img_filepath=img_filepath)

    def _tile_images(self, tiled_img_filename, marl=0, marr=0, mart=0, marb=0, padding=0):
        imgs = []
        for structure_index in range(len(self.structure_filepaths)):
            image_filename = os.path.join(self._tmpdir, 'img%d.png' % structure_index)
            image = Image.open(image_filename)
            # image = image.resize((self.width, self.height), Image.ANTIALIAS)
            label = self.labels[structure_index]
            if label is not None:
                fontsize = 10
                # font = ImageFont.truetype("HelveticaNeueLight.ttf", fontsize)
                font = ImageFont.load_default()
                draw = ImageDraw.Draw(image)
                print(label)
                draw.text((fontsize, 0), label, (0,0,0), font=font)
            imgs.append(image)

        # Calculate the size of the output image, based on the
        #  photo thumb sizes, margins, and padding
        marw = marl + marr
        marh = mart + marb

        padw = (self.ncols-1)*padding
        padh = (self.nrows-1)*padding
        isize = (self.ncols*self.width+marw+padw, self.nrows*self.height+marh+padh)

        # Create the new image. The background doesn't have to be white
        bgcolor = (255, 255, 255)
        image = Image.new('RGB', isize, bgcolor)

        # Insert each thumb:
        for irow in range(self.nrows):
            for icol in range(self.ncols):
                left = marl + icol*(self.width+padding)
                right = left + self.width
                upper = mart + irow*(self.height+padding)
                lower = upper + self.height
                bbox = (left, upper, right, lower)
                try:
                    img = imgs.pop(0)
                except:
                    break
                image.paste(img, bbox)

        image.save(tiled_img_filename)
