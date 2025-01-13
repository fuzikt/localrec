#!/usr/bin/env python3

# **************************************************************************
# *
# * Authors:  Serban Ilca (serban@strubi.ox.ac.uk)
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *           J.M. de la Rosa Trevin
# *           Tibor Fuzik
# *
# * Oxford Particle Imaging Centre,
# * Division of Structural Biology, University of Oxford
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# **************************************************************************


import os
import sys
import textwrap

from distutils import spawn
import argparse

from lib.localrec import *
from lib.pyrelion import MetaData


class LocalizedReconstruction():

    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
                    Localized reconstruction of subparticles. Normally the script is run in four steps:
                        1. Prepare particles.
                        2. Create subparticles.
                        3. Extract subparticles.
                        4. Reconstruct subparticles.
            '''))

        # Parameter groups
        general = self.parser.add_argument_group('General parameters')
        steps = self.parser.add_argument_group('Steps', 'Several steps can be combined in one run.')
        prepareParticles = self.parser.add_argument_group('Prepare particles')
        createSubparticles = self.parser.add_argument_group('Create subparticles')
        extractSubparticles = self.parser.add_argument_group('Extract subparticles')
        reconstructSubparticles = self.parser.add_argument_group('Reconstruct subparticles')
        expertOptions = self.parser.add_argument_group('Expert options')

        # Shortcuts
        add = general.add_argument
        adds = steps.add_argument
        addpp = prepareParticles.add_argument
        addcs = createSubparticles.add_argument
        addes = extractSubparticles.add_argument
        addrs = reconstructSubparticles.add_argument
        addeo = expertOptions.add_argument

        # General parameters
        add('input_star', nargs='?', help="Input STAR filename with particles.")
        add('--output', default='subparticles',
              help="Output root for results.")

        # Parameters for "Steps" group
        adds('--prepare_particles', action='store_true',
              help="Prepare particles for extracting subparticles.")
        adds('--create_subparticles', action='store_true',
              help="Calculate the cooridnates and Euler angles for the subparticles.")
        adds('--extract_subparticles', action='store_true',
              help="Extract subparticles from particle images.")
        adds('--reconstruct_subparticles', action='store_true',
              help="Calculate a reconstruction of the subunit from the subparticles. "
                 "Subparticle coordinates are read from "
                 "[output]_subparticles.star and [output]_subparticles_subtracted.star")

        # Parameters for "Prepare particles" group
        addpp('--masked_map',
              help="Create another set of particles with partial signal subtraction using this map.")
        addpp('--extract_from_micrographs', action='store_true',
              help="Extract subparticles directly from micrographs instead from particle images. Cannot be used with signal subtraction enabled.")

        # Parameters for "Create subparticles" group
        addcs('--angpix', type=float, default=1, help="Pixel size (A). (default: 1)")
        addcs('--sym', default="C1", help="Symmetry of the particle. (default: C1)")
        addcs('--randomize', action='store_true',
              help="Randomize the order of the symmetry matrices. "
                   "Useful for preventing preferred orientations.")
        addcs('--relax_symmetry', action='store_true',
              help="Create one random subparticle for each particle ")
        addcs('--vector', default="0,0,1", help="Vector defining the location of the subparticle. (default: 0,0,1)")
        addcs('--align_subparticles', action='store_true',
              help="Align subparticles to the standard orientation.")
        addcs('--length',
              help="Alternative length of the vector. Use to adjust the "
                 "subparticle center (A). (default: length of the given vector)")
        addcs('--cmm',
              help="A CMM file defining the location(s) of the subparticle(s) "
                 "(use instead of --vector). Coordinates should be in Angstrom.")
        addcs('--unique', type=float, default=-1,
              help="Keep only unique subparticles within angular distance "
                 "(useful to remove overlapping subparticles on symmetry axis).")
        addcs('--mindist', type=float, default=-1,
              help="Minimum distance between the subparticles in the image "
                 "(all overlapping ones will be discarded; pixels).")
        addcs('--side', type=float, default=-1,
              help="Keep only particles within specified angular distance from "
                 "side views (all others will be discarded; degrees).")
        addcs('--top', type=float, default=-1,
              help="Keep only particles within specified angular distance from "
                 "top views (all others will be discarded; degrees).")

        # Parameters for "Extract subparticles" group
        addes('--subparticle_size', type=int,
            help="Size of the subparticle box (pixels).")
        addes('--rescale_size', type=int, default=-1,
              help="Rescale (bin) subparticle to this box (pixels). (Default: -1 => no rescale)")
        addes('--np', type=int, default=1, help="Number of MPI procs. (default: 1)")
        addes('--only_extract_unfinished', action='store_true',
              help="Extract only particles not extracted in previous run.")
        addes('--invert_contrast', action='store_true',
              help="Use this option when extracting from micrographs that were not inverted (i.e. black particles on white background).")
        addes('--normalize', action='store_true',
              help="Normalize the extracted particles. Useful when extracting from micrographs.")

        # Parameters for "Reconstruct subparticles" group
        addrs('--do_halves', action='store_true',
              help="Make reconstruction of the random halves.")
        addrs('--j', type=int, default=1, help="Number of threads.")
        addrs('--maxres', type=float, help="Maximum resolution of the reconstruction (A). (default: Nyquist)")
        addrs('--subsym', default="C1", help="Symmetry of the subparticle. (default: C1)")

        # Expert parameters
        addeo('--library_path', type=str, default='',
              help="define LD_LIBRARY_PATH used. Default: empty")


    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("\nError: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if len(sys.argv) == 1:
            self.usage()
        # Check that required software is in PATH

        if not (spawn.find_executable("relion_refine")):
            self.error("Relion not found.",
                       "Make sure Relion programs are in $PATH.")

        # Check that required parameters are given for each mode
        if args.prepare_particles or args.create_subparticles or args.extract_subparticles:
            if len(sys.argv) == 1:
                self.error("No input particles STAR file given.")
            if not os.path.exists(args.input_star):
                self.error("\nInput file '%s' not found." % args.input_star)

        if args.masked_map:
            if args.extract_from_micrographs:
               self.error("You cannot use signal subtraction and extracting from micrographs option together. Please, choose one of them.")

        if args.extract_subparticles:
            if not args.subparticle_size:
                self.error("Parameter --subparticle_size is required.")

        if args.reconstruct_subparticles:
            if not args.output:
                self.error("Parameter --output not specified. Cannot find STAR file with subparticles.")

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        # Validate input arguments and required software
        self.validate(args)

        try: 
            os.makedirs(args.output)
        except OSError:
            if not os.path.isdir(args.output):
                raise


        md_in = MetaData(args.input_star)

        if md_in.version == "3.1":
            if args.extract_from_micrographs:
                mrcFilename = md_in.data_particles[0].rlnMicrographName
                if hasattr(md_in.data_optics[0], 'rlnMicrographOriginalPixelSize'):
                    apix = md_in.data_optics[0].rlnMicrographOriginalPixelSize
                    print("Using micrograph pixel size from star file: %s A/px" % apix)
                else:
                    apix = args.angpix
                    print("Micrograph pixel size not found in star file, using user defined value: %s A/px" % apix)
            else:
                mrcFilename = md_in.data_particles[0].rlnImageName.split("@")[1]
                apix = md_in.data_optics[0].rlnImagePixelSize
                print("Using particle pixel size from star file: %s A/px" % apix)
        else:
            apix = args.angpix
            print("Using user defined pixel size: %s A/px" % apix)
            if args.extract_from_micrographs:
                mrcFilename = md_in.data_[0].rlnMicrographName
            else:
                mrcFilename = md_in.data_[0].rlnImageName.split("@")[1]


        mrcFile = open(mrcFilename, "rb")
        particleImageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
        particleImageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcFile.close()
        if args.extract_from_micrographs:
            print("Micrograph image size:")
        else:
            print("Particle image size:")
        print("%s, %s pixels\n" % (particleImageSizeX, particleImageSizeY))

        if args.prepare_particles:
            print("Preparing particles for extracting subparticles.")
            create_initial_stacks(args.input_star, apix, args.masked_map, args.output, args.extract_from_micrographs, args.library_path)
            print("\nFinished preparing the particles!\n")

        if args.create_subparticles:
            # Load subparticle vectors either from Chimera CMM file or from
            # command line (command and semi-colon separated)
            # Distances can also be specified to modify vector lengths
            subparticle_vector_list = load_vectors(args.cmm, args.vector,
                                                   args.length, apix)

            particles_star = args.output + "/particles.star"

            if not os.path.exists(args.output + "/particles.star"):
                self.error("Input file '%s not found. "
                           "Run the script first with --prepare_particles option."
                           % particles_star)

            md = MetaData(particles_star)

            print("Creating subparticles...")

            # Generate symmetry matrices with Relion convention
            symmetry_matrices = matrix_from_symmetry(args.sym, args.library_path)

            # Initialize progress bar
            progressbar = ProgressBar(width=60, total=len(md))

            # Define some conditions to filter subparticles
            filters = load_filters(radians(args.side), radians(args.top), args.mindist)

            # Compute all subparticles (included subtracted if masked_map given)
            mdOut = MetaData()
            mdOutSub = MetaData()

            if md.version == "3.1":

                #if extracted from micrographs set the image apix same as micrograph
                if args.extract_from_micrographs:
                    for optic_group_nr in range(md.size('data_optics')):
                        md.data_optics[optic_group_nr].rlnImagePixelSize = apix

                mdOut.version = "3.1"
                mdOut.addDataTable("data_optics", True)
                mdOut.addLabels("data_optics", md.getLabels("data_optics"))
                mdOut.addData("data_optics", getattr(md, "data_optics"))
                mdOut.addDataTable("data_particles", True)

                mdOutSub.version = "3.1"
                mdOutSub.addDataTable("data_optics", True)
                mdOutSub.addLabels("data_optics", md.getLabels("data_optics"))
                mdOutSub.addData("data_optics", getattr(md, "data_optics"))
                mdOutSub.addDataTable("data_particles", True)

                particleTableName = "data_particles"
            else:
                mdOut.addDataTable("data_", True)
                mdOutSub.addDataTable("data_", True)
                particleTableName = "data_"

            for particle in md:
                subparticles, subtracted = create_subparticles(particle,
                                                       symmetry_matrices,
                                                       subparticle_vector_list,
                                                       particleImageSizeX,
                                                       particleImageSizeY,
                                                       args.randomize,
                                                       args.output,
                                                       args.unique,
                                                       len(mdOut),
                                                       args.align_subparticles,
                                                       args.masked_map,
                                                       args.extract_from_micrographs,
                                                       True,
                                                       filters,
                                                       apix)

                mdOut.addData(particleTableName, subparticles)
                mdOutSub.addData(particleTableName, subtracted)

                progressbar.notify()

            md.removeLabels(particleTableName, 'rlnOriginZ', 'rlnOriginalName')

            write_output_starfiles(md.getLabels(particleTableName), mdOut, mdOutSub, args.output)


            #if args.extract_from_micrographs:
            unique_micrographs(md, args.output, 'micrographs')

            if len(mdOutSub):
                unique_micrographs(mdOutSub, args.output, 'micrographs_subtracted')

            print("\nFinished creating the subparticles!\n")

        if args.extract_subparticles:
            print("Extracting subparticles...")
            if args.extract_from_micrographs:
                if not args.create_subparticles:
                    md = MetaData(args.output + "/particles.star")                    
                extract_subparticles(args.subparticle_size, args.rescale_size, args.np, args.masked_map, args.output, args.library_path, args.only_extract_unfinished, args.invert_contrast, args.normalize, False, getattr(md,particleTableName)[-2].rlnMicrographName.split('/').pop(0))
            else:
                extract_subparticles(args.subparticle_size, args.rescale_size, args.np, args.masked_map, args.output, args.library_path, args.only_extract_unfinished, args.invert_contrast, args.normalize, True, args.output)
            print("\nFinished extracting the subparticles!\n")

        if args.reconstruct_subparticles:
            print("Reconstructing subparticles...")
            reconstruct_subparticles(args.j, args.output, args.maxres, args.subsym, apix, args.do_halves, args.library_path)
            print("\nFinished reconstructing the subparticles!\n")

        print("\nAll done!\n")

if __name__ == "__main__":
    LocalizedReconstruction().main()
