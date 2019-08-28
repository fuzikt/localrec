# **************************************************************************
# *
# * Authors:  Serban Ilca
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *           J.M. de la Rosa Trevin
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

import math
import random
import shutil
import struct
from itertools import izip
from glob import glob
from matrix3 import *
from vector3 import *
from euler import *
from os.path import splitext
from os.path import basename
from os.path import exists
from os.path import isdir
from os.path import islink
from pyrelion import MetaData


def cleanPath(*paths):
    """ Remove a list of paths, either folders or files"""
    for p in paths:
        if exists(p):
            if isdir(p):
                if islink(p):
                    os.remove(p)
                else:
                    shutil.rmtree(p)
            else:
                os.remove(p)


def cleanPattern(pattern):
    """ Remove all files that match the pattern. """
    files = glob(pattern)
    cleanPath(*files)


def within_mindist(p1, p2, mindist):
    """ Returns True if two particles are closer to each other
    than the given distance in the projection. """

    x1 = p1.rlnCoordinateX
    y1 = p1.rlnCoordinateY
    x2 = p2.rlnCoordinateX
    y2 = p2.rlnCoordinateY
    distance_sqr = (x1 - x2) ** 2 + (y1 - y2) ** 2
    mindist_sqr = mindist ** 2

    return distance_sqr < mindist_sqr


def within_unique(p1, p2, unique):
    """ Returns True if two particles are closer to each other
    than the given angular distance. """

    v1 = vector_from_two_eulers(p1.rlnAnglePsi, p1.rlnAngleTilt)
    v2 = vector_from_two_eulers(p2.rlnAnglePsi, p2.rlnAngleTilt)

    dp = dot_product(v1, v2) / (v1.length() * v2.length())

    if dp < -1:
        dp = -1.000

    if dp > 1:
        dp = 1.000

    angle = math.acos(dp)

    return angle <= math.radians(unique)


def filter_unique(subparticles, subpart, unique):
    """ Return True if subpart is not close to any other subparticle
        by unique (angular distance).
        For this function we assume that subpart is not contained
        inside."""
    for sp in subparticles:
        if within_unique(sp, subpart, unique):
            return False

    return True


def filter_mindist(subparticles, subpart, mindist):
    """ Return True if subpart is not close to any other subparticle
    by mindist. """
    for sp in subparticles:
        if (sp.rlnImageName[:6] != subpart.rlnImageName[:6] and
                within_mindist(sp, subpart, mindist)):
            return False

    return True


def filter_side(subpart, side):
    return (abs(abs(subpart.rlnAngleTilt) - radians(90)) < side)


def filter_top(subpart, top):
    return (abs(abs(subpart.rlnAngleTilt) - radians(180)) < top)


def filter_subparticles(subparticles, filters):
    return [sp for sp in subparticles
            if all(f(subparticles, sp) for f in filters)]


def angles_to_radians(particle):
    """ Convert the particle angles to radians. """
    particle.rlnAnglePsi = math.radians(particle.rlnAnglePsi)
    particle.rlnAngleTilt = math.radians(particle.rlnAngleTilt)
    particle.rlnAngleRot = math.radians(particle.rlnAngleRot)


def angles_to_degrees(particle):
    """ Convert the particle angles to radians. """
    particle.rlnAnglePsi = math.degrees(particle.rlnAnglePsi)
    particle.rlnAngleTilt = math.degrees(particle.rlnAngleTilt)
    particle.rlnAngleRot = math.degrees(particle.rlnAngleRot)


def create_subparticles(particle, symmetry_matrices, subparticle_vector_list,
                        part_image_sizeX, part_image_sizeY, randomize, output,
                        unique, subparticles_total, align_subparticles,
                        subtract_masked_map, extract_from_micrographs, do_create_star, filters, ang_pix):
    """ Obtain all subparticles from a given particle and set
    the properties of each such subparticle. """

    part_filename = splitext(basename(particle.rlnImageName))[0]

    # We convert the particle angles to radian for further computations
    angles_to_radians(particle)

    # Euler angles that take particle to the orientation of the model

    rot = particle.rlnAngleRot
    tilt = particle.rlnAngleTilt
    psi = particle.rlnAnglePsi

    matrix_particle = matrix_from_euler(rot, tilt, psi)

    subparticles = []
    subtracted = []
    subpart_id = 1
    subparticles_total += 1

    symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

    if randomize:
        # randomize the order of symmetry matrices, prevents preferred views
        random.shuffle(symmetry_matrix_ids)

    for subparticle_vector in subparticle_vector_list:
        matrix_from_subparticle_vector = subparticle_vector.matrix()

        for symmetry_matrix_id in symmetry_matrix_ids:
            # symmetry_matrix_id can be later written out to find out
            # which symmetry matrix created this subparticle
            symmetry_matrix = symmetry_matrices[symmetry_matrix_id - 1]

            subpart = particle.clone()

            m = matrix_multiply(matrix_particle, (matrix_multiply(matrix_transpose(symmetry_matrix), matrix_transpose(matrix_from_subparticle_vector))))

            if align_subparticles:
                rotNew, tiltNew, psiNew = euler_from_matrix(m)
            else:
                m2 = matrix_multiply(matrix_particle, matrix_transpose(symmetry_matrix))
                rotNew, tiltNew, psiNew = euler_from_matrix(m2)

            # save Euler angles that take the model to the orientation of the subparticle

            subpart.rlnAngleRot = rotNew
            subpart.rlnAngleTilt = tiltNew
            subpart.rlnAnglePsi = psiNew

            # subparticle origin
            d = subparticle_vector.distance()
            x = -m.m[0][2] * d + particle.rlnOriginX
            y = -m.m[1][2] * d + particle.rlnOriginY
            z = -m.m[2][2] * d

            # modify the subparticle defocus paramaters by its z location
            if hasattr(particle, 'rlnDefocusU'):
                z_ang = ang_pix * z
                subpart.rlnDefocusU = particle.rlnDefocusU + z_ang
                subpart.rlnDefocusV = particle.rlnDefocusV + z_ang

            # save the subparticle coordinates (integer part) relative to the
            # user given image size and as a small shift in the origin (decimal part)
            x_d, x_i = math.modf(x)
            y_d, y_i = math.modf(y)

            if not extract_from_micrographs:
                subpart.rlnCoordinateX = int(part_image_sizeX / 2) - x_i
                subpart.rlnCoordinateY = int(part_image_sizeY / 2) - y_i
            else:
                subpart.rlnCoordinateX = particle.rlnCoordinateX - x_i
                subpart.rlnCoordinateY = particle.rlnCoordinateY - y_i

            subpart.rlnOriginX = x_d
            subpart.rlnOriginY = y_d

            overlaps = (unique >= 0 and
                        not filter_unique(subparticles, subpart, unique))

            if not overlaps:
                subpart.rlnImageName = "%06d@%s/%s.mrcs" % (subpart_id, output, part_filename)
                subpart.rlnParticleName = str(subparticles_total)

                if not extract_from_micrographs:
                    subpart.rlnMicrographName = "%s/%s.mrc" % (output, part_filename)
                    subparticles.append(subpart)
                    subpart_id += 1
                    subparticles_total += 1
                else:
                    subpart.rlnMicrographName = particle.rlnMicrographName
                    if not ((subpart.rlnCoordinateX<0) or (subpart.rlnCoordinateY<0) or (subpart.rlnCoordinateX>part_image_sizeX) or (subpart.rlnCoordinateY>part_image_sizeY)):
                        subparticles.append(subpart)
                        subpart_id += 1
                        subparticles_total += 1


    if subtract_masked_map:
        subtracted = clone_subtracted_subparticles(subparticles, output)

    if filters:
        subparticles = filter_subparticles(subparticles, filters)

        if subtract_masked_map:
            subtracted = clone_subtracted_subparticles(subparticles, output)

    return subparticles, subtracted


def create_symmetry_related_particles(particle, symmetry_matrices,
                                      keep_one=False):
    """ Return all related particles from the given symmetry matrices.
    If keep_one is True, randomly select only one of these equivalent
    particles. NOTE: Input particle should already contains angles in radians. """
    new_particles = []

    rot = particle.rlnAngleRot
    tilt = particle.rlnAngleTilt
    psi = particle.rlnAnglePsi

    matrix_particle = matrix_from_euler(rot, tilt, psi)

    for symmetry_matrix in symmetry_matrices:
        m = matrix_multiply(matrix_particle, matrix_transpose(symmetry_matrix))
        rotNew, tiltNew, psiNew = euler_from_matrix(m)

        new_particle = particle.clone()

        new_particle.rlnAngleRot = rotNew
        new_particle.rlnAngleTilt = tiltNew
        new_particle.rlnAnglePsi = psiNew
        angles_to_degrees(new_particle)
        new_particles.append(new_particle)

    if keep_one:
        new_particles = random.sample(new_particles, 1)

    return new_particles


def clone_subtracted_subparticles(subparticles, output):
    subparticles_subtracted = []

    for sp in subparticles:
        sp_new = sp.clone()
        sp_new.rlnImageName = add_suffix(sp.rlnImageName)
        sp_new.rlnMicrographName = add_suffix(sp.rlnMicrographName)
        subparticles_subtracted.append(sp_new)

    return subparticles_subtracted


def add_suffix(filename, output='particles'):
    return filename.replace('%s_' % output,
                           '%s_subtracted_' % output)


def create_star(subparticles, star_filename):
    """ Create a Relion style STAR file for extracting (using relion_preprocess)
    all the subparticles for a given particle. """

    md = MetaData()
    md.addLabels("rlnCoordinateX", "rlnCoordinateY")
    md.addData(subparticles)
    md.write(star_filename)


def load_vectors(cmm_file, vectors_str, distances_str, angpix):
    """ Load subparticle vectors either from Chimera CMM file or from
    a vectors string. Distances can also be specified for each vector
    in the distances_str. """

    if cmm_file:
        subparticle_vector_list = vectors_from_cmm(cmm_file, angpix)
    else:
        subparticle_vector_list = vectors_from_string(vectors_str)

    if distances_str:
        # Change distances from A to pixel units
        subparticle_distances = [float(x) / angpix for x in
                                 distances_str.split(',')]

        if len(subparticle_distances) != len(subparticle_vector_list):
            raise Exception("Error: The number of distances does not match "
                            "the number of vectors!")

        for vector, distance in izip(subparticle_vector_list,
                                     subparticle_distances):
            if distance > 0:
                vector.set_distance(distance)
            else:
                vector.compute_distance()
    else:
        for vector in subparticle_vector_list:
            vector.compute_distance()

    print("Using vectors:")

    for subparticle_vector in subparticle_vector_list:
        print("Vector: ")
        subparticle_vector.normalize()
        subparticle_vector.compute_matrix()
        subparticle_vector.print_vector()
        print("")
        print("Length: %.2f pixels" % subparticle_vector.distance())
    print("")

    return subparticle_vector_list


def load_filters(side, top, mindist):
    """ Create some filters depending on the conditions imposed by the user.
    Each filter will return True if the subparticle will be kept in the
    subparticles list.
    """
    filters = []

    if side > 0:
        filters.append(lambda x, y: filter_side(y, side))

    if top > 0:
        filters.append(lambda x, y: filter_top(y, top))

    if mindist > 0:
        filters.append(lambda x, y: filter_mindist(x, y, mindist))

    return filters

def splitMrcStack(stackFile, outFile):
    # extract filename, image index
    imageIndex, mrcsFilename = stackFile.split("@")
    imageIndex = int(imageIndex)

    # read in mrcs file, open mrc file for output
    mrcsFile = open(mrcsFilename, "rb")
    mrcFile = open(outFile, 'wb+')

    # get image size
    imageSize = int(struct.unpack('i', mrcsFile.read(4))[0])

    # write header
    mrcsFile.seek(0, 0)
    chunkSize = 1024
    mrcHeader = mrcsFile.read(chunkSize)
    mrcFile.write(mrcHeader)

    # change Z dimension to 1
    mrcFile.seek(8, 0)
    mrcFile.write(b"\x01")
    mrcFile.seek(1024, 0)

    # write mrc data file
    mrcsFile.seek(imageSize ** 2 * 4 * (imageIndex - 1) + 1024, 0)
    chunkSize = imageSize ** 2 * 4
    mrcImage = mrcsFile.read(chunkSize)
    mrcFile.write(mrcImage)

    mrcsFile.close()
    mrcFile.close()

def split_particle_stacks(extract_from_micrographs, inputStar, inputStack, output, filename_prefix, deleteStack):
    """ Read a STAR file with particles and write as individual images.
    If a stack of images is given, use these instead of the images from the STAR file.
    Also write a new STAR file pointing to these images.
    This function requires that the script is run within Scipion Python environment. """

    md = MetaData(inputStar)
    md.addLabels('rlnOriginalName')

    # Initialize progress bar
    progressbar = ProgressBar(width=60, total=len(md))

    for i, particle in enumerate(md, start=1):
        outputImageName = '%s/%s_%06d.mrc' % (output, filename_prefix, i)

        # split stacks only if not extracting directly from micrographs
        if not extract_from_micrographs:
            if inputStack:
                splitMrcStack('%06d@%s' %(i, inputStack), outputImageName)
                particle.rlnOriginalName = '%s/%06d@%s' %(output, i, inputStack)
            else:
                splitMrcStack(particle.rlnImageName, outputImageName)
                particle.rlnOriginalName = particle.rlnImageName
            particle.rlnImageName = outputImageName
            particle.rlnMicrographName = outputImageName
        else:
            particle.rlnOriginalName = particle.rlnImageName
            particle.rlnImageName = '%s/%s' % (output, particle.rlnMicrographName.split('/').pop())

        progressbar.notify()

    print("\n")
    md.write("%s/%s.star" % (output, filename_prefix))

    if inputStack and deleteStack:
        cleanPath(inputStack)



def create_initial_stacks(input_star, angpix, masked_map, output, extract_from_micrographs, library_path):
    """ Create initial particle stacks (and STAR files) to be used
    for extraction of subparticles. If the masked_map is passed,
    another stack with subtracted particles will be created. """

    print(" Creating particle images from which nothing is subtracted...")
    split_particle_stacks(extract_from_micrographs, input_star, None, output, 'particles', deleteStack=False)

    if masked_map:
        print(" Creating particle images from which the projections of the masked "
              "particle reconstruction have been subtracted...")

        subtractedStackRoot = "%s/particles_subtracted" % output

        md = MetaData(input_star)
        args = " --i %s --o %s --ang %s --subtract_exp --angpix %s "
        if "rlnDefocusU" in md.getLabels():
            args = args + "--ctf "
        else:
            print ("\nWarning: no CTF info found in %s!\n"
                   "The subtraction will be performed without CTF correction.\n" % input_star)
        run_command("relion_project" + args %
                    (masked_map, subtractedStackRoot, input_star, angpix),"",library_path)
        run_command("mv %s.star %s_orig.star" % (subtractedStackRoot, subtractedStackRoot),"",library_path)

        subtractedStack = subtractedStackRoot + '.mrcs'

        split_particle_stacks(extract_from_micrographs, input_star, subtractedStack, output, 'particles_subtracted', deleteStack=True)


def extract_subparticles(subpart_size, np, masked_map, output, library_path, only_extract_unfinished, invert_contrast, normalize, deleteParticles, outDir):
    """ Extract subparticles images from each particle
    (Using 'relion_preprocess' as if the particle was a micrograph.
    Notice that this command line works only in Relion 1.4, not 2.0"""

    additional_parameters = ''

    if np == 1:
        cmd = 'relion_preprocess '
    else:
        cmd = 'mpirun -np %s relion_preprocess_mpi ' % np
        additional_parameters='--max_mpi_nodes %s' % np

    if only_extract_unfinished:
        additional_parameters = additional_parameters+' --only_extract_unfinished'
    if invert_contrast:
        additional_parameters = additional_parameters+' --invert_contrast'
    if normalize:
        bg_radius = int(0.75*subpart_size/2)
        additional_parameters = additional_parameters+' --norm --bg_radius '+str(bg_radius)

    def run_extract(suffix=''):
        args = ('--i "%s/micrographs%s.star" --part_star subparticles_tmp.star --extract --extract_size %s --reextract_data_star '
                '"%s%s.star" %s') % (output, suffix, subpart_size, output, suffix, additional_parameters)
        run_command(cmd + args, "", library_path)

        print(" Cleaning up temporary files...")

        f = open('subparticles_tmp.star','r')
        filedata = f.read()
        f.close()
        f = open('%s%s.star' % (output,suffix),'w')
        f.write(filedata.replace('Particles/%s' % outDir, output))
        f.close()

        try:
            os.remove('%s.star' % 'subparticles_tmp')
        except OSError:
            pass

        if deleteParticles:
            cleanPattern('%s/particles%s_??????.mrc' % (output, suffix))

    run_extract()  # Run extraction without subtracted density

    if masked_map:
        run_extract('_subtracted')

    print(" Moving subparticles to the output directory...")

    files = os.listdir('Particles/%s' % outDir)
    for f in files:
        os.rename('Particles/%s/%s' % (outDir, f), '%s/%s' % (output, f))
    os.rmdir('Particles/%s' % outDir)



def write_output_starfiles(labels, mdOut, mdOutSub, output):

    labels.extend(['rlnCoordinateX', 'rlnCoordinateY', 'rlnMicrographName'])

    print("\nWriting output STAR files.")

    starfile1 = output + ".star"
    print(" Subparticles (without subtraction):\t\t%s" % starfile1)
    # We convert back angles to degrees and write subparticles star file
    def _writeMd(md, starfile):
        for subpart in md:
            angles_to_degrees(subpart)
        md.addLabels(labels)
        md.write(starfile)

    _writeMd(mdOut, starfile1)

    if len(mdOutSub):
        starfile2 = starfile1.replace('.star', '_subtracted.star')
        print(" Subparticles (with subtraction):\t\t%s" % starfile2)
        _writeMd(mdOutSub, starfile2)


def split_star_to_random_subsets(inputStarRoot):
    inputStarName = inputStarRoot+'.star'
    md = MetaData(inputStarName)
    mdHalf1 = MetaData()
    mdHalf2 = MetaData()

    half1StarRoot = inputStarRoot+'_half1'
    half2StarRoot = inputStarRoot+'_half2'

    half1StarName = half1StarRoot+'.star'
    half2StarName = half2StarRoot+'.star'

    particlesHalf1 = []
    particlesHalf2 = []

    for particle in md:
        if particle.rlnRandomSubset % 2 == 1:
            particlesHalf1.append(particle.clone())
        if particle.rlnRandomSubset % 2 == 0:
            particlesHalf2.append(particle.clone())

    mdHalf1.addData(particlesHalf1)
    mdHalf2.addData(particlesHalf2)
    labels = md.getLabels()
    mdHalf1.addLabels(labels)
    mdHalf2.addLabels(labels)
    mdHalf1.write(half1StarName)
    mdHalf2.write(half2StarName)

    return half1StarRoot, half2StarRoot


def reconstruct_subparticles(threads, output, maxres, sym, angpix, library_path):
    """ Reconstruct subparticles. Also create two half maps using random subsets. """

    def run_reconstruct(input, suffix='', extraArgs=''):
        if threads>1 :
            cmd = ('mpirun -np %s relion_reconstruct_mpi ') % (threads)
        else:
            cmd = ('relion_reconstruct ')

        args = ('--sym %s %s --o %s%s.mrc --i %s.star --angpix %s') % (sym, extraArgs, input, suffix, input, angpix)
        run_command(cmd + args,"", library_path)

    for input in [output, output+'_subtracted']:
        args = ""
        inputStarName = input + '.star'
        if os.path.exists(inputStarName):
            md = MetaData(inputStarName)

            if "rlnDefocusU" in md.getLabels():
                args = args + "--ctf "
            else:
                print ("\nWarning: no CTF info found in %s!\n"
                       "The reconstruction will be performed without CTF correction.\n" % inputStarName)

            # reconstruct random halves to Nyquist frequency
            if "rlnRandomSubset" in md.getLabels():
                half1Star, half2Star = split_star_to_random_subsets(input)
                run_reconstruct(half1Star, "_class001_unfil", args)
                run_reconstruct(half2Star, "_class001_unfil", args)

            # reconstruct the map to maximum resolution given
            if maxres:
                args = args + "--maxres %s " % maxres
            run_reconstruct(input, '_class001', args)


def run_command(command, output, library_path):
    import subprocess

    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = env['LD_LIBRARY_PATH']+":"+library_path

    if not output:
        print("+++ " + command)
        sys.stdout.flush()
        proc = subprocess.Popen(command, shell=True, env=env)
        proc.wait()
    else:
        print("+++ " + command)
        sys.stdout.flush()
        with open(output, "wb", 0) as out:
            proc = subprocess.Popen(command, shell=True, stdout=out, env=env)
            proc.wait()

def unique_micrographs(md):
    particles = []
    seen = set()
    unique_micrographs = []

    for particle in md:
            particles.append(particle)

    for particle in particles:
        if particle.rlnMicrographName not in seen:
            unique_micrographs.append(particle)
            seen.add(particle.rlnMicrographName)
    return unique_micrographs

class ProgressBar():
    """ Implements a simple command line progress bar
    with a big fish catching a small fish.
    Works nicely only if the number of iterations is larger than the width. """

    def __init__(self, width, total):
        # hide cursor
        sys.stdout.write("\033[?25l")
        # setup progressbar
        self.width = width
        sys.stdout.write("%s><^>" % ("~" * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * width)
        self.count = 0  # total count
        self.c = 0  # progress count
        self.n = 0 # total count
        self.percent = 1.0/width
        self.timer = self.percent
        self.total = total

    def notify(self):
        self.n += 1
        if self.count == int(self.total * self.timer):
            sys.stdout.write("\b" * (self.c + 8))
            sys.stdout.write("~" * self.c)
            sys.stdout.write("><))^)>")
            sys.stdout.flush()
            self.timer += self.percent
            self.c += 1
        if self.n == int(self.total):
            # restore cursor
            sys.stdout.write("\033[?25h")
            sys.stdout.flush()

        self.count += 1
