#!/usr/bin/env python3

# **************************************************************************
# *
# * Authors:  Serban Ilca (serban@strubi.ox.ac.uk)
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *           Tibor Fuzik (tibor.fuzik@ceitec.muni.cz)
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
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
from lib.localrec import *


class CreateSymmetryRelatedParticles():
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Rotates all particles in the input STAR file.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('input_star', help="Input STAR filename with particles.")
        add('--vector', default=None, help="Vector defining the additional rotation of the particles (x, y, z).")
        add('--angles', default=[],
            help="Euler angles defining the additional rotation of the particles (rot, tilt, psi).")
        addr('--output', required=True, help="Output STAR filename.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if len(sys.argv) == 1:
            self.error("Error: No input file given.")

        if not os.path.exists(args.input_star):
            self.error("Error: Input file '%s' not found."
                       % args.input_star)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        print("Creating rotated particles...")

        md = MetaData(args.input_star)

        if args.vector and len(args.angles) != 0:
            print("Please only provide a vector or a triplet of Euler angles for the particle rotation.")
            sys.exit(0)
        elif args.vector:
            vectors = load_vectors(None, args.vector, None, 1)
            rot_matrix = vectors[0].matrix()
        elif len(args.angles) != 0:
            angles = args.angles.split(',')
            if len(angles) != 3:
                print("Please provide exactly 3 Euler angles for the particle rotation.")
                sys.exit(0)
            else:
                rot_rot = math.radians(float(angles[0]))
                rot_tilt = math.radians(float(angles[1]))
                rot_psi = math.radians(float(angles[2]))
                rot_matrix = matrix_from_euler(rot_rot, rot_tilt, rot_psi)
        else:
            print("Please provide a vector or a triplet of Euler angles for the particle rotation.")
            sys.exit(0)

        new_particles = []

        for particle in md:
            new_particle = particle.clone()

            angles_to_radians(particle)

            rot = particle.rlnAngleRot
            tilt = particle.rlnAngleTilt
            psi = particle.rlnAnglePsi

            matrix_particle = matrix_from_euler(rot, tilt, psi)

            m = matrix_multiply(matrix_particle, matrix_transpose(rot_matrix))
            rotNew, tiltNew, psiNew = euler_from_matrix(m)

            new_particle.rlnAngleRot = rotNew
            new_particle.rlnAngleTilt = tiltNew
            new_particle.rlnAnglePsi = psiNew

            angles_to_degrees(new_particle)
            new_particles.append(new_particle)

        if md.version == "3.1":
            mdOut = md.clone()
            dataTableName = "data_particles"
            mdOut.removeDataTable(dataTableName)
        else:
            mdOut = MetaData()
            dataTableName = "data_"

        mdOut.addDataTable(dataTableName, md.isLoop(dataTableName))
        mdOut.addLabels(dataTableName, md.getLabels(dataTableName))
        mdOut.addData(dataTableName, new_particles)

        mdOut.write(args.output)

        print("All done!")
        print(" ")


if __name__ == "__main__":
    CreateSymmetryRelatedParticles().main()
