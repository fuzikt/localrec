# Localized Reconstruction

Localized reconstruction of flexible subunits from macromolecular complexes.

Version: 1.5.0

RELION 3.1+ compatible particle extraction.

NO Scipion needed to run!!!

### New features added

- compatible with Relion 3.1 star files, backward compatible with Relion 3.0 star files (be aware to use the right version of relion for extraction)

- for Relion 3.1 star files automatically reads the apix from star file

- for Relion 3.1 star files binned particle star file can be used for extraction from micrographs (for Relion 3.0 star files the star file needs to be "unbinned first")

Changes in 1.5.0 version
- Originating particle of the subparticle is stored in rlnOriginalParticle in the output star file.
- The run command with all parameters is stored in the output star file as comment at the beginning of the file.

New options for extraction:

--rescale_size
Rescale (bin) subparticle to this box (pixels). (Default: -1 => no rescale)

--float16
Use float16 format for the output MRC files.

--no_ramp
Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background."


Changes in 1.3.0 version:

RELION 3 compatible particle extraction.
NO Scipion needed to run!!!

--extract_from_micrographs  
Extract subparticles directly form micrographs instead from particle images. Cannot be used with signal subtraction.

--only_extract_unfinished  
Extract only particles not extracted in previous run.

--invert_contrast  
Use this option when extracting from micrographs that were not inverted (i.e. black particles ons on white background).

--normalize  
Normalize the extracted particles. Useful when extracting from micrographs.

--do-halves
Make reconstruction of the random halves only if this option is set.

--library_path  
define LD_LIBRARY_PATH used by Relion. Default: empty

--particle_size
!!! not used anymore. Autodetected from particle/micrograph image. Support for non-square micrographs (e.g. K2, K3)


### It also requires the following software:
* Relion 3

## Reference

If localized reconstruction is useful in your work, please cite:

Ilca SL, Kotecha A, Sun X, Poranen MM, Stuart DI & Huiskonen JT (2015).
Localized reconstruction of subunits from electron cryomicroscopy images of macromolecular complexes.
Nat Commun 6, 8843. [doi:10.1038/ncomms9843](http://dx.doi.org/10.1038/ncomms9843)

## Instructions

See: https://github.com/OPIC-Oxford/localrec/wiki

## Acknowledgements

The authors acknowledge funding from the European Research Council under the European Union’s Horizon 2020 research and innovation programme (649053) and from the Wellcome Trust PhD Programme in Structural Biology.
