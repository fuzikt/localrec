# Localized Reconstruction

Localized reconstruction of flexible subunits from macromolecular complexes.

Version: 1.3.0

### New features added

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

--library_path  
define LD_LIBRARY_PATH used by Relion. Default: empty


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
