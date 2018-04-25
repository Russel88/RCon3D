RCon3D: Analyzing 3D confocal images of microbial biofilms/communities
----------------------------------------------------------------------

This is an R package for various 3D analyses on confocal images of
microbial biofilms, microcolonies and communities.

### Citation

If you use functions in the notebooks **Loading and preparing images**
or **Quantification of images** please cite: [Liu et al. (2017)
Low-abundant species facilitates specific spatial organisation that
promotes multispecies biofilm formation. *Envir.
Microbiol.*](http://onlinelibrary.wiley.com/doi/10.1111/1462-2920.13816/abstract)

If you use functions in the notebook **Co-aggregation and relative
positioning** please cite: [Liu et al. (2018) Micro-scale intermixing; a
requisite for stable and synergistic co-establishment in a four-species
biofilm. *ISMEJ*](https://rdcu.be/LSda)

### Install package

    # devtools has to be installed (install.packages("devtools"))
    devtools::install_github("Russel88/RCon3D")

### Notebooks on different analyses:

[Loading and preparing images, required for all
analyses](https://nbviewer.jupyter.org/github/Russel88/RCon3D/blob/master/Notebooks/Loading.html)

[Quantification of
images](https://nbviewer.jupyter.org/github/Russel88/RCon3D/blob/master/Notebooks/Quant.html)

[Co-aggregation and relative
positioning](https://nbviewer.jupyter.org/github/Russel88/RCon3D/blob/master/Notebooks/Coagg.html)

[Dynamic sectioning of images in
z-direction](https://nbviewer.jupyter.org/github/Russel88/RCon3D/blob/master/Notebooks/Section.html)

### Multi-threading

Most functions can run in parallel by setting the `cores` argument. This
will highly speed up analyses, but will also use more memory.

### All functions

`loadIMG` Load .tif files and turn them into arrays, and save them as
RDS files ready for downstream analysis

`findIMG` Find already loaded images in a set directory

`smoothIMG` Median smooth images

`morphIMG` Apply mathematical morphology kernel on image

`merge_channels` Create arrays based on combinations of channels. E.g.
intersects, unions or subtractions

`quant` Quantifies the number of pixels in each image, for each channel
at each layer.

`layer_stand` Standardize layers based on fill. Relevant if the bottom
of the specimen is the layer with highest fill

`layer_split` Splits quantification in Top, Middle and Bottom based on
fill and/or a set number of layers

`occupancy` Estimates the proportion a (target) channel occupy around a
(focal) channel

`co_agg` Estimates the co-aggregation between two channels. An
undirected version of `occupancy`

`cross_ratio` Estimates the ratio between two channels (targets), at
some distance from a focal channel

`create_random` Create random images for testing

`create_nulls` Create image files for the empty spaces in an image, such
that this can be used in the analysis. For example, it can then be
calculated how much empty space there is a around a certain channel with
`occupancy`

`extract_layers` With the output from `layer_split` it makes a list of
what layers to analyze in `occupancy`, `clumps`, `co_agg` and
`cross_ratio`

`xy_splits` Splits the image in each xy position and then run an
analysis

#### Acknowledgment note:

An internal function, `tiff_to_array`, is partly borrowed from
<https://github.com/rmnppt/iMage>. Furthermore, some of the algorithmic
framework for the `co_agg`,`occupancy` and `cross_ratio` analysis is
also borrowed from this repository.
