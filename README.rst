``CenteR`` - centering and frame selection
======================================

.. image:: https://img.shields.io/badge/R--cran-3.5.1-brightgreen
    :target: https://www.r-project.org/

.. image:: https://img.shields.io/badge/license-MIT-blue
    :target: https://github.com/Nioc-Godoy/CenteR/main/LICENSE

Introduction
------------

``CenteR`` is a code designed for VLT/NACO data reduction using vortex cononagraph  AGPM. The pipeline was designed initially to reduced NACO-ISPY observations, but it should work with any similar observing strategy. It is based on a new centering technique using the images themselves for the detection of the AGPM in thermal emission, to later determine the position of the stars behind the coronagraph. In addition, it uses a new approach for frame selection based on the parameters obtained during centering and a decomposition based on pseudo-Zernike moments.

Requirements
------------
``CenteR`` is based on R-cran, and use PynPoint for the principal component analysis and de-rotate and stack the frames. ``CenteR`` pipeline provides the registration adn celection of frames, and both the star location and the AGPM position. Different stategies are carry out for the principal component analysis and de-roate and stack the frames, using different centers in both steps, and the use of the entire datacube or the option of using frame selection.
The R-cran, Python and PynPoint instalation and the R-cran packages instalation are summary in Pipeline_and_codes/README_files.txt 

Using ``CenteR``
------------
``CenteR`` was tested in ISPY observations and other NACO obsertations with similar observing strategies. However, the right working of the pipeline can be affected by the signal of the frames, windowing and strategies implemented. 

An explame of how to run ``CenteR`` and the accopled PynPoint from python is presented at Pipeline_and_codes/Master_file.R, with the necessary inputs and definitions necesary for. The pipeline works running this file and do not require any human intervetion during the processing, just the inition inputs. the code can be run open R from the terminal writing R+ENTER, and writing source('path_to_/Master_file.R').

