
Description of core files:
----------

``CenteR`` core files used in the main code sequence:

- Processing_Data_v10.R / Processing_Data_v11.R
   Program to correct the files by cosmetics. The program do not recognize both photometric and science observations with more than 2 different exposure times.
   The program follow the next path: reduced the data for cosmetics, construct a master file containing the information of the FITS files, compute an aproximation of the circular aperture center, compute the sky subtraction (frames closest in time are selected for), and finish with the pseudo-Zernike moments decomposition (noralize and unormalize torus S/N and Homogeneity), constructing a final file called 'Total_Final_Info_Red.txt'


- Fitting_2D_Nor_v8.R
   Contain the program to obtain both the AGPM and circular aperture centers in the sky frames.
   In addition, determine the center of the circular aperture in the science frames.
   - The input are the FITS files arealdy corrected by current dark, flat field, bad pixels, both with and without sky subtraction, and the text file with the summary of the FITS files.
   - The outputs are a serie of text files with the information about the fitting for the both cases.

- Finding_Star_Pos_Parallel_v12.R / Finding_Star_Pos_Parallel_v12b.R
   Program to fit the negative and positive 2D gaussians, modeling the torus shape. 
   - The inputs are all the outputs of previous steps.
   - The output corresponds to the file 'Final_table_complete_information_v7.txt', containing the information of the science frames, the modeling and the result of the fitting.
   The code require a certain numbers of cores, due to the process is optimized and wroten to run in parallel processes.


External codes to compute the PCA using PynPoint. These codes are constructed to use the full dataset with different combination of centers:
- PIC_v12_mod_DB_coro.R
   Use in both the PCA and de-rotate and stack the frames, the AGPM center
- PIC_v12_mod_DB_coro__STACK.R
   The same than the previous one, but stacking sequence of frames.
- PIC_v12_mod_DB_only_star.R
   Use in both cases the position of the star.
- PIC_v12_mod_DB_star.R
   Use first the AGPM position to performce PCA, and then the star location to de-rotate and stack all the frames.

   External codes to compute the PCA using PynPoint. These codes are constructed to use frame selection with different combination of centers:
- PI_v14_DB_coro.R
   Using for the PCA and de-rotate and stack the frames the AGPM center.
- PI_v14_DB_only_star.R
   Using for the PCA and de-rotate and stack the frames the star location.
- PI_v14_DB_star.R
   Using for the PCA the AGPM center an then to de-rotate and stack the frames the star location.


   The 'Master_file.R' corresponds to an example of how to use the codes together, and the variables that must be specified for.


Requiriments:
----------

``CenteR`` require the previous instalation of the next:
- R-cran installed
- Python 3 installed

R-cran can be installed follow the instructions in the link: https://cran.r-project.org/

For example, using the terminal in Linux should be:

.. code-block:: bash

    $ pip install vip_hci



.. code-block:: bash

  $ sudo apt update
  $ sudo apt install r-base
  $ sudo apt install build-essential



for Python, you can see the download/instalation here: https://www.python.org/downloads/

CenteR R-cran codes need the next packages:
- FITSio
- astro
- fields
- astroFns
- abind
- IM
- mvtnorm
- plot3D
- doParallel
- h5
- rPython

It is possible to install the packages directly from R-cran or downloading the source-package.
From R-cran, open a terminal and then write R + ENTER:
You will see in the terminal something like:

.. code-block:: R
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.
  Natural language support but running in an English locale
R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.
Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

Then copy the follow commands one by one:


install.packages('FITSio', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('astro', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('fields', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('astroFns', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('abind', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('IM', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('mvtnorm', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('plot3D', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('doParallel', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('h5', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('rPython', dependencies=TRUE, repos='http://cran.rstudio.com/')

or via source-package, downloading from: https://cran.r-project.org/web/packages/available_packages_by_name.html
Then, writen in R-cran terminal:
install.packages(path_to_file, repos = NULL, type="source")
with path_to_file the source of the respective pakcage.

some packages (astro, h5, doParallel) could have some problems with the instalation.
In general, using the source-package solves the problem, or using:
install.packages('astro')
then, selecting the 'old' mirror [0] in the displayed window.
There are more options to install the packages. Always R-cran will indicate the problem in the terminal, for example, a specific packages in needed for the instalation of one specific package. In that case, you need to install this before continue.

From Python, you need to install PynPoint. You can follow the steps from the oficial webpage: https://pynpoint.readthedocs.io/en/latest/installation.html
In addition, you must have installed:
- matplotlib
- ephem
- numpy
You can install in Linux, for exmaple, using the following commands:
sudo pip install matplotlib ephem numpy

The code use therminal commands, so it is strongly recommended to run CenteR pipeline in Linux.

