library(FITSio)
library(astro)
library(fields)
#library(rgl)
library(astroFns)
library(abind)
library(IM)
library(mvtnorm)
library(plot3D)
library(doParallel)
library(h5)

# Path to the directory that contain the CenteR core files
path.pip<-"/path/to/centeR/"
# Path to the raw data directory
path.data2<-"/data/with/raw/calibs/"
# Path to the work-area directory, where CenteR will generates all the outputs and pre-processes files
work.area2<-"/work/area/"

# Coordinates of the target (for example, in ICRS J2015). This coordinates are used to validate the
# identification of files (sky/science), or to identify both (in cases of older FITS files)
ra2<- "00:00:00.0000" # RA
dec2<-"-00:00:00.0000" # DEC

path.data<-path.data2
work.area<-work.area2

################################################################################
####################           CenteR codes path:           ####################
################################################################################

#-------------------------------------------------------------------------------
# Inputs:

# Pre-processing dataset
   Data_reduction<-"TRUE"  # TRUE or FALSE, with capital letters: correction for master dark, master flat and bad pixels
     NODES<-20             # from 1 up to the number of availables cores for the parallel process
     OLD<-"FALSE"          # just in particular cases you can change this to TRUE
   Master_File<-"TRUE"     # TRUE or FALSE. Generate the master file with the information about the FITS files
   Center_File<-"TRUE"     # TRUE or FALSE. Estimate the approximation location of the circular aperture center
   Sky_Subs<-"TRUE"        # TRUE or FALSE. Compute the sky subtraction, saving new FITS files WITHOUT previous headers.
   EvA_File<-"TRUE"        # TRUE or FALSE. Compute the pseudo-Zerkine moments decomposition and calculate the torus S/N and H.

 Proc_Data<-"TRUE"         # TRUE or FALSE. Run the entire sequence of codes above or not. 

# Calibrating the AGPM detector location 
 # The program find the AGPM center in the sky frames, and the circular aperture center in both science and sky frames
 Fit_Nor<-"TRUE"           # TRUE or FALSE. 

# Finding the star position
 # Fit the negative and positive 2D gaussians to model the torus shape.
   Cluster_number<-30      # from 1 up to the number of availables cores for the parallel process

 Star_Pos<-"TRUE"         # TRUE or FALSE

#-------------------------------------------------------------------------------

# Every step of CenteR code will or will not run depending on the defined parameters above.
# Please check the proper version of the code for every step fo the program.  
if(Proc_Data=="TRUE"){ source(paste0(path.pip,"Processing_Data_v11.R"))}
if(Fit_Nor=="TRUE"){ source(paste0(path.pip,"Fitting_2D_Nor_v7.R"))}
if(Star_Pos=="TRUE"){ source(paste0(path.pip,"Finding_Star_Pos_Parallel_v12b.R"))}

################################################################################

################################################################################
###########     External code to compute the final image via PCA     ###########
################################################################################


# PCA using a combination of codes between PynPoint, Python and R-cran

# Below an example of inputs for the 7 cases and files (4 for the full dataset and 3 for the frame selection).
# The complete dataset and the frame selection are separated in 3 cases, depedning the combination of centers:
# - only using the AGPM for both the PCA and de-rotate and stack the datacube
# - using the AGPM for the PCA, and then the star to de-rotate and stack the frames
# - only using the star location for both the PCA and de-rotate and stack the datacube


#################-----------      Full dataset      -----------#################

# only AGPM
   Cluster_number<-20 # from 1 up to the number of availables nodes/cores for the parallel process used to construct the h5 file.
   Cluster_number_p2<-10 # number of nodes/cores to use in the parallel process in the PCA
   Object_hdf5_tot<-"Object_file__CT_39x39-TF_on-AGPM.hdf5"# name of the H5 file.
   Xl<-39  # Dimention of the cube along the X-axis (the cube will be cropped following this value)
   Yl<-39  # Dimention of the cube along the Y-axis. Both must be an odd number.
   re<-0.5 # The rho value for the frame selection. Request in the program but not used
   SiG<-3  # The phi value for the frame selection. Request in the program but not used
   Radius_pix<-2  # The radius of the inner mask, in pixels
   Re_Size<-"FALSE"  # TRUE or FALSE. Actually, this is use only for internal tests. Please set FALSE
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100 # % pixels w/r size image
   HDF5_file_tot<-"TRUE" # TRUE or FALSE. Construct the HDF5 file with the frames and PA.
   PAng<-"New"  # Old or New. 
                   # PAng='Old' compute the Parallactic Angle just taking the first and last value in the header
                   # and then generating a sequence of values to providea PA to every frame in the FITS file
                   # PAng='New' compute the PA in  ... 
   HDF5_folder_tot<-"HDF5_file__CT_39x39-TF_CORO_" # Folder name to save the final products
   HDF5_path<-'/debris_data/ngodoy/PynPoint_products/'  # Path where the PCA will run. It is recommended do not run the parallel
                   # process in a NAS, for example.
   COMP<-unique(c(seq(1,55, by=1)))  # Number of principal components to compute for.
   PCA_run_tot<-"TRUE"  # TRUE or FALSE, run the PCA code.

  PI_ac<-"TRUE" # Run the codes (create the h5 file and the PCA procedure)

if(PI_ac=="TRUE"){ source(paste0(path.pip,"PIC_v12_mod_DB_coro.R"))}
# PIC_v12_mod_DB_coro.R is the file that use only the AGPM position as center

# star+AGPM
   Cluster_number<-20
   Cluster_number_p2<-10
   Object_hdf5_tot<-"Object_file__CT_39x39-TF_on-AGPM.hdf5"
   Xl<-39
   Yl<-39
   re<-0.5
   SiG<-3
   Radius_pix<-2
   Re_Size<-"FALSE"
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100 # % pixels w/r size image
   HDF5_file_tot<-"FALSE"
   PAng<-"New"
   HDF5_folder_tot<-"HDF5_file__CT_39x93-TF_STAR_"
   COMP<-unique(c(seq(1,55, by=1)))
   PCA_run_tot<-"TRUE"

  PI_ac<-"TRUE"

if(PI_ac=="TRUE"){ source(paste0(path.pip,"PIC_v12_mod_DB_star.R"))}
# PIC_v12_mod_DB_star.R is the file that use both the star and AGPM centers


# PIC_v8.R ## Preparing_Input (PI)
   Cluster_number<-20
   Cluster_number_p2<-10
   Object_hdf5_tot<-"Object_file__CT_39x39-TF_on-STAR.hdf5"
   Xl<-39
   Yl<-39
   re<-0.5
   SiG<-3
   Radius_pix<-2
   Re_Size<-"FALSE"
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100 # % pixels w/r size image
   HDF5_file_tot<-"TRUE"
   PAng<-"New"
   HDF5_folder_tot<-"HDF5_file__CT_39x39-TF_onlySTAR_"
   COMP<-unique(c(seq(1,55, by=1)))
   PCA_run_tot<-"TRUE"

  PI_ac<-"TRUE"

if(PI_ac=="TRUE"){ source(paste0(path.pip,"PIC_v12_mod_DB_only_star.R"))}
# PIC_v12_mod_DB_only_star.R is the file that use only the star position



##############----------        Frame selection        ----------###############

# Only AGPM
   Cluster_number<-20 # from 1 up to the number of availables nodes/cores for the parallel process used to construct the h5 file.
   Cluster_number_p2<-10 # number of nodes/cores to use in the parallel process in the PCA
   Object_hdf5<-paste0("Object_file__FS-CT-TF_39x39_on-AGPM.hdf5") # name of the H5 file.
   Xl<-39  # Dimention of the cube along the X-axis (the cube will be cropped following this value)
   Yl<-39  # Dimention of the cube along the Y-axis. Both must be an odd number.
   re<-0.5 # The rho value for the frame selection. We recommend to use re=0.5
   SiG<- 0 # The phi value for the frame selection. We recommend values between -1 to +1
   Re_Size<-"FALSE" # TRUE or FALSE. Actually, this is use only for internal tests. Please set FALSE
   Radius_pix<-2 # The radius of the inner mask, in pixels
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100  # % pixels w/r size image
   HDF5_file<-"TRUE" # TRUE or FALSE. Construct the HDF5 file with the frames and PA.
   PAng<-"New"  # Old or New. 
                   # PAng='Old' compute the Parallactic Angle just taking the first and last value in the header
                   # and then generating a sequence of values to providea PA to every frame in the FITS file
                   # PAng='New' compute the PA considering the rotation of the sky and more parameters in the header.
                   # This one is more accurate. 
   HDF5_folder<-paste0('HDF5_file__FS-CT_39x39-TF_CORO_') # Folder name to save the final products
   HDF5_path<-'/debris_data/ngodoy/PynPoint_products/'  # Path where the PCA will run. It is recommended do not run the parallel
                   # process in a NAS, for example.
   COMP<-unique(c(seq(1,55, by=1)))  # Number of principal components to compute for.
   PCA_run<-"TRUE" # TRUE or FALSE, run the PCA code.

  PI_a<-"TRUE" # Run the codes (create the h5 file and the PCA procedure)

if(PI_a=="TRUE"){ source(paste0(path.pip,"PI_v14_DB_coro.R"))}
# PI_v14_DB_coro.R is the file that use only the AGPM position as center

# AGPM + star
   Cluster_number<-20
   Cluster_number_p2<-10
   Object_hdf5<-paste0("Object_file__FS-CT_39x39_on-AGPM_.hdf5") 
   Xl<-39  
   Yl<-39 
   re<-0.5 
   SiG<- 0 
   Re_Size<-"FALSE" 
   Radius_pix<-2 
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100 
   HDF5_file<-"FALSE" # TRUE or FALSE. Construct the HDF5 file with the frames and PA.
            # If you already construct the h5 file (for example, when you run the code for 'only AGPM'), you should set FALSE
   PAng<-"New" 
   HDF5_folder<-paste0("HDF5_file__FS-CT_39x39-TF_STAR_") 
   HDF5_path<-'/debris_data/ngodoy/PynPoint_products/'  
   COMP<-unique(c(seq(1,55, by=1))) 
   PCA_run<-"TRUE"

  PI_a<-"TRUE" 

if(PI_a=="TRUE"){ source(paste0(path.pip,"PI_v14_DB_star.R"))} 
# PI_v14_DB_star.R is the file that use both the star and AGPM centers


# Only star
   Cluster_number<-20
   Cluster_number_p2<-10
   Object_hdf5<-paste0("Object_file__FS-CT-TF_39x39_on-STAR.hdf5")
   Xl<-39
   Yl<-39
   re<-0.5
   SiG<- AuxV[IUp]
   Re_Size<-"FALSE"
   Radius_pix<-2
   Mask_radius<- round((Radius_pix*27.2/1000)*100)/100  # % pixels w/r size image
   HDF5_file<-"TRUE"
   PAng<-"New"
   HDF5_folder<-paste0("HDF5_file__FS-CT_39x39-TF_onlySTAR_")
   COMP<-unique(c(seq(1,55, by=1)))
   PCA_run<-"TRUE"

  PI_a<-"TRUE"

if(PI_a=="TRUE"){ source(paste0(path.pip,"PI_v14_DB_only_star.R"))}
# PI_v14_DB_only_star.R is the file that use only the star position






