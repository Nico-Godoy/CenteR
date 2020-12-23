import pynpoint
from pynpoint import Pypeline
from pynpoint import Hdf5ReadingModule
from pynpoint import PSFpreparationModule
from pynpoint import PcaPsfSubtractionModule
from matplotlib import pyplot as plt


pipeline = Pypeline(PATH1,
                    PATH2,
                    PATH1)

reading_dict = {"im_arr": "im_arr"}

reading = Hdf5ReadingModule(name_in="hdf5_reading",
                            tag_dictionary=reading_dict)

pipeline.add_module(reading)


preparation = PSFpreparationModule(name_in="prep",
                                   image_in_tag="im_arr",
                                   image_out_tag="prep",
                                   mask_out_tag=None,
                                   norm=False,
                                   resize=None,
                                   cent_size=RADIUS,
                                   edge_size=None)

pipeline.add_module(preparation)

subtraction = PcaPsfSubtractionModule(COMP, 
                                   name_in="PSF_subtraction",
                                   images_in_tag="prep",
                                   reference_in_tag="prep",
                                   res_median_tag="result")

pipeline.add_module(subtraction)

pipeline.run()

result = pipeline.get_data("result")

#pixscale = pipeline.get_attribute("stack", "PIXSCALE")


