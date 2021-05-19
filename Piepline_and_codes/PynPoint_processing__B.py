import pynpoint
from pynpoint import Pypeline
from pynpoint import Hdf5ReadingModule
from pynpoint import PSFpreparationModule
from pynpoint import PcaPsfSubtractionModule
from pynpoint import DerotateAndStackModule
from matplotlib import pyplot as plt


pipeline = Pypeline(PATH1,
                    PATH2,
                    PATH1)

reading_dict = {"for_derotate": "for_derotate"}

reading = Hdf5ReadingModule(name_in="hdf5_reading",
                            tag_dictionary=reading_dict)

pipeline.add_module(reading)


derotating = DerotateAndStackModule(name_in="for_derotate",
	                            image_in_tag="for_derotate",
                                    image_out_tag="derotated_im",
                                    derotate=True,
                                    stack="median")

pipeline.add_module(derotating)


pipeline.run()

result = pipeline.get_data("derotated_im")

#pixscale = pipeline.get_attribute("stack", "PIXSCALE")


