import pynpoint
from pynpoint import Pypeline
from pynpoint import Hdf5ReadingModule
from pynpoint import PcaPsfSubtractionModule
from matplotlib import pyplot as plt
from pynpoint import FalsePositiveModule

pipeline = Pypeline(PATH1,
                    PATH2,
                    PATH1)

reading_dict = {TAG_IM: TAG_IM}

reading = Hdf5ReadingModule(name_in="hdf5_reading",
                            tag_dictionary=reading_dict)

pipeline.add_module(reading)

SNR_mod = FalsePositiveModule( position=( XPOS , YPOS),
                     aperture= APERT_RAD,
                     ignore=True,
                     image_in_tag=TAG_IM,
                     name_in="SNR_mod",
                     snr_out_tag="SNR_result")


pipeline.add_module(SNR_mod)

pipeline.run()

resutl = pipeline.get_data("SNR_result")


