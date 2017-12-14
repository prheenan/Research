# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../../../../")
sys.path.append("../../")
from Util import IoUtil


def run(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    input_dir =  IoUtil.data_dir(in_dir)
    output_dir = IoUtil._traces_dir(in_dir)
    IoUtil.read_images(input_dir,cache_dir=output_dir,force=Truew)
    
if __name__ == "__main__":
    run(IoUtil.get_directory_command_line())
