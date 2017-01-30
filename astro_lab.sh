#!/bin/bash

export PYTHONPATH=/user/phstf/md0046/tools/photutils/build/lib.linux-x86_64-2.7/:/user/phstf/md0046/tools/pyds9/build/lib.linux-x86_64-2.7:$PYTHONPATH
export PATH=/user/phstf/md0046/tools/xpa/build/bin:$PATH #/user/phstf/md0046/tools/ds9:$PATH
export LD_LIBRARY_PATH=/user/phstf/md0046/tools/libs:$LD_LIBRARY_PATH

ipython -i astro_lab.py
