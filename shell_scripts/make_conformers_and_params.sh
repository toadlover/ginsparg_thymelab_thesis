#!/bin/bash

#this is a simplified script that is used to create conformers and Rosetta params files
#this script makes 

set -e 

tar -xzf split_new_named_$1.sdf.tar.gz
python make_conformator_and_params.py  split_new_named_$3.sdf
