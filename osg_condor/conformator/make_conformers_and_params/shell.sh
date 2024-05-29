#!/bin/bash

set -e 

tar -xzvf split_new_named_$3.sdf.tar.gz
python make_conformator_and_params.py  split_new_named_$3.sdf
