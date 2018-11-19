#!/bin/bash

export FHICL_DIR=${MRB_SOURCE}/analysistree/analysistree
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
cd -
export FHICL_FILE_PATH=${FHICL_DIR}:${FHICL_FILE_PATH}
