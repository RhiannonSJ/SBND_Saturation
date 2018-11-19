#!/bin/bash

OUTPUT_DIR=${MRB_SOURCE}/saturation/saturation/output_dir
WORKING_DIR=${MRB_SOURCE}/saturation/saturation

TEST=/pnfs/sbnd/persistent/users/rsjones/g4_sam_tests/g4-cdffe4ba-0bb0-426d-8901-8e0fc17883e9.root
TEST_LIST=/pnfs/sbnd/scratch/users/rsjones/g4_file_list.txt
TEST_RECO_LIST=/pnfs/sbnd/scratch/users/rsjones/reco_file_list.txt

rm -rf $OUTPUT_DIR/*

cd $OUTPUT_DIR

###lar -c run_saturation.fcl -s $TEST
###lar -c run_saturation.fcl -S $TEST_LIST
lar -c run_saturation.fcl -S $TEST_RECO_LIST

cd $WORKING_DIR
