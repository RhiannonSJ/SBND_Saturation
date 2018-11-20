#!/bin/bash

OUTPUT_DIR=${MRB_SOURCE}/saturation/saturation/output_dir

TEST=/pnfs/sbnd/persistent/users/rsjones/g4_sam_tests/g4-cdffe4ba-0bb0-426d-8901-8e0fc17883e9.root
TEST_CORSIKA=/pnfs/sbnd/mc/reco/artroot/pre-production/MCP0.9/prodoverlay_corsika_cosmics_cmc_genie_nu_spill_gsimple-configd-v1_tpc/sbndcode/v07_07_00_2_MCP0_9/run_number/00/00/00/01/subrun_number/00/00/00/58/reco-cd71ec1d-d2f6-4ca8-b147-e8a8e33c8b50.root
TEST_LIST=/pnfs/sbnd/scratch/users/rsjones/g4_file_list.txt
TEST_RECO_LIST=/pnfs/sbnd/scratch/users/rsjones/reco_file_list.txt

rm -rf $OUTPUT_DIR/*

cd $OUTPUT_DIR

###lar -c run_saturation.fcl -s $TEST
lar -c run_saturation.fcl -s $TEST_CORSIKA
###lar -c run_saturation.fcl -S $TEST_LIST
###lar -c run_saturation.fcl -S $TEST_RECO_LIST

cd $WORKING_DIR
