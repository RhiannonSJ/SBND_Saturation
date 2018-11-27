#!/bin/bash

OUTPUT_DIR=${MRB_SOURCE}/saturation/saturation/output_dir
WORKING_DIR=${MRB_SOURCE}/saturation/saturation

TEST=/pnfs/sbnd/persistent/users/rsjones/g4_sam_tests/g4-cdffe4ba-0bb0-426d-8901-8e0fc17883e9.root
TEST_CORSIKA=/pnfs/sbnd/mc/g4/artroot/pre-production/MCP0.9/prodoverlay_corsika_cosmics_cmc_genie_nu_spill_gsimple-configd-v1_tpc/sbndcode/v07_07_00_2_MCP0_9/run_number/00/00/00/01/subrun_number/00/00/00/53/g4-744cbbe3-8d8e-4e0e-ba56-3c98a4f0188a.root
TEST_LIST=/pnfs/sbnd/scratch/users/rsjones/g4_file_list.txt
TEST_RECO_LIST=/pnfs/sbnd/scratch/users/rsjones/reco_file_list.txt
TEST_CORSIKA_LIST=/pnfs/sbnd/scratch/users/rsjones/corsika_reco_files.txt

rm -rf $OUTPUT_DIR/*

cd $OUTPUT_DIR

###lar -c run_saturation.fcl -s $TEST
lar -c run_saturation.fcl -s $TEST_CORSIKA
###lar -c run_saturation.fcl -S $TEST_LIST
###lar -c run_saturation.fcl -S $TEST_RECO_LIST

cd $WORKING_DIR
