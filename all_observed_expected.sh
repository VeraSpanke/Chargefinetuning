#!/usr/bin/env bash
#run in docker
#pwd: same as before
#docker pull immcantation/suite:4.3.0
#docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.3.0 bash 

for dir in ./by_vgene_project/OE_heavy_IGHV*; do 
cd $dir 
Rscript ../observed_expected_chargegainloss_heavy.R 
cd ../..
done 

for dir in ./by_vgene_project/OE_lambda_IGLV*; do
cd $dir 
Rscript ../observed_expected_chargegainloss_lambda.R
cd ../..
done

for dir in ./by_vgene_project/OE_kappa_IGKV*; do
cd $dir
Rscript ../observed_expected_chargegainloss_kappa.R
cd ../..
done
