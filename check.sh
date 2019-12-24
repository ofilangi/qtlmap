#!/bin/bash

function test_exit {
  $1
  cr=$?;if [ $cr -ne 0 ] ;then echo "## Error `pwd`:$1 ##" ;exit $cr;fi
}


echo " ** Execution ** "
pushd testdev/execution/ld
test_exit ./test.sh
popd

pushd testdev/execution/output_files
test_exit ./test.sh
test_exit ./test_multi.sh
#snp
test_exit ./test.sh "--snp"
test_exit ./test_multi.sh "--snp"
#simulation...
test_exit ./test.sh "--nsim=10"
test_exit ./test_multi.sh "--nsim=10"
test_exit ./test.sh "--nsim=10 --permute"
test_exit ./test_multi.sh "--nsim=10 --permute"
#confidence intervals
test_exit ./test.sh "--ci=1,2,3,4 --ci-nsim=10"
test_exit ./test_multi.sh "--ci=1,2,3,4 --ci-nsim=10"
popd

pushd testdev/execution/calcul_corcd
test_exit ./test.sh "--calcul-cd"
test_exit ./test_multi.sh "--calcul-cd"
test_exit ./test.sh "--calcul-cd --nsim=10"
test_exit ./test_multi.sh "--calcul-cd --nsim=10"
test_exit ./test.sh "--calcul-cd --nsim=10 --permute"
test_exit ./test_multi.sh "--calcul-cd --nsim=10 --permute"
test_exit ./test.sh "--calcul-cd --nsim=10 --snp"
test_exit ./test_multi.sh "--calcul-cd --nsim=10 --snp"
popd

pushd testdev/execution/transcript
test_exit ./test.sh
test_exit ./test_multi.sh
#snp
test_exit ./test.sh "--snp"
test_exit ./test_multi.sh "--snp"
#simulation...
test_exit ./test.sh "--nsim=10"
test_exit ./test_multi.sh "--nsim=10"
test_exit ./test.sh "--nsim=10 --permute"
test_exit ./test_multi.sh "--nsim=10 --permute"
#confidence intervals
test_exit ./test.sh "--ci=1,2,3,4 --ci-nsim=10"
test_exit ./test_multi.sh "--ci=1,2,3,4 --ci-nsim=10"
popd
pushd testdev/execution/ld
test_exit ./test.sh "--snp"
popd
