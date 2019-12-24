#!/bin/bash
#qtlmap analyse 1 trait


function qtlmap_exec {
               echo "**************************************"
	       echo "DIR:`pwd`"
	       echo "ARGS= $*"
	       echo "**************************************"
	       echo
               qtlmap p_analyse --data-transcriptomic $*
	       cr=$?;if [ $cr -ne 0 ] ;then echo "## Error : --data-transcriptomic  $* ##" ;exit $cr;fi
       }

options=$*

opt="--calcul=1 ${options}"
qtlmap_exec ${opt}

opt="--calcul=1 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=3 ${options}"
qtlmap_exec ${opt}

opt="--calcul=3 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=4 ${options}"
qtlmap_exec ${opt}

opt="--calcul=4 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=21 ${options}"
qtlmap_exec ${opt}

opt="--calcul=21 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=25 ${options}"
qtlmap_exec ${opt}

opt="--calcul=25 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=26 ${options}"
qtlmap_exec ${opt}

opt="--calcul=26 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=27 ${options}"
qtlmap_exec ${opt}

opt="--calcul=27 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=28 ${options}"
qtlmap_exec ${opt}

opt="--calcul=28 --qtl=2 ${options}"
qtlmap_exec ${opt}



