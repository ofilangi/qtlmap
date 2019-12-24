#!/bin/bash
#qtlmap analyse 1 trait


function qtlmap_exec {
               echo "**************************************"
	       echo "DIR:`pwd`"
	       echo "ARGS= $*"
	       echo "**************************************"
	       echo
               qtlmap p_analyse_n --data-transcriptomic $*
	       cr=$?;if [ $cr -ne 0 ] ;then echo "## Error :--data-transcriptomic $* ##" ;exit $cr;fi
       }

options=$*

opt="--calcul=5 ${options}"
qtlmap_exec ${opt}

opt="--calcul=6 ${options}"
qtlmap_exec ${opt}





