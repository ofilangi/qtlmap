#!/bin/bash
#qtlmap analyse 1 trait


function qtlmap_exec {
               qtlmap p_analyse_n $1
	       cr=$?;if [ $cr -ne 0 ] ;then echo "## Error :$1 ##" ;exit $cr;fi
       }

options=$*

opt="--calcul=5 ${options}"
qtlmap_exec ${opt}

opt="--calcul=6 ${options}"
qtlmap_exec ${opt}



