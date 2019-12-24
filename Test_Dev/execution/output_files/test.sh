#!/bin/bash
#qtlmap analyse 1 trait


function qtlmap_exec {
               qtlmap p_analyse $1
	       cr=$?;if [ $cr -ne 0 ] ;then echo "## Error :$1 ##" ;exit $cr;fi
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

opt="--calcul=3 --qtl=3 ${options}"
qtlmap_exec ${opt}

opt="--calcul=4 ${options}"
qtlmap_exec ${opt}

opt="--calcul=4 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=4 --qtl=3 ${options}"
qtlmap_exec ${opt}

opt="--calcul=21 ${options}"
qtlmap_exec ${opt}

opt="--calcul=21 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=25 ${options}"
qtlmap_exec ${opt}

opt="--calcul=25 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=25 --qtl=3 ${options}"
qtlmap_exec ${opt}

opt="--calcul=26 ${options}"
qtlmap_exec ${opt}

opt="--calcul=26 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=26 --qtl=3 ${options}"
qtlmap_exec ${opt}

opt="--calcul=27 ${options}"
qtlmap_exec ${opt}

opt="--calcul=27 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=27 --qtl=3 ${options}"
qtlmap_exec ${opt}
opt="--calcul=28 ${options}"
qtlmap_exec ${opt}

opt="--calcul=28 --qtl=2 ${options}"
qtlmap_exec ${opt}

opt="--calcul=28 --qtl=3 ${options}"
qtlmap_exec ${opt}



