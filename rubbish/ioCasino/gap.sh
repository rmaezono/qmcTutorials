#!/bin/sh

tegs=`grep 'Total energy ' 3_gs_dmc/out | grep + | awk '{print $4}'`
ebgs=`grep 'Total energy ' 3_gs_dmc/out | grep + | awk '{print $6}'`
tees=`grep 'Total energy ' 5_es_dmc/out | grep + | awk '{print $4}'`
ebes=`grep 'Total energy ' 5_es_dmc/out | grep + | awk '{print $6}'`

gap=`echo "($tees - $tegs)*27.21" | bc | awk '{printf "%2.1f", $1}'`
err=`echo "27.21*sqrt($ebgs^2 + $ebes^2)" | bc | awk '{printf "%2.1f", $1}'`
echo GAP = $gap +/- $err eV