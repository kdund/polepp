#!/bin/tcsh
#
if ( ${#argv} <= 0 ) then
    echo "Usage log2construct <log-file>"
    exit 1
endif

set input      = $1
set tmpout     = "`mktemp /tmp/$1_output.XXXXXX`"
set output     = "construct.dat"

echo "shyp/F:N/F:RL/F:P/F:norm/F" > $tmpout
echo "# Construct data created `date` by user <`whoami`>" >> $tmpout
echo "# shyp  : signal hypothesis" >> $tmpout
echo "# N     : number of observed events" >> $tmpout
echo "# RL    : likelihood ratio" >> $tmpout
echo "# P     : probability P(N|shyp)" >> $tmpout
echo "# norm  : normalisation" >> $tmpout

grep CONSTRUCT ${input} | cut -d ":" -f 2 >> $tmpout

if ( $status ) then
    echo "Log file does not contain CONSTRUCT: lines!"
    rm -f $tmpout
    exit 1
endif

mv ${tmpout} $output
echo "* Construct data file in $output ."
echo "* Run root plotconst.C to view the construct ."
