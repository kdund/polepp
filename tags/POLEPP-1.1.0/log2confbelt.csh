#!/bin/tcsh
#
if ( ${#argv} <= 0 ) then
    echo "Usage log2confbelt <log-file>"
    exit 1
endif

set input      = $1
set tmpout     = "`mktemp /tmp/$1_output.XXXXXX`"
set output     = "confbelt.dat"

echo "shyp/F:N1/F:N2/F:P/F" > $tmpout
echo "# Confidence belt data created `date` by user <`whoami`>" >> $tmpout
echo "# shyp  : signal hypothesis" >> $tmpout
echo "# N1,N2 : range in N for the given shyp" >> $tmpout
echo "# P     : probability for range" >> $tmpout
grep CONFBELT ${input} | cut -d ":" -f 2 >> $tmpout

if ( $status ) then
    echo "Log file does not contain CONFBELT: lines!"
    rm -f $tmpout
    exit 1
endif


mv ${tmpout} $output
echo "* Confidence belt data file in $output ."
echo "* Run root plotconst.C to view the confidence belt ."
