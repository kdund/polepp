#!/bin/tcsh
#
if ( ${#argv} <= 0 ) then
    echo "Usage log2confbelt <log-file>"
    exit 1
endif

set input      = $1
set tmpout     = "`mktemp /tmp/$1_output.XXXXXX`"
set output     = "confbelt.root"

grep CONFBELT ${input} | cut -d ":" -f 2 > $tmpout

if ( $status ) then
    echo "Log file does not contain CONFBELT: lines!"
    exit 1
endif

ascii2root $tmpout
mv ${tmpout}.root $output
echo "* Confidence belt root file in $output ."
echo "* Run root plotconst.C to view the confidence belt ."
