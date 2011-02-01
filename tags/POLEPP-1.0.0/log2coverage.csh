#!/bin/tcsh
#
if ( ${#argv} <= 0 ) then
    echo "Usage log2coverage <log-file>"
    exit 1
endif

set input      = $1
set tmpout     = "`mktemp /tmp/$1_output.XXXXXX`"
set output     = "coverage.dat"

echo "shyp/F:effmean/F:effsigma/F:bkgmean/F:bkgsigma/F:cov/F:coverr/F:nlok/I:nltot/I:time/F" > $tmpout
echo "# Coverage data created `date` by user <`whoami`>" >> $tmpout

grep COVERAGE ${input} | cut -d ":" -f 2 >> $tmpout

if ( $status ) then
    echo "Log file does not contain COVERAGE: lines!"
    rm -f $tmpout
    exit 1
endif

mv ${tmpout} $output
echo "* Coverage data file in $output ."
echo "* Run root plotcov.C to view the coverage ."
