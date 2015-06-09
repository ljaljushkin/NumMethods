#!/usr/bin/env bash

rm -f bin/experiments

for i in `seq 0 10`;
do
    cp ../test_data/input${i} bin/input
    cd bin; ./solution input output time

    echo "Time #${i} = $(<time)"
    cat time >> experiments
    echo >> experiments
    cd ../
done

#echo "Time results: "
#cat bin/experiments
