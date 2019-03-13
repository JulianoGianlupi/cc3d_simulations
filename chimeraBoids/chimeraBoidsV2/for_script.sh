#!/bin/bash/

for i in {1..400};
do
    echo $i
    qsub submit_simpleBoids.sh
done