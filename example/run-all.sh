#!/bin/bash

for method in pc-grid-i pc-i fix-i
do 
    ./higher-order-example -m $method -i 10 --sigma 7 --eta 7 -g 1 medium-test
done
