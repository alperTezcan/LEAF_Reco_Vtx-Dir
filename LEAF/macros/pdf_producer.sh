#!/bin/bash +x

args=("$@")

/sps/t2k/atezcan/HK/LEAF/macros/ProducePDF_new -f ${args[0]} -o ${args[1]}
