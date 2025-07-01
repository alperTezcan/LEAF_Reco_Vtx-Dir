#!/bin/bash +x

args=("$@")

/sps/t2k/atezcan/HK/LEAF/macros/ProducePDF_DIR -f ${args[0]} -o ${args[1]} -r ${args[2]}
