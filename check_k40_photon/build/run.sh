#!/bin/bash
source /home/hufan/env_hailing.sh
make -j8
./main ../config/config_10.yaml
mv data.root 10/

./main ../config/config_20.yaml
mv data.root 20/