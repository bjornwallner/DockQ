#!/bin/bash
mkdir -p testdata_tmp
sed s/'diff'/'cp'/g run_test.sh  | sed s/'testdata'/'testdata_tmp'/g | bash

