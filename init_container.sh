#!/bin/bash

chmod +x /home/dnanexus/ukbb_analysis*/run.sh
docker run -it -v /home/dnanexus:/home/dnanexus -w /home/dnanexus ensorp/ukbb_analysis /bin/bash

