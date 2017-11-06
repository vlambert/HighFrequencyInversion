#!/bin/bash
 matlab -nodisplay -nosplash -nodesktop -r "try, run('MultiArray_HFI_Optimize.m'), catch, exit, end, exit;"
echo "matlab exit code: $?"
