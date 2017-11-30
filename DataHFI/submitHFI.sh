#!/bin/bash
# matlab -nodisplay -nosplash -nodesktop -r "try, run('MultiArray_HFI_Optimize.m'), catch, exit, end, exit;"
 matlab -nodisplay -nosplash -nodesktop -r "try, run('MultiFreq_Optimize.m'), catch err, fid=fopen('errorFile','w+'),fprintf(fid,'%s',err.getReport('extended','hyperlinks','off')),fclose(fid), exit, end, exit;"
echo "matlab exit code: $?"
