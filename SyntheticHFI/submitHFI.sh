#!/bin/bash
 matlab -nodisplay -nosplash -nodesktop -r "try, run('MultiFreq_Complex.m'), catch err, fid=fopen('errorFile','w+'),fprintf(fid,'%s',err.getReport('extended','hyperlinks','off')),fclose(fid), exit, end, exit;"
echo "matlab exit code: $?"
