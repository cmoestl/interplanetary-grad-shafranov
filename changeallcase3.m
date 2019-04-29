function changeallcase3(dirstring, writefile, polyorder,control,control2)


fid=fopen([dirstring,'\ALLCASE.par'], 'wt');
 fprintf(fid, '%c',dirstring);
 fprintf(fid, '\n',dirstring);
 fprintf(fid, '%i\n',writefile);
 fprintf(fid, '%i\n',polyorder);
 fprintf(fid, '%i %i %i\n',control);
 fprintf(fid, '%i %f %f\n',control2);
fclose(fid);


