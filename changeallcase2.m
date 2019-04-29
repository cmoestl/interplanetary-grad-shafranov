function changeallcase2(dirstring, writefile, polyorder,control,control2)


fid=fopen(strcat('ALLCASE.par'), 'wt');
 fprintf(fid, '%c',dirstring);
 fprintf(fid, '\n',dirstring);
 fprintf(fid, '%i\n',writefile);
 fprintf(fid, '%i\n',polyorder);
 fprintf(fid, '%i %i %i\n',control(1:3));
 fprintf(fid, '%i %f %f\n',control(4),control2(1:2));
fclose(fid);


