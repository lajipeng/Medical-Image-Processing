function write_MHD3D(fname, data, info)
%-- Author£ºWang, Peng. Fudan University 
%-- Date£º2019/8/13
%-- Fuction£ºSave data and info in raw and mhd
%-- fname:where to save the mhd
%-- data: the image voxel data
%-- info£ºthe mhd info of image 
%%
% Transpose
x = permute(data,[2 1 3]);
% Write .raw
f = fopen(strcat([fname,'.raw']), 'wb');
fwrite(f, x(:), 'short');
fclose(f);
%%
% Write .mhd
fid=fopen(strcat([fname,'.mhd']),'w');
str = sprintf('ObjectType = %s',info.ObjectType{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('NDims = %d',info.NDims);
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('BinaryData = %s',info.BinaryData{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('BinaryDataByteOrderMSB = %s',info.BinaryDataByteOrderMSB{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('CompressedData = %s',info.CompressedData{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('TransformMatrix = %d %d %d %d %d %d %d %d %d',info.TransformMatrix(1),info.TransformMatrix(2),info.TransformMatrix(3),info.TransformMatrix(4),info.TransformMatrix(5),info.TransformMatrix(6),info.TransformMatrix(7),info.TransformMatrix(8),info.TransformMatrix(9));
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('Offset = %f %f %f' ,info.Offset(1),info.Offset(2),info.Offset(3));
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('CenterOfRotation = %d %d %d',info.CenterOfRotation(1),info.CenterOfRotation(2),info.CenterOfRotation(3));
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('AnatomicalOrientation = %s',info.AnatomicalOrientation{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('ElementSpacing = %f %f %f',info.ElementSpacing(1),info.ElementSpacing(2),info.ElementSpacing(3));
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('DimSize = %d %d %d',info.DimSize(1),info.DimSize(2),info.DimSize(3));
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('ElementType =  %s',info.ElementType{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
str = sprintf('ElementDataFile = %s',info.ElementDataFile{1});
fwrite(fid,str,'char');
fprintf(fid,'\n');
fclose(fid);