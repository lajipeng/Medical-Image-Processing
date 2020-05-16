clear all;
close all;
clc;

imagepath='E:\image\';
path_save='E:\manual\';
datafiles = dir(imagepath);
file_no = size(datafiles,1);
num=0;
for k=3:file_no
    
    filename = datafiles(k).name;
    datapath=[imagepath,filename];
    i=k-2;
    Img= dicomread(datapath);
%     segment by hand
    [BW3,xi,yi]= roipoly(Img);
    num=num+1;
    XY{num}=[xi,yi];
    BW(:,:,num)=BW3;
    IMG(:,:,num)=Img;
%     figure(1);
imshow(Img,[0,600]);
% hold on;contour(BW3,[0.1 0.1],'.w','linewidth',1);
    name=[path_save,num2str(i),'.tif'];
    saveas(figure(1),name);
pause(0.03)
end

save data_2 IMG XY BW


 