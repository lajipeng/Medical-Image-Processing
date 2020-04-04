function [img] = windw_adjust(img,info)
%-- 作者：王鹏 复旦大学
%-- 时间：2019年
%-- 项目：肝癌疾病分类
%-- 调整DICOM图像窗宽、窗位
%-- img和info分别为用dicomread和dicominfo函数获取到的图像信息
WindowCenter = info.WindowCenter(1)/info.RescaleSlope - info.RescaleIntercept; %info.WindowCenter=40
WindowWidth = info.WindowWidth(1)/info.RescaleSlope;
I = mat2gray(img,[WindowCenter-(WindowWidth/2),WindowCenter+(WindowWidth/2)]); % window
img = double(im2uint8(I));
%fd.Hu = img/255*info.WindowWidth(1)+WindowCenter(1)-(WindowWidth(1)/2)+info.RescaleIntercept; % CT intensity (Hu)
end