function [img] = img_contour(Image_dicom,dim)
%-- 作者：王鹏 复旦大学
%-- 时间：2019年
%-- 项目：肝癌疾病分类
%-- 读取DICOM图像轮廓
Image_dicom=im2bw(Image_dicom);
contour = bwperim(Image_dicom); 
contour = double(im2uint8(contour));
if dim == 3
    img(:,:,1)=contour;
    img(:,:,2)=0;
    img(:,:,3)=0;
else
    img=contour;
end
end


