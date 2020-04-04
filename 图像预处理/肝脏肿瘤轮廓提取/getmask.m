function [outputArg1,outputArg2] = getmask(fname,mhd_path1,raw_path1,mhd_path2,raw_path2)
%%
% Read the mhd and decide what the size of one picture and what type of the picture 

[Image_Series,info] = read_mhd_raw(mhd_path1,raw_path1);
%%
% Extract the contour of tumor and show in other phases
count = size(Image_Series);
img_contour = zeros(count);
se = strel('square',2);
contour = zeros(count);
figure(1),
for i = 1:1:count(3)
    img_contour(:,:,i) = img2bw(Image_Series(:,:,i));
    bw = imerode(img_contour(:,:,i),se);
    bw = imerode(bw,se);
    bw = imerode(bw,se);
    contour(:,:,i) = bwperim(bw,8);
    contour(:,:,i) = imfill(contour(:,:,i),'holes');
    contour(:,:,i) = bwareaopen(contour(:,:,i),50,8);
end
%%
[Dicom_Series,mhd_info] = read_mhd_raw(mhd_path2,raw_path2);

Result_Series = zeros(count);

Result_Series = Dicom_Series.*contour;

write_MHD3D(fname, Result_Series, mhd_info)

end


