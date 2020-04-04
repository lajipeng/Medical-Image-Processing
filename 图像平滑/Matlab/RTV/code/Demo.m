% Demo script
% Uncomment each case to see the results
% I1 = imread('imgs/262.png');
I1 = imread('C:/Github/RTV_Smooth-master/#2_50_TVSmooth.png');
% imwrite(I,'imgs/261.jpg')
S_1 = tsmooth(I1,0.015,3.0,0.02,4);
% imwrite(S_1,'imgs/212.jpg')
figure, 
subplot(1,2,1)
imshow(I1);
subplot(1,2,2)
imshow(S_1);
% subplot(2,2,3)
% imshow(I)
% subplot(2,2,4)
% imshow(S_2)

% I = (imread('imgs/graffiti.jpg'));
% S = tsmooth(I,0.015,3);
% figure, imshow(I), figure, imshow(S);

% I = (imread('imgs/crossstitch.jpg'));
% S = tsmooth(I,0.015,3);
% figure, imshow(I), figure, imshow(S);


% I = (imread('imgs/mosaicfloor.jpg'));
% S = tsmooth(I, 0.01, 3, 0.02, 5);
% figure, imshow(I), figure, imshow(S);






