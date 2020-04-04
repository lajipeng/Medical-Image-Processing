function I = ImageRotate(filename,ang,isSameSize)
% Author：lajipeng
% Date: 2019/11/2
% Fuction: ImageRotate for image registration
rad = pi/180*ang;
OldImage = imread(filename);
figure(1),
imshow(OldImage);
[nrows,ncols] = size(OldImage);
OldWidth = nrows;
OldHeight = ncols;
if(isSameSize==0)
    %计算原图像四个角的坐标，以图像中心为原点
    oldX1=-(OldWidth-1)/2;
    oldY1=(OldHeight-1)/2;
    oldX2=(OldWidth-1)/2;
    oldY2=(OldHeight-1)/2;
    oldX3=-(OldWidth-1)/2;
    oldY3=-(OldHeight-1)/2;
    oldX4=(OldWidth-1)/2;
    oldY4=-(OldHeight-1)/2;
    %计算新图像四个角的坐标，以图像中心为原点
    newX1=oldX1*cos(rad)+oldY1*sin(rad);
    newY1=-oldX1*sin(rad)+oldY1*cos(rad);
    newX2=oldX2*cos(rad)+oldY2*sin(rad);
    newY2=-oldX2*sin(rad)+oldY2*cos(rad);
    newX3=oldX3*cos(rad)+oldY3*sin(rad);
    newY3=-oldX3*sin(rad)+oldY3*cos(rad);
    newX4=oldX4*cos(rad)+oldY4*sin(rad);
    newY4=-oldX4*sin(rad)+oldY4*cos(rad);
    %计算旋转后的图像宽度和高度
    NewWidth = round(max(abs(newX4-newX1),abs(newX3-newX2))+0.5);
    NewHeight = round(max(abs(newY4-newY1),abs(newY3-newY2))+0.5);
    %旋转前中心坐标
    a = round((OldWidth-1)/2+0.5);
    b = round((OldHeight-1)/2+0.5);
    %旋转后中心坐标
    c = round((NewWidth-1)/2+0.5);
    d = round((NewHeight-1)/2+0.5);
else
    a = round((OldWidth-1)/2+0.5);
    b = round((OldHeight-1)/2+0.5);
    c=a;
    d=b;
    NewWidth=OldWidth;
    NewHeight=OldHeight;
end

t1 = [1,0,0;0,1,0;-a,-b,1];
t2 = [cos(rad),-sin(rad),0;sin(rad),cos(rad),0;0,0,1];
t3 = [1,0,0;0,1,0;c,d,1];
T = t1*t2*t3;
tform = maketform('affine',T);

tx = zeros(NewWidth,NewHeight);
ty = zeros(NewWidth,NewHeight);

for i = 1:NewWidth
    for j = 1:NewHeight
        tx(i,j) = i;
    end
end

for i = 1:NewWidth
    for j = 1:NewHeight
        ty(i,j) = j;
    end
end

[w,z] = tforminv(tform,tx,ty);
NewImage = uint8(zeros(NewWidth,NewHeight));

for i = 1:NewWidth
    for j = 1:NewHeight
        source_x = w(i,j);
        source_y = z(i,j);
        if(source_x>=OldWidth-1||source_y>=OldHeight-1||double(uint8(source_x))<=0||double(uint8(source_y))<=0)
            NewImage(i,j) = 0;
        else
            if(source_x/double(uint16(source_x))==1.0&&source_y/double(uint16(source_y))==1.0)
                NewImage(i,j) = OldImage(int16(source_x),int16(source_y));
            else
                a = double(uint8(source_x));
                b = double(uint8(source_y));
                x11 = double(OldImage(a,b));
                x12 = double(OldImage(a,b+1));
                x21 = double(OldImage(a+1,b));
                x22 = double(OldImage(a+1,b+1));
                NewImage(i,j) = uint8((b+1-source_y)*((source_x-a)*x21+(a+1-source_x)*x11)+(source_y-b)*((source_x-a)*x22+(a+1-source_x)*x12));
            end
        end
    end
end
I = NewImage; 
figure(2),
imshow(I);
end
