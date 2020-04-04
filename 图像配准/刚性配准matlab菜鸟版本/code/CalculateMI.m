function mi = CalculateMI(a,b)
% Author£ºlajipeng
% Date: 2019/11/3
% Fuction: Calculate Mutual Information of image a and image b
a = double(a);
b = double(b);
[M,N]= size(a);
hab = zeros(256,256);
if max(max(a))~=min(min(a))
    a = (a-min(min(a)))/(max(max(a)-min(min(a))));
else 
    a = zeros(M,N);
end
if max(max(b))~=min(min(b))
    b = (b-min(min(b)))/(max(max(b)-min(min(b))));
else 
    b = zeros(M,N);
end
a = double(int16(a*255))+1;
b = double(int16(b*255))+1;

for i = 1:M
    for j = 1:N
        index_x = a(i,j);
        index_y = b(i,j);
        hab(index_x,index_y) = hab(index_x,index_y)+1;
    end
end
% Compute entropy Hab
habsum = sum(sum(hab));
index = find(hab~=0);
pab = hab/habsum;
Hab = sum(sum(-pab(index).*log2(pab(index))));

% Compute entropy Ha 
pa = sum(pab'); % sum along a dim 
index = find(pa~=0);
Ha = sum(sum(-pa(index).*log2(pa(index))));
% Compute entropy Hb
pb = sum(pab);  % sum along b dim 
index = find(pb~=0);
Hb = sum(sum(-pb(index).*log2(pb(index))));

mi = Ha + Hb - Hab;
end



