function O_new=refine_grid(O_trans)
% Refine image transformation grid of 1D b-splines with use of spliting matrix
%
%     O_new = refine_grid(O_trans)
%
% Function is written by D.Kroon University of Twente (September 2008)


Msplit=(1/8)*[4 4 0 0;1 6 1 0; 0 4 4 0;0 1 6 1; 0 0 4 4];

if(ndims(O_trans)==3)
    Olx=((size(O_trans,1)-2)*2-1)+2;
    O_newA=zeros(Olx,size(O_trans,2),2);
    for j=1:size(O_trans,2)
        for i=1:size(O_trans,1)-3
            for h=1:2
                P0=O_trans(i+0,j,h); P1=O_trans(i+1,j,h); P2=O_trans(i+2,j,h); P3=O_trans(i+3,j,h);
                Pnew=Msplit*[P0 P1 P2 P3]';
                O_newA(1+((i-1)*2),j,h)=Pnew(1);
                O_newA(2+((i-1)*2),j,h)=Pnew(2);
                O_newA(3+((i-1)*2),j,h)=Pnew(3);
                O_newA(4+((i-1)*2),j,h)=Pnew(4);
                O_newA(5+((i-1)*2),j,h)=Pnew(5);
            end
        end
    end
    
    Oly=((size(O_newA,2)-2)*2-1)+2;
    O_newB=zeros(size(O_newA,2),Oly,2);

    for j=1:size(O_newA,1)
        for i=1:size(O_newA,2)-3
            for h=1:2
                P0=O_newA(j,i+0,h); P1=O_newA(j,i+1,h); P2=O_newA(j,i+2,h); P3=O_newA(j,i+3,h);
                Pnew=Msplit*[P0 P1 P2 P3]';
                O_newB(j,1+((i-1)*2),h)=Pnew(1);
                O_newB(j,2+((i-1)*2),h)=Pnew(2);
                O_newB(j,3+((i-1)*2),h)=Pnew(3);
                O_newB(j,4+((i-1)*2),h)=Pnew(4);
                O_newB(j,5+((i-1)*2),h)=Pnew(5);
            end
        end
    end
    O_new=O_newB;
else
    
    Olx=((size(O_trans,1)-2)*2-1)+2;
    O_newA=zeros(Olx,size(O_trans,2),size(O_trans,3),3);
    for k=1:size(O_trans,3)
        for j=1:size(O_trans,2)
            for i=1:size(O_trans,1)-3
                for h=1:3
                    P0=O_trans(i+0,j,k,h); P1=O_trans(i+1,j,k,h); P2=O_trans(i+2,j,k,h); P3=O_trans(i+3,j,k,h);
                    Pnew=Msplit*[P0 P1 P2 P3]';
                    O_newA(1+((i-1)*2),j,k,h)=Pnew(1);
                    O_newA(2+((i-1)*2),j,k,h)=Pnew(2);
                    O_newA(3+((i-1)*2),j,k,h)=Pnew(3);
                    O_newA(4+((i-1)*2),j,k,h)=Pnew(4);
                    O_newA(5+((i-1)*2),j,k,h)=Pnew(5);
                end
            end
        end
    end
    
    Oly=((size(O_newA,2)-2)*2-1)+2;
    O_newB=zeros(size(O_newA,1),Oly,size(O_newA,3),3);

    for k=1:size(O_newA,3)
        for j=1:size(O_newA,1)
            for i=1:size(O_newA,2)-3
                for h=1:3
                    P0=O_newA(j,i+0,k,h); P1=O_newA(j,i+1,k,h); P2=O_newA(j,i+2,k,h); P3=O_newA(j,i+3,k,h);
                    Pnew=Msplit*[P0 P1 P2 P3]';
                    O_newB(j,1+((i-1)*2),k,h)=Pnew(1);
                    O_newB(j,2+((i-1)*2),k,h)=Pnew(2);
                    O_newB(j,3+((i-1)*2),k,h)=Pnew(3);
                    O_newB(j,4+((i-1)*2),k,h)=Pnew(4);
                    O_newB(j,5+((i-1)*2),k,h)=Pnew(5);
                end
            end
        end
    end
    
    Olz=((size(O_newB,3)-2)*2-1)+2;
    O_newC=zeros(size(O_newB,1),size(O_newB,2),Olz,3);

    for k=1:size(O_newB,2)
        for j=1:size(O_newB,1)
            for i=1:size(O_newB,3)-3
                for h=1:3
                    P0=O_newB(j,k,i+0,h); P1=O_newB(j,k,i+1,h); P2=O_newB(j,k,i+2,h); P3=O_newB(j,k,i+3,h);
                    Pnew=Msplit*[P0 P1 P2 P3]';
                    O_newC(j,k,1+((i-1)*2),h)=Pnew(1);
                    O_newC(j,k,2+((i-1)*2),h)=Pnew(2);
                    O_newC(j,k,3+((i-1)*2),h)=Pnew(3);
                    O_newC(j,k,4+((i-1)*2),h)=Pnew(4);
                    O_newC(j,k,5+((i-1)*2),h)=Pnew(5);
                end
            end
        end
    end
    O_new=O_newC;
end

