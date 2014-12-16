function [E,segCnt] = segmentDetect(I)

%Detect segments in image

%Takes as input an image after a Felzenswalb segmentation
I = rgb2gray(I);

[h,w] = size(I);

segCnt = 2;
E = ones(h,w);

for r=2:h
    if (I(r,1) == I(r-1,1))
        E(r,1) = E(r-1,1);
    else
        E(r,1) = segCnt;
        segCnt = segCnt + 1;
    end
end

for c=2:w
    if (I(1,c) == I(1,c-1))
        E(1,c) = E(1,c-1);
    else
        E(1,c) = segCnt;
        segCnt = segCnt + 1;
    end
end

for r=2:h
    for c=2:w
        if (I(r,c) == I(r,c-1))
            E(r,c) = E(r,c-1);
        elseif (I(r,c) == I(r-1,c))
            E(r,c) = E(r-1,c);
        elseif (I(r,c) == I(r-1,c-1))
            E(r,c) = E(r-1,c-1);
        else 
             E(r,c) = segCnt;
             segCnt = segCnt + 1;
        end
    end
end

segCnt = segCnt - 1;

end


