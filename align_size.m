function img2 = align_size(img1,Sy2,Sx2,Sz2, mode, padValue)

if(nargin <= 5)
    padValue = 0; 
end

if(nargin<=4)
    mode = 1; % center pad
end

[Sy1,Sx1,Sz1] = size(img1);

Sx = max(Sx1,Sx2);
Sy = max(Sy1,Sy2);
Sz = max(Sz1,Sz2);
img2 = ones(Sy,Sx,Sz)*padValue;

if mode == 1
    imgTemp = single(ones(Sy,Sx,Sz)*padValue);
    Sox = round((Sx-Sx1)/2)+1;
    Soy = round((Sy-Sy1)/2)+1;
    Soz = round((Sz-Sz1)/2)+1;
 %   img2(Soy:Soy+Sy1-1,Sox:Sox+Sx1-1,Soz:Soz+Sz1-1) = img1;
    imgTemp(Soy:Soy+Sy1-1,Sox:Sox+Sx1-1,Soz:Soz+Sz1-1) = img1;
     Sox = round((Sx-Sx2)/2)+1;
     Soy = round((Sy-Sy2)/2)+1;
     Soz = round((Sz-Sz2)/2)+1;
     img2 = imgTemp(Soy:Soy+Sy2-1,Sox:Sox+Sx2-1,Soz:Soz+Sz2-1);
else   % left top side pad
    img2 = img2 + padValue;
    img2(1:Sy1,1:Sx1,1:Sz1) = img1;
end

end


