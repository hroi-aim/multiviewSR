function I_SIM = Shrink(data, phase_number, initial_phase)
[ny, nx, nz] = size(data);
I_shrink = zeros(2*ny, nx, nz);
lines = floor(ny/phase_number)-4;

hp = 1; % hlaf width of the pinhole
hw = 3; % slit width or set to 2 that matched to real experiments, but will be more fluctuation
sigma = (2*hp+1)/2.35;
for y=1:2*hw+1
     slit(y)  = exp(-1/2*((y-hw-1)/sigma)^2);
end

%slit = slit - min(slit(:));
SlitImage = (repmat(slit, [nx, 1]))';

for k=1:lines
    for pp=1:phase_number    
        index = initial_phase(pp);      
        c0 = round(index) + k*phase_number;
        SubImage = data(c0-hw:c0+hw,:,pp);
        SubImage = SubImage.*SlitImage;
        c1 = c0*2; 
        I_shrink(c1-hw:c1+hw,:,pp) = SubImage;
    end   
end

I_SIM = imresize(squeeze(sum(I_shrink,3)),[2*ny, 2*nx],'bilinear');

end