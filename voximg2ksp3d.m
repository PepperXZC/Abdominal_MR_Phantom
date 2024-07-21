function mixsamp = voximg2ksp3d(imPall,cmap,nval,opts)

[~,~,npar,wfcomp,nframe] = size(imPall);
[~,~,~,nc] = size(cmap);
% [~,~,nc,~] = size(cmap);
im2kspmask = true(opts.nr,opts.np,npar);
    
mixsamp = zeros(opts.nr,opts.np,npar,nc,'single');

tstart = tic;
imW = repmat(imPall(:,:,:,1),[1 1 1 nc]).*cmap;
imF = repmat(imPall(:,:,:,2),[1 1 1 nc]).*cmap;

for c = 1:nc
    if strcmp(opts.trajectory,'cartesian')
        ksp2DW = squeeze(fft2n(imW(:,:,c,:),1,2)); % one slice/partition
        ksp2DF = squeeze(fft2n(imF(:,:,c,:),1,2));
    else
        ksp2DW(:,:,:,c) = embed(opts.G*imW(:,:,:, c),im2kspmask);
        ksp2DF(:,:,:,c) = embed(opts.G*imF(:,:,:, c),im2kspmask);
    end
    
    % Apply phase accumulation
    mixksp = ksp2DW(:,:,:,c)+ksp2DF(:,:,:,opts.c).*repmat(exp(1i*2*pi*opts.FWshift*30*10^-6*(0:opts.nr-1))',[1 opts.np npar]);
    
    % Add noise
    mixsamp(:,:,:,c) = mixksp+nval*(randn(opts.nr,opts.np,npar)+1i*randn(opts.nr,opts.np,npar));
end

disp(['Simulate current frame in ' num2str(toc(tstart)) ' sec']);
