%% main_phantom_generation
%
% Preparation for Abdominal_MR_Phantom generation:
% (1) user-defined file name
% (2) sequence type
% (3) sampling trajectory
% (4) FOV and spatial resolution
% (5) sampling lines
% (6) respratory motion file and temporal resolution
% (7) coil sensitivity file and coil numbers
% (8) Fat-water chemical shift frequency
%
%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define tissue properties and respiratory motion.
% Generate voxelized volumetric mask with indexes.
% ------------------------------------------------
%
%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sequence parameters.
% Generate signal evolution using Bloch simulation.
% Generate k-space data for each time frame.
% ------------------------------------------------
%
%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data reconstruction and gridding
% ------------------------------------------------
%
% -----------------------------------------------------------------------------------------
% Realistic 4D abdominal phantom for magnetic resonance imaging
% Wei-Ching (Tina) Lo
% wxl317@case.edu
% Case Western Reserve University
% July 2018
% -----------------------------------------------------------------------------------------
 
clear; close all;

% Define file name
% savename = 'voxim15';
savename = 'voxim23';

% Define sequence parameters:
% 'SpoiledGradientEcho' => T1-weighted image
% 'SpoiledGradientEchoWithFatSat' => Perfusion-weighted image
% 'InversionRecoveryLookLocker' => T1 mapping
% 'SingleSpinEcho' => T2 mapping
% 'SingleSpinEchoWithFatSat' => Diffusion-weighted image
% 'MultiEchoSpoiledGradientEcho' => Proton density fat fraction
sigtype = 'SpoiledGradientEcho';

% Define sampling trajectory:
% 'cartesian'
% 'radial'
% 'spiral'
% samptraj = 'cartesian';
samptraj = 'radial';

% Define 3D resolution and FOV
fov = 420; % field-of-view (mm)
mtx = 256; % matrix size
npar = 64; % # of partitions
slthick = 3; % slice thickness (mm)
nset = 1; % # of sets
type = 2; % 2 = 2D, 3 = 3D

% Define number of sampling lines
np_cartesian = 256; % # of Cartesian lines
np_radial = 21; % # of projections/spokes
np_spiral = 48; % # of spiral arms
sampmode = 'simple'; % sampling mode: 'demo', 'simple', 'eachEcho' 
% Note that "sampmode" will affect the simulation time (mins to hours, depending on the number of Echoes).
% Use 'eachEcho' option only when Echo-by-Echo simulation is nessasary

% Define 3D respiratory motion curve and temporal resolution (eq. each TR or each volume)
respmotion = 'respmov.mat';
tempres = 400; % temporal resolution (ms)
tempdur = 4000; % duration of each respiratory cycle (ms)
SImov = 13; % largest superior-inferior (SI) excursion (mm)
APmov = 6.5; % largest anterior-posterior (AP) excursion (mm)
LRmov = 2; % largest left-right (LR) excursion (mm)

% Define coil sensitivity map
coilmap = 'origcmap.mat';
nc = 36; % # of coils

% Define fat-water chemical shift
FWshift = 220; % water and fat are separated by approximately 440Hz in a 3T static field

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath = fileparts(mfilename('fullpath'));
cd(datapath)
addpath(genpath(datapath));

% Load defined tissue properties: T1, T2, ADC, PDFF
tissueprop = tissueproperty;

% respiratory motion pattern
load(respmotion); % respmov
tframe = floor(tempdur/tempres);
linv = genresp(respmov,tframe,[SImov APmov LRmov]);

% Convert mesh models to voxels (This step will take mins to hours, depending on the number of time frames)
xpts = -round((fov/420)*mtx/2)+1:round((fov/420)*mtx/2);
ypts = -round((fov/420)*mtx/2)+1:round((fov/420)*mtx/2); % x matrix size = y matrix size
zpts = linspace(round(-npar*slthick/3+35),round(npar*slthick/3+35),npar);
phanimg = mesh2model(tissueprop,linv,xpts,ypts,zpts);

% Show phantom images
showimg(phanimg);title('Indexed phantom images: axial plane');

save([savename '_phanimg.mat'],'phanimg','-v7.3')
fprintf('Phantom generation done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare FFT/NUFFT
if strcmp(samptraj,'cartesian')
    opts.np = np_cartesian;
elseif strcmp(samptraj,'radial')
    opts.np = np_radial;
elseif strcmp(samptraj,'spiral')
    opts.np = np_spiral;
end
np = opts.np;

%% Simulation

% Sequence parameters
[seqparam,defseq] = setseqparam(sigtype,[np npar nset],sampmode);

% Bloch simulation
sigevo = gensigevo(tissueprop,seqparam);
%% coil sensitivity maps
% cmap = gencmap([mtx mtx npar],nc);
cmap = mri_sensemap_sim('chat', 0, 'nx', mtx, ...
	'rcoil', [], 'ncoil', nc, 'coil_distance', 1.2);
rss = sqrt(sum(cmap.*conj(cmap), 3));
cnorm = cmap ./ rss;
cnorm(isnan(cnorm)) = 0;
cmap = cnorm;
if 0
    figure;
    imshow(sum(cmap.*conj(cmap), 3),[]);
    colorbar;
end
save([savename, '_cmap.mat'],'cmap','-v7.3')

% load(coilmap);    
% cmap = gencmap([mtx mtx nnp, traj_sort = 'linear_sorted';par],nc,origcmap);
%% Convert voxels to phantom images (This step will take mins to hours, depending on the number of time frames)
nt = length(defseq.demosig);
nphase = size(phanimg, 4);
nr = 2*mtx;
opts.nr = nr;
% refimg = zeros(mtx, mtx, npar, nt, 2);
refimg = zeros(mtx, mtx, npar, nt, nphase, 2);
for iphase = 1:nphase
    for itp = 1:nt
        % imPall = model2voximg(phanimg(:,:,:,mod(defseq.demosig(itp)-1,tframe)+1),sigevo(defseq.demosig(itp),:,:)); % Ground truth images
        imPall = model2voximg(phanimg(:,:,:,iphase),sigevo(defseq.demosig(itp),:,:)); % Ground truth images
        if itp == 1
            nval = calcnoiselvl(imPall, cmap);
        end
        refimg(:,:,:,itp, iphase, :) = imPall;
    end
end

save([savename '_refimg.mat'],'refimg','-v7.3')
fprintf('Data acquisition done\n');


%% ktraj, dcf, gmri prepare

if type == 2
    opts = prepareNUFFT(nt, mtx,npar,np,samptraj,'linear_GA','fov',fov,'FWshift',FWshift,'sltk',slthick);
else
    dcf_path = strcat(savename, '_dcf.mat');
    opts = radial3dsos(mtx,npar,np,samptraj, ...
        'linear_GA','FWshift',FWshift);
    if ~exist(dcf_path)
        tic
        wi = ir_mri_density_comp(opts.kspace,'pipe', 'G', opts.G.Gnufft, ...
            'arg_pipe', {'fov', [opts.ig.fovs], 'niter', 60});
        toc
        save(dcf_path,'wi','-v7.3')
    else
        wi = load(dcf_path);
    end
    opts.wib = wi.wi;
end
wi_max = 1.05 / prod([fov, fov]);
%% save opt info
save([savename '_opts.mat'],'opts','-v7.3')
%% plot traj
im plc 1 2
if 1
    nsample = size(opts.kspace{1},1);
    im subplot
    plot([1 nsample], wi_max * [1 1], 'm-', ...
        1:nsample, opts.wib{1}, '.');
    titlef('DCF')
    xlim([1 nsample]), xtick([1 nsample])
    legend('wi max', ['DCF ' 'voronoi'], ...
        'location', 'southeast'), drawnow
end
if 1
    dx = opts.ig.deltas;
	omega = opts.kspace{1} .* (2*pi .* dx(1:2));
	im subplot
	plot(omega(:,1), omega(:,2), '.')
	titlef('radial with %d samples', size(omega,1))
	axis_pipi, axis square
end

%% transform to kspace
if type == 2
    % mixsamp = zeros(nr,np,npar,nc,nt,'single'); % 3D
    mixsamp = zeros(nr,np,nc,npar,nt,'single'); % 2D
    for itp = 1:nt
        tic;
        mixsamp(:,:,:,:,itp) = voximg2ksp(squeeze(refimg(:,:,:,itp, :)),itp,cmap,nval,opts); % 2D  k-space + noise
        timeElapsed=toc;    fprintf('Time for recon one frame: %f seconds.\n', timeElapsed);
    end
else
    cmap_pmt = permute(cmap,[1,2,4,3]);
    mixsamp = zeros(nr,np,npar,nc,nt,'single'); % 2D
    for itp = 1:nt
        tic;
        mixsamp(:,:,:,:,itp) = voximg2ksp3d(squeeze(refimg(:,:,:,itp, :)),cmap_pmt,nval,opts); % 3D k-space + noise
        timeElapsed=toc;    fprintf('Time for recon one frame: %f seconds.\n', timeElapsed);
    end
end

save([savename '_mixsamp.mat'],'mixsamp','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert 4D phantom k-space to images
if type == 2
    tic;reconimg = ksp2img(mixsamp,opts,cmap);toc
else
    tic;reconimg = ksp2img3d(mixsamp,opts,cmap_pmt);toc
end
save([savename '_reconimg.mat'],'reconimg','-v7.3')
fprintf('Data reconstruction done\n');

%% Show phantom images
showimg(reconimg(:,:,round(npar/2),:));colormap(gray);title('Reconstructed phantom images: axial plane')
showimg(imrotate(squeeze(reconimg(round(mtx/2),:,:,:)),90));colormap(gray);title('Reconstructed phantom images: frontal plane')


