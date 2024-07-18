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
savename = 'voxim19';

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

% Define number of sampling lines
np_cartesian = 256; % # of Cartesian lines
np_radial = 200; % # of projections/spokes
np_spiral = 48; % # of spiral arms
sampmode = 'demo'; % sampling mode: 'demo', 'simple', 'eachEcho' 
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
nc = 20; % # of coils

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
    np = np_cartesian;
elseif strcmp(samptraj,'radial')
    np = np_radial;
elseif strcmp(samptraj,'spiral')
    np = np_spiral;
end
% opts = prepareNUFFT3d(mtx,npar,np,slthick,samptraj,'goldenAngle_sorted_180','fov',fov,'FWshift',FWshift);
opts = prepareNUFFT(mtx,np,samptraj,'goldenAngle_sorted_180','fov',fov,'FWshift',FWshift);
% load(coilmap);
% cmap = gencmap([mtx mtx nnp, traj_sort = 'linear_sorted';par],nc,origcmap);
cmap = gencmap([mtx mtx npar],nc);
% cmap = permute(cmap, [1,2,4,3]);
% Sequence parameters
[seqparam,defseq] = setseqparam(sigtype,[np npar nset],sampmode);

% Bloch simulation
sigevo = gensigevo(tissueprop,seqparam);

% Convert voxels to phantom images (This step will take mins to hours, depending on the number of time frames)
%%
nt = length(defseq.demosig);
[nr,~] = size(opts.kx);
mixsamp = zeros(nr,np,nc,npar,nt,'single');
refimg = zeros(mtx, mtx, npar, nt, 2, 'single');
for itp = 1:nt
    imPall = model2voximg(phanimg(:,:,:,mod(defseq.demosig(itp)-1,tframe)+1),sigevo(defseq.demosig(itp),:,:)); % Ground truth images
    if itp == 1
        nval = calcnoiselvl(imPall,cmap);
    end
    mixsamp(:,:,:,:,itp) = voximg2ksp(imPall,cmap,nval,opts); % k-space + noise
    refimg(:,:,:,itp, :) = imPall;
end

save([savename '_mixsamp.mat'],'mixsamp','-v7.3')
save([savename '_refimg.mat'],'refimg','-v7.3')
fprintf('Data acquisition done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert 4D phantom k-space to images
reconimg = ksp2img(mixsamp,opts,cmap);
save([savename '_reconimg.mat'],'reconimg','-v7.3')
save([savename '_cmap.mat'],'cmap','-v7.3')
fprintf('Data reconstruction done\n');

%% Show phantom images
% showimg(reconimg(:,:,round(npar/2),:));colormap(gray);title('Reconstructed phantom images: axial plane')
% showimg(imrotate(squeeze(reconimg(round(mtx/2),:,:,:)),90));colormap(gray);title('Reconstructed phantom images: frontal plane')
figure
sli = 10;
subplot(1,2,1)
imshow(abs(reconimg(:,:,sli,1)),[]);
colorbar;
subplot(1,2,2)
imshow(abs(refimg(:,:,sli,1,1)+refimg(:,:,sli,1,2)),[]);
colorbar;

%% phantom check
addpath("utils/");
np2 = 64;
tic
op2 = prepareNUFFT3d(mtx,np2,np,slthick,samptraj,'goldenAngle_sorted_180','fov',fov,'FWshift',FWshift);
toc
imp = phantom(256,256, np2);
imp = repmat(imp, 1,1,np2);
kspsize = true([512, 200, np2]);
imsize = true([256, 256, np2]);
ksp = embed(op2.G * imp, kspsize);
im2 = embed(op2.G' * (op2.wib .* ksp(:)), imsize);
%%
figure
sli = 10;
subplot(1,2,1)
imshow(abs(imp(:,:,sli)),[]);
colorbar;
subplot(1,2,2)
imshow(abs(im2(:,:,sli)),[]);
colorbar;