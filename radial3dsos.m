function opts = radial3dsos(N, Nz ,np,trajectory,viewOrder,varargin)
    opts = struct();
    opts.trajectory = trajectory;
    opts.viewOrder = viewOrder;
    opts.interleaves = [];
    opts.N = N;
    opts.sltk = 3;
    opts.fov = 420; % default field-of-view is 300mm2
    opts.readShift = 0;
    opts.phaseShift = 0;
    opts.correctionFactor =1;
    opts.G = [];
    opts.wib = [];
    opts.kx = [];
    opts.ky = [];
    opts.gridReadLow = []; % smallest readout point to grid and reconstruct (for spiral)
    opts.gridReadHigh = []; % lst readout point to grid and reconstruct (for spiral)
    opts.FWshift = 0;
    dcf_path = [];

    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'trajname'
                trajname = varargin{i+1};
            case 'fov' 
                opts.fov = varargin{i+1};
            case 'correctionFactor'
                opts.correctionFactor = varargin{i+1};
            case 'readShift'
                opts.readShift = varargin{i+1};
            case 'phaseShift'
                opts.phaseShift = varargin{i+1};
            case 'gridLow'
                opts.gridReadLow = varargin{i+1};
            case 'gridHigh'
                opts.gridReadHigh = varargin{i+1};
            case 'gridMiddle'
                opts.gridReadMiddle = varargin{i+1};
            case 'interleaves'
                opts.interleaves = varargin{i+1};
            case 'kspall'
                kspall = varargin{i+1};
            case 'FWshift'
                opts.FWshift = varargin{i+1};
            otherwise % skip it
        end
    end
    % opts.fovs = [opts.fov, opts.fov, opts.sltk*Nz];
    opts.ig = image_geom('nx', N, 'nz', Nz, ...
		'offsets', 'dsp', ... % (-n/2:n/2-1) for mri
		'fov', [opts.fov, opts.fov, opts.sltk*Nz]); % 20 cm transaxial FOV
    opts.dim = [N,N,Nz];

    switch opts.trajectory
        case 'radial'; prepare_radial();
        % case 'spiral'; prepare_spiral();
        % case 'cartesian'; prepare_cartesian();
        otherwise; error('Trajectory not supported')
    end
    function prepare_radial()
        
        golden_ratio = (sqrt(5)+1)/2;
        golden_angle = 180/golden_ratio;
        
        switch opts.viewOrder
            
            case 'linear_sorted'
                ang = 0:180/np:180-180/np;
            case 'linear_360'
                ang = 0:360/np:360-360/np;
            case 'linear_GA'
                ang = 0:golden_angle:golden_angle*(np-1);
                ang = rem(ang,180);
                [~, angix] = sort(ang);
                opts.angix = angix;
                ang = 0:180/np:180-180/np;
    %                 ang = ang(angix);
            case 'goldenAngle_sorted_180'
                ang = 0:golden_angle:golden_angle*(np-1);
                ang = rem(ang,180);
                ang = sort(ang);
                [~, angix] = sort(ang);
                opts.angix = angix;
            case 'goldenAngle_sorted_360'
                ang = 0:golden_angle:golden_angle*(np-1);
                ang = rem(ang,360);
                ang = sort(ang);
                [~, angix] = sort(ang);
                opts.angix = angix;
            case 'goldenAngle_180'
                ang = 0:golden_angle:golden_angle*(np-1);
                ang = rem(ang,180);
                [~, angix] = sort(ang);
                opts.angix = angix;
            case 'goldenAngle_360'
                ang = 0:golden_angle:golden_angle*(np-1);
                ang = rem(ang,360);
                [~, angix] = sort(ang);
                opts.angix = angix;
            case 'interleaved'
                temp = 0:180/np:180-180/np;
                af = opts.interleaves;
                assert(~isempty(af),'Must specify number of interleaves');
                assert(af == round(af),'Number of projections is not divisible by number of interleaves');
                ang = zeros(np/af,af);
                for u=1:af
                    ang(:,u) = temp(u:af:end);
                end
                ang = ang(:).';
                
            otherwise
                error('View ordering scheme not supported')
        end
        
        li = -N/2 : 0.5 : N/2-0.5;
        ky = li'*sind(ang);
        kx = li'*cosd(ang);
        
        omega = mri_trajectory_radial(opts.dim(1:2), opts.ig.fovs(1:2), ang);
        % convert to physical units
        omega3d = mri_trajectory_stack(omega, Nz);
        kspace = zeros(size(omega3d), 'single');
        for id = 1:length(opts.dim)
            dx = opts.ig.fovs(id) / opts.dim(id);
            kspace(:,id) = omega3d(:,id) / (2*pi) / dx;
        end
        opts.kspace = kspace;
        if 1 % plot trajectory
            omega3d = kspace .* (2*pi * opts.ig.deltas);
            im subplot
            plot3(omega3d(:,1), omega3d(:,2), omega3d(:,3), '.')
            titlef('%s with %d samples', trajectory, size(omega3d,1))
            axis_pipi, axis square
        end
        opts.mask3d = true(N,N,Nz);
        opts.maskSize = size(opts.mask3d);
        
        nufft_args = {opts.maskSize, [6,6,6], 2*opts.maskSize, opts.maskSize/2, 'table', 2^12, 'minmax:kb'};
        opts.G = Gmri(opts.kspace, opts.mask3d, 'fov', opts.ig.fovs, 'basis', {'rect'}, 'nufft', nufft_args);
        % wi = mri_density_comp(kspace,'voronoi','fix_edge',0,'G',G.Gnufft);
        
    end

end

function omega = mri_trajectory_stack(omega2, N3)
if mod(N3,2) == 0 % iseven
	o3 = ((0:(N3-1))/N3 - 0.5);
else % isodd
	o3 = (-(N3-1)/2 : (N3-1)/2) / N3;
end
o3 = single(o3 * 2 * pi);
o3 = repmat(o3, nrow(omega2), 1); % [N12,N3]
omega = repmat(omega2, N3, 1); % [N12*N3,2]
omega = [omega o3(:)]; % [N12*N3,3]
end


function omega = mri_trajectory_radial(dimxy, fovs, angles);
na_nr = 2*pi;	% default ensures proper sampling at edge of k-space
na = [];		% angular spokes (default: na_nr * nr)
% nr = max(opts.dimxy)/2;	% radial samples per (outward) spoke
nr = max(dimxy) * 2; % 512
ir = [];		% default: -nr/4:nr/4
omax = pi;		% maximum omega
if isempty(ir), ir = -nr/2:nr/2-1; end
if isempty(na), na = 4*ceil(na_nr * nr/4); end % mult of 4
om = ir/nr * (2*pi);
% ang = (0:na-1)/na * 2*pi;
[om, angles] = ndgrid(om, angles); % [nr+1, na]
omega = [col(om.*cosd(angles)) col(om.*sind(angles))];
end
