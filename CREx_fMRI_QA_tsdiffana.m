function CREx_tsdiffana
    % Author: Valérie Chanoine, Research Engineer at Brain and Language
    % Institute (http://www.blri.fr/)
    % From Matthew Brett website - 'Time series diagnostics' script (see
    % http://imaging.mrc-cbu.cam.ac.uk/imaging/DataDiagnostics)
    % 
    
%     w.dataDir           = 'F:\IRM\Intermod2\Data\LOCA';            % root directory
%     w.prefix            = 'LOCA';
%     w.sessions          = {'func01'};                             % session directory (parent=functional)

    
    w.dataDir           = 'F:\IRM\Intermod2\Data\AUDIO';            % root directory
    w.prefix            = 'AUDIO';                                  
    w.sessions          = {'func01', 'func02', 'func03', 'func04'}; % session directory (parent=functional)

    w.subjects          = {'01LA', '02CE'};                         % list of subjects 
    w.funcDir           =  'Functional';                            % functional directory (parent=subject)    
    w.dummy 			=  0;                                       % number of dummy files
    
    w.subjects          = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK',...
        '10TL', '11CL', '12SP', '13CC', '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'};  % subject directory (parent=dataDir)
      

    % Loop on subjects
    for iS=1:numel(w.subjects)  
 
        w.subName          = w.subjects{iS};
        w.subPath          = fullfile (w.dataDir,  w.subjects{iS}); 
        w.funcPath         = fullfile (w.subPath, w.funcDir);
        
        tsdiffana(w)
    end
    
end

function tsdiffana(w)

    % Loop for sessions
    for j=1:numel(w.sessions)
          
        % Get EPI Sliced files without dummy files
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName  '.*'  w.prefix '.*\.nii$']);        
        nScans = get_nii_frame(f);
        imgs = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName  '.*'  w.prefix '.*\.nii$'], w.dummy+1:nScans);        
    
        [td globals slicediff] = timediff(imgs,'mv');
        tdfn = fullfile(w.funcPath, w.sessions{j},'timediff.mat');
        save(tdfn, 'imgs', 'td', 'globals', 'slicediff');
        %%
        tsdiffplot(tdfn,[], '');
        saveas(gcf, fullfile(w.funcPath, w.sessions{j}, [w.prefix '_session_' num2str(j)]), 'png');
        %%
    end  
end

function [imdiff, g, slicediff] = timediff(imgs, flags)
    % Analyses slice by slice variance across time series
    % FORMAT [imdiff, g, slicediff] = timediff(imgs, flags)
    %
    % imgs   - string or cell or spm_vol list of images
    % flags  - specify options; if contains:
    %           m - create mean var image (vmean*), max slice var image
    %               (vsmax*) and scan to scan variance image (vscmean*)
    %           v - create variance image for between each time point
    %
    % imdiff - mean variance between each image in time series
    % g      - mean voxel signal intensity for each image
    % slicediff - slice by slice variance between each image
    %
    % Matthew Brett 17/7/00

    imgs = spm_vol(char(imgs));
    V1 = imgs(1);
    Vr = imgs(2:end);

    ndimgs = numel(imgs)-1;
    Hold = 0;

    if any(flags == 'v') % create variance images
        for i = 1:ndimgs
            vVr(i) = makevol(Vr(i),'v',16); % float
        end
    end
    if any(flags == 'm') % mean /max variance
        mVr = makevol(V1,'vmean',16);
        sVr = makevol(V1,'vscmean',16);
        xVr = makevol(V1,'vsmax',16);
    end

    [xydim zno] = deal(V1.dim(1:2),V1.dim(3));

    p1 = spm_read_vols(V1);
    slicediff = zeros(ndimgs,zno);
    g = zeros(ndimgs,1);
    for z = 1:zno % across slices
        M = spm_matrix([0 0 z]);
        pr = p1(:,:,z); % this slice from first volume
        if any(flags == 'm')
            [mv sx2 sx mxvs]  = deal(zeros(size(pr)));
        end
        % SVD is squared voxel difference (usually a slice of same)
        % MSVD is the mean of this measure across voxels (one value)
        % DTP is a difference time point (1:T-1)
        cmax = 0; % counter for which slice has the largest MSVD
        % note that Vr contains volumes 2:T (not the first)
        for i = 1:ndimgs % across DTPs
            c = spm_slice_vol(Vr(i),M,xydim,Hold); % get slice from this time point
            v = (c - pr).^2; % SVD from this slice to last
            slicediff(i,z) = mean(v(:)); % MSVD for this slice
            g(i) = g(i) + mean(c(:)); % simple mean of data
            if slicediff(i,z)>cmax  % if this slice has larger MSVD, keep
                mxvs = v;
                cmax = slicediff(i,z);
            end
            pr = c; % set current slice data as previous, for next iteration of loop
            if any(flags == 'v') % write individual SVD slice for DTP
                vVr(i) = spm_write_plane(vVr(i),v,z);
            end
            if any(flags == 'm')
                mv = mv + v; % sum up SVDs for mean SVD (across time points)
                sx = sx + c; % sum up data for simple variance calculation
                sx2 = sx2 + c.^2; % sum up squared data for simple variance
                % calculation
            end
        end
        if any(flags == 'm') % mean variance etc
            sVr = spm_write_plane(sVr,mv/(ndimgs),z); % write mean of SVDs
            % across time
            xVr = spm_write_plane(xVr,mxvs,z); % write maximum SVD
            mVr = spm_write_plane(mVr,(sx2-((sx.^2)/ndimgs))./(ndimgs-1),z);
            % (above) this is the one-pass simple variance formula
        end
    end

    g = [mean(p1(:)); g/zno];
    imdiff = mean(slicediff,2);
end

function Vo = makevol(Vi, prefix, datatype)
    Vo = Vi;
    fn = Vi.fname;
    [p f e] = fileparts(fn);
    Vo.fname = fullfile(p, [prefix f e]);
    switch lower(spm('ver'))
        case {'spm5','spm8','spm8b','spm12','spm12b'}
            Vo.dt = [datatype 0];
            Vo = spm_create_vol(Vo);
        otherwise
            error('What ees thees version "%s"', spm('ver'));
    end
end

function h = tsdiffplot(tdfn,fg,flags, varargin)
    % tsfiffplot - plots image difference etc info
    % FORMAT tsdiffplot(tdfn,fg,flags, varargin)
    % 
    % tdfn       - time difference file name - mat file with diff parameters
    % fg         - figure handle of figure to display in [spm graphics]
    % flags      - zero or more of 
    %              'r' - display realignment parameters
    %              These are either passed as filename in next arg or
    %              collected from the same directory as the above .mat file
    %              or selected via the GUI if not present

    if nargin < 1
      tdfn = spm_get(1, 'timediff.mat', 'Select time diff information');
    end
    if nargin < 2
      fg = [];
    end
    if isempty(fg)
      fg = spm_figure('GetWin', 'Graphics'); % use SPM figure window
      spm_figure('Clear');
    end
    if nargin < 3
      flags = '';
    end
    if isempty(flags)
      flags = ' ';
    end

    if any(flags == 'r')
      % plot realignment params
      if ~isempty(varargin)  % file name has been passed (maybe)
        fs(1).name = varargin{1};
      else
        % need to get realignment parameter file
        rwcard = 'realignment*.txt';
        [pn fn ext] = fileparts(tdfn);
        % check for realignment file in directory
        fs = dir(fullfile(pn, rwcard));
        if length(fs) > 1 || isempty(fs)
          % ask for realignment param file
          rfn = spm_get([0 1], rwcard, 'Realignment parameters');
          if ~isempty(rfn)
        fs(1).name = rfn;
          end
        end
      end
      if ~isempty(fs)
        % we do have movement parameters
        mparams = spm_load(fs(1).name);
        subpno = 5;
      else
        % we don't
        mparams = [];
        subpno = 4;
      end
    else
      subpno = 4;
    end

    load(tdfn)
    imgno = size(slicediff,1)+1;
    zno =   size(slicediff,2);
    mom = mean(globals);
    sslicediff = slicediff/mom;

    figure(fg);

    h1 = axes('position', [.1 1-.7/subpno .8 .65*1/subpno]);
    h2 = plot(2:imgno,td/mom);
    axis([0.5 imgno+0.5 -Inf Inf]);
    dxt = unique(round(get(h1, 'xtick')));
    dxt = dxt(dxt > 1 & dxt <= imgno);
    set(h1,'xtick',dxt);
    h3 = xlabel('Difference between image and reference (1)');
    h4 = ylabel('Scaled variance');
    h  = [h1; h2; h3; h4];

    h1 = axes('position', [.1 1-1.7/subpno .8 .65*1/subpno]);
    h2 = plot(2:imgno,sslicediff, 'x');
    axis([0.5 imgno+0.5 -Inf Inf]);
    set(h1,'xtick',dxt);
    h3 = xlabel('Difference between image and reference (1)');
    h4 = ylabel('Slice by slice variance');
    h  = [h; h1; h2; h3; h4];

    h1 = axes('position', [.1 1-2.7/subpno .8 .65*1/subpno]);
    h2 = plot(globals/mom);
    axis([0.5 imgno+0.5 -Inf Inf]);
    xt = unique(round(get(h1, 'xtick')));
    xt = xt(xt > 0 & xt <= imgno);
    set(h1,'xtick',xt);
    h3 = xlabel('Image number');
    h4 = ylabel('Scaled mean voxel intensity');
    h  = [h; h1; h2; h3; h4];

    h1 = axes('position', [.1 1-3.7/subpno .8 .65*1/subpno]);
    mx = max(sslicediff);
    mn = min(sslicediff);
    avg = mean(sslicediff);
    h2 = plot(avg, 'k');
    axis([0 zno+1 -Inf Inf]);
    hold on
    h3 = plot(mn, 'b');
    h4 = plot(mx, 'r');
    hold off
    h5 = xlabel('Slice number');
    h6 = ylabel('Slice variance');
    h7 = legend('Mean','Min','Max',0);
    h  = [h; h1; h2; h3; h4; h5; h6; h7];

    % realignment params
    if any(flags == 'r')
      h1 = axes('position', [.1 1-4.7/subpno .8 .65*1/subpno]);
      h2 = plot(mparams(:,1:3));
      h3 = legend('x translation','y translation','z translation',0);
      h4 = xlabel('image');
      h5 = ylabel('translations in mm');
      h  = [h; h1; h2; h3; h4; h5];
    end

    % and label with first image at bottom
    cp = get(gca,'Position');
    h1 =  axes('Position', [0 0 1 1], 'Visible', 'off');
    img1  = deblank(imgs(1,:));
    h2 = text(0.5,cp(2)/2.5,{'First image:',img1},'HorizontalAlignment','center');
    h  = [h; h1; h2];

end

