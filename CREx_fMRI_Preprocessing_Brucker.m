function SPM12_Preprocessing
    % Example of Preprocessing implementation for SPM
    % Author: Valérie Chanoine, Research Engineer at Brain and Language
    % Institute (http://www.blri.fr/)
    % Co-authors: Jean-Luc Anton, Bruno Nazarian and Pascal Belin from fMRI
    % platform (fMRI Center, Timone Hospital, Marseille, France)
    % Date: June 15, 2015
    
    
    % Initialise SPM
    spm('defaults','fmri');  
    spm_jobman('initcfg');
    
    % Define a general struture 'w'
    
    w.dataDir           = 'F:\blog_CREX\fMRI\BLRIproject';  % root directory
    w.groups            = {'Temoins'};                      % group directory if necessary (parent=dataDir)
    w.subjects          = {'S01','S02'};                    % subject directory (parent=group)
     
    w.funcDir           =  'Functional';                    % functional directory (parent=subject)
    w.structDir         =  'Structural';                    % structural directory (parent=subject)
    w.sessions          = {'func01', 'func02'};             % session directory (parent=functional)

    w.T1Dir             =  'anat01';                        % T1 directory (parent=structural)
    w.fieldMapDir       =  'fieldmap01';                    % fieldmap directory  (parent=structural)
    w.stimDir           =  'Stim';                          % stimulation directory (parent=subject)
    w.prefix            =  'prefixName';                    % prefix of all raw data files
    w.dummy 			=  5;								% number of dummy files
    w.nSlices           = 30;
    w.RT                = 2;
 
    % Loop on groups
    for iG=1:numel(w.groups) 

        % Loop on subjects
        for iS=1:numel(w.subjects)  
      
            w.subPath          = fullfile (w.dataDir, w.groups{iG},  w.subjects{iS}); 
            w.funcPath         = fullfile (w.subPath, w.funcDir);
            w.structPath       = fullfile (w.subPath, w.structDir);
            w.stimPath         = fullfile (w.subPath, w.stimDir);
            w.fieldmapPath     = fullfile (w.structPath, w.fieldMapDir);
            w.T1Path           = fullfile (w.structPath, w.T1Dir);
            
            
            % Do Preprocessing step by step
            DoFieldMap(w);
            DoRealignUnwarp(w);
            DoCoregister(w);
            DoSliceTiming(w);
            DoSegmentSPM12(w);
            DoNormalise(w);
            DoSmooth(w);
            
        end
    end
end

function DoFieldMap (w)

    % Define Fieldmap Parameters 
    SHORT_ECHO_TIME = 3.7;        
    LONG_ECHO_TIME = 8.252;   
    READOUT_TIME =  33.178;

    % Get Fieldmap files  
    shortreal   =   spm_select('FPList', w.fieldmapPath,['^' w.prefix '.*_fieldmap01_echo01_real\.nii$']); 
    shortimag   =   spm_select('FPList', w.fieldmapPath,['^' w.prefix '.*_fieldmap01_echo01_imag\.nii$']); 
    longreal    =   spm_select('FPList', w.fieldmapPath,['^' w.prefix '.*_fieldmap01_echo02_real\.nii$']); 
    longimag    =   spm_select('FPList', w.fieldmapPath,['^' w.prefix '.*_fieldmap01_echo02_imag\.nii$']);    
    
    % Get the T1 template  
    template	=   fullfile(spm('Dir'),'templates', 'T1.nii');

    % Get T1 structural file
    anatFile    =   spm_select('FPList', w.T1Path, ['^' w.prefix '.*anat01\.nii$']);      

    % Get the fisrt EPI file removing dummy files
    f = cellstr(spm_select('FPList',  fullfile(w.funcPath, w.sessions{1}), ['^' w.prefix '.*\.nii$']));    
    EPIfile = f{w.dummy,1};


    clear matlabbatch
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.shortreal  =   cellstr(shortreal);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.shortimag  =   cellstr(shortimag);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.longreal   =   cellstr(longreal);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.longimag   =   cellstr(longimag);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.et = [SHORT_ECHO_TIME LONG_ECHO_TIME];
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.maskbrain = 1; % Magnitude Image is choosed to generate Mask Brain
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.blipdir = 1;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.tert = READOUT_TIME; % readout time
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.epifm = 0;  % non-EPI bases fieldmap
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.ajm = 0;    % Jacobian use do not use
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.uflags.method = 'Huttonish';
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.uflags.pad = 15;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.template = {template};
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.defaults.defaultsval.mflags.reg = 0.02;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.session.epi = cellstr(EPIfile);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.matchvdm = 0;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.sessname = 'session';
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.anat = cellstr(anatFile);
    matlabbatch{1}.spm.tools.fieldmap.realimag.subj.matchanat = 0; 

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch); 
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_1_Fieldmap.mat'),'matlabbatch'); 
       
end

function DoRealignUnwarp(w)
        
    vdm = spm_select('FPList', w.FieldmapPath, '^vdm.*\.nii$');  % vdm (voxel depplacement map) file or phase map
  
    P = [];
    % Loop for sessions
    for j=1:numel(w.sessions)
        % Get EPI Sliced files
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}),['^' w.prefix '.*\.nii$']);    
        newF = f(w.dummy:end,:) % Remove dummy scans
        P = [P; newF ] ;               

    end        

    clear matlabbatch
    matlabbatch{1}.spm.spatial.realignunwarp.data.scans     = cellstr(P);  
    matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan    = cellstr(vdm);
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 3;  % NORMALEMENT 4 PAR DEFAUT
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 1; % Trilinear interpolation
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];      
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 1;    %  modifié
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_2_Realign&Unwarp.mat'),'matlabbatch');
end   

function DoCoregister(w)

    % Get T1 structural file
    a = spm_select('FPList', w.T1Path, ['^' w.prefix '.*anat01\.nii$']); 
           
    % Get mean realigned and unwarped EPI
    meanUnwarp = spm_select('FPList', fullfile(w.funcPath, w.sessions{1}), ['^meanu' w.prefix '.*\.nii$']); 
    

    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(meanUnwarp);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(a);
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_3_Coregister.mat'),'matlabbatch');
    
end

function DoSliceTiming(w)

    clear matlabbatch;
    
    % Loop for sessions
    matlabbatch{1}.spm.temporal.st.scans = {};
    for j=1:numel(w.sessions)

        % Get EPI raw files
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^u' w.prefix '.*\.nii$']); 
        matlabbatch{1}.spm.temporal.st.scans{j} = cellstr(f);  

    end

    sliceOrder = 1:1:w.nSlices; % (Slice Order = ascending)
    refslice = sliceOrder(floor ( w.nSlices/2));

    matlabbatch{1}.spm.temporal.st.nslices = w.nSlices;
    matlabbatch{1}.spm.temporal.st.tr = w.RT;
    matlabbatch{1}.spm.temporal.st.ta = w.RT-(w.RT/ w.nSlices); 
    matlabbatch{1}.spm.temporal.st.so = sliceOrder; 
    matlabbatch{1}.spm.temporal.st.refslice = refslice;      
    matlabbatch{1}.spm.temporal.st.prefix = 'a'; 
     
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_4_SliceTiming.mat'),'matlabbatch');          
end

function DoSegmentSPM12(w)

    % Select the T1 3D image coregistered
    coregAnat = cellstr(spm_select('FPList', w.T1Path, ['^' w.prefix '.*anat01\.nii$'])); 


    % Get template of each cerebral tissue
    tmpGM           = {fullfile(spm('Dir'),'tpm', 'TPM.nii,1')};
    tmpWM           = {fullfile(spm('Dir'),'tpm', 'TPM.nii,2')};
    tmpCSF          = {fullfile(spm('Dir'),'tpm', 'TPM.nii,3')};
    tmpBone         = {fullfile(spm('Dir'),'tpm', 'TPM.nii,4')};      
    tmpSoftTissue   = {fullfile(spm('Dir'),'tpm', 'TPM.nii,5')};               
    tmpAirBck       = {fullfile(spm('Dir'),'tpm', 'TPM.nii,6')};        

    clear matlabbatch;                             
    matlabbatch{1}.spm.spatial.preproc.channel.vols = coregAnat;
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; % Bias regularisation light
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;  % 60 mm cutoff
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];     
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = tmpGM;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];  %native
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];       
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = tmpWM;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];  %native
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];      
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = tmpCSF;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];  %native
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];      
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = tmpBone;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = tmpSoftTissue;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];        
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = tmpAirBck;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];     
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];   % Inverse + Forward

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_5_Segment.mat'),'matlabbatch');     
end

function DoNormalise(w)
        
	% Get Field Deformation image
	forwardDeformation = spm_select('FPList', w.T1Path, ['^y_' w.prefix '.*anat01\.nii$']); 
    
    % Get coregistered structural image    
    coregAnat = spm_select('FPList', w.T1Path, ['^' w.prefix '.*anat01\.nii$']); 
    
    % Get Sliced EPI images of all runs
    P = [];
    % Loop on sessions
    for j=1:numel(w.sessions)
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}),'^au.*\.nii$');  
        P = [P; f ];               

    end  
    clear matlabbatch; 
        
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {coregAnat};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
                                                      90 90 108];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};



    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(P);    
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
                                                      90 90 108];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 1;    

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_6_Normalize.mat'),'matlabbatch');   
end

function DoSmooth(w)

    clear matlabbatch;
    P = [];
    % Get Normalized EPI files of all sessions
    for j=1:numel(w.sessions)
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}),'^wau.*\.nii$');  
        P = [P; f ];               
    end
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(P);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_7_Smooth.mat'),'matlabbatch');  
end

