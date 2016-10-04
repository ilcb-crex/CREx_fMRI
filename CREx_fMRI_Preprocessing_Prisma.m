function CREx_fMRI_Preprocessing_Prisma
    % Example of Preprocessing implementation for SPM
    % Author: Valérie Chanoine, Research Engineer at Brain and Language
    % Institute (http://www.blri.fr/)
    % Co-authors from BLRI: Samuel Planton and Chotiga Pattadimalok
    % Co-authors from fMRI platform:   Julien Sein, Jean-Luc Anton, Bruno Nazarian and Pascal Belin from fMRI
    % platform (fMRI Center, Timone Hospital, Marseille, France)
    % Date: June 29, 2016
    
     
    
    % Initialise SPM
    spm('defaults','fmri');  
    spm_jobman('initcfg');
    
    %%========= Define parametres in a general struture 'w' ===============
    w.dataDir           = 'F:\IRM\Intermod2\Data\PREPROC';      % root directory 
    w.subjects          = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK','10TL', '11CL', '12SP', '13CC',...
                        '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'};  % subject directory (parent=dataDir)
           
    w.funcDir           =  'Functional';                                            % functional directory (parent=subject)
    w.structDir         =  'Structural';                                            % structural directory (parent=subject)
    w.sessions          = {'func01', 'func02', 'func03', 'func04', 'func05'};   	% session directory (parent=functional)
    w.T1Dir             =  'anat01';                                                % T1 directory (parent=structural)
    w.fieldMapDir       =  'fieldmap01';                                            % fieldmap directory  (parent=structural)
    w.dummy 			=  0;                                               % number of dummy files
    w.anat_ref          = 'T1_MPRAGE';                                              % part of anatomical filename
   
    
    %% Parameters from 'info.txt'of EPI files	  
    w.nSlices           = 54;     % Number of slices
    w.TR                = 1.224;  % Repetition time (s)
    w.thickness         = 2.5;    % Slice thickness (mm)
    
   	%% Parameters from Fieldmap files 
    w.SHORT_ECHO_TIME   = 4.92;        
    w.LONG_ECHO_TIME    = 7.38;   
    w.READOUT_TIME      =  42; % from EPI info file 
    
    
    %%===================================================================== 
    % Loop on subjects
    for iS=1:numel(w.subjects)
        
        fprintf('==================================================================\n');
        fprintf([w.subjects{iS} ' Preprocessing...\n']);       
        
         
        % Session order for the first level analysis
        w.sessions_First = {};
        switch (w.subjects{iS})
            
            case {'02CE', '04AC', '06EC', '08NB', '10TL', '12SP', '14EP','16SM','18EF', '20ET', '22SR', '24IA'}
                w.sessions_First    = {'func02', 'func01', 'func04', 'func03'}; 
                
            
            case {'01LA', '03NC', '05GV','07MR',  '09BK','11CL',  '13CC', '17NR', '19AK',  '21LP', '23LG', 'pilote2'}
                w.sessions_First	= {'func01', 'func02', 'func03', 'func04'};      
            
        end 
        
        % Labview prefix of behavioral files
        w.AudioLabview = {};
        switch (w.subjects{iS})
            case {'01LA','02CE'}
                w.AudioLabview  = {'COMPR1', 'PERCE1', 'COMPR2', 'PERCE2'};           
            
            case {'04AC', '06EC', '08NB', '10TL', '12SP', '14EP','16SM','18EF', '20ET', '22SR', '24IA'}
                w.AudioLabview = {'RUN2_COMPR', 'RUN1_PERCE', 'RUN4_COMPR', 'RUN3_PERCE'};
                
            case {'03NC', '05GV','07MR',  '09BK','11CL',  '13CC', '17NR', '19AK',  '21LP', '23LG', 'pilote2'}  
                w.AudioLabview = {'RUN1_COMPR', 'RUN2_PERCE', 'RUN3_COMPR', 'RUN4_PERCE'};         
        end 
  
        
        % Fieldmap prefix 
        w.prefix.fieldmap_mag   = '';   
        w.prefix.fieldmap_phase	= '';  
        switch (w.subjects{iS})
              
            case {'pilote2'} 
                   w.prefix.fieldmap_mag    =  '006_fieldmap';  
                   w.prefix.fieldmap_phase  =  '007_fieldmap'; 

            case {'06EC', '08NB', '10TL', '11CL', '12SP', '14EP', '16SM', '17NR', '18EF', '19AK', '20ET', '21LP', '22SR','23LG', '24IA'}
                   w.prefix.fieldmap_mag    =  '010_Fieldmap'; 
                   w.prefix.fieldmap_phase  =  '011_Fieldmap'; 

            case {'01LA','02CE','07MR'}               
                   w.prefix.fieldmap_mag    =  '011_Fieldmap'; 
                   w.prefix.fieldmap_phase  =  '012_Fieldmap'; 

            case {'04AC', '05GV'}
                   w.prefix.fieldmap_mag    =  '012_Fieldmap'; 
                   w.prefix.fieldmap_phase  =  '013_Fieldmap';             

            case {'13CC'} 
                   w.prefix.fieldmap_mag    =  '014_Fieldmap'; 
                   w.prefix.fieldmap_phase  =  '015_Fieldmap'; 

            case {'03NC'}  
                   w.prefix.fieldmap_mag    =  '022_Fieldmap';
                   w.prefix.fieldmap_phase  =  '023_Fieldmap'; 

            case {'09BK'}  
                   w.prefix.fieldmap_mag    =  '027_Fieldmap';      %%%%%%%%%%%%% 2 Fieldmap files !!!
                   w.prefix.fieldmap_phase  =  '028_Fieldmap'; 
        end
       
        w.subName          = w.subjects{iS};
        w.subPath          = fullfile (w.dataDir,  w.subjects{iS}); 
        w.funcPath         = fullfile (w.subPath, w.funcDir);
        w.structPath       = fullfile (w.subPath, w.structDir);
        w.fieldmapPath     = fullfile (w.structPath, w.fieldMapDir);
        w.T1Path           = fullfile (w.structPath, w.T1Dir);
        
  
        if strcmp(w.subName, 'pilote2')
            w.anat_ref       = 't1_mprage'; 
            w.READOUT_TIME 	 =  21; % from EPI info file  
        end      
         
        % Get parameters from slice onset matrix
        slice_times     = load('slice_onsets_Intermod2_VI_1225ms.mat');
 
        w.slice_times 	= slice_times.slice_onsets_Intermod2;
        w.refslice    	= w.slice_times(size(w.slice_times,2)/2);       
        
        cd(w.subPath);
           
        
        %%===== Do Preprocessing step by step  ============================     
        DoFieldMap(w, iS); 
        DoSliceTiming(w);
        DoRealignUnwarp(w);
        DoCoregister(w);
        DoSegment(w);
        DoNormalise(w);
        DoSmooth(w);
        DoExplicitMask(w);
	DoART_Outliers_regression(w);
        %%================================================================= 
    end    
end

function DoFieldMap(w, iS)

    %% magnitude    
    shortmag = spm_select('ExtFPList', w.fieldmapPath,[w.prefix.fieldmap_mag '.*\.nii$'], 1:1); 
    
    %% phase difference in radian
    phase_diff = spm_select('FPList', w.fieldmapPath,[w.prefix.fieldmap_phase '.*\.nii$']); 
        
	%% Get the T1 template  
    template	=   which('T1.nii');

    %% Get T1 structural file
    anatFile    =   spm_select('FPList', w.T1Path, ['^' w.subName   '.*' w.anat_ref '.*\.nii$']); 
    
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = cellstr(phase_diff);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(shortmag);   
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [w.SHORT_ECHO_TIME w.LONG_ECHO_TIME];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0; % 1= Magnitude Image is choosed to generate Mask Brain
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = w.READOUT_TIME; % readout time
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;  % non-EPI bases fieldmap
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;    % Jacobian use do not use
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Huttonish';    
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 15;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = cellstr(template);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;  
    
    for j=1:numel(w.sessions)
        
        switch (w.subName)
            case {'09BK'}       
                % Only for RUN1 and RUN2 
                if (j <=2)

                %% magnitude
                shortmag = spm_select('ExtFPList', w.fieldmapPath,['010_Fieldmap.*\.nii$'], 1:1); 

                %% phase difference in radian 
                phase_diff = spm_select('FPList', w.fieldmapPath,['011_Fieldmap.*\.nii$']); 
            end
        end
       
        %% Get the fisrt EPI file removing dummy files
        % f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName  '.*\.nii$']);       
        EPIfile = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName  '.*\.nii$'], w.dummy+1:w.dummy+1); 

        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(j).epi = cellstr(EPIfile);  
    end
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = cellstr(anatFile);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_1_Fieldmap.mat'),'matlabbatch'); 
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);        
end

function DoSliceTiming(w)

    clear matlabbatch;
    
    % Loop for sessions
    matlabbatch{1}.spm.temporal.st.scans = {};
    for j=1:numel(w.sessions)

        %% Get EPI raw files           
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName   '.*\.nii$']);        
        nScans = get_nii_frame(f);
        f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^' w.subName   '.*\.nii$'], w.dummy+1:nScans); 
 
        matlabbatch{1}.spm.temporal.st.scans{j} = cellstr(f); 
        %%
    end

    matlabbatch{1}.spm.temporal.st.nslices = w.nSlices;
    matlabbatch{1}.spm.temporal.st.tr = w.TR;
    matlabbatch{1}.spm.temporal.st.ta = 0; 
    matlabbatch{1}.spm.temporal.st.so = w.slice_times; 
    matlabbatch{1}.spm.temporal.st.refslice = w.refslice;      
    matlabbatch{1}.spm.temporal.st.prefix = 'a'; 
     
    save(fullfile(w.subPath, 'SPM12_matlabbatch_2_SliceTiming.mat'),'matlabbatch');     
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    %%
      
end

function DoRealignUnwarp(w)
	clear matlabbatch
  
	EPI = {};
    % Loop for sessions
    for j=1:numel(w.sessions)
        expReg = ['^vdm5.*' num2str(j) '\.nii$'];
        vdm = cellstr(spm_select('FPList', w.fieldmapPath, expReg));  % vdm (voxel depplacement map) file or phase map
        
        % Get EPI Sliced files without dummy files   
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^a' w.subName  '.*'  '.*\.nii$']);        
        nScans = get_nii_frame(f);
        newF = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^a' w.subName  '.*'  '.*\.nii$'], w.dummy+1:nScans);  
        EPI = cellstr(newF); 
        
        matlabbatch{1}.spm.spatial.realignunwarp.data(j).scans     = EPI; 
        matlabbatch{1}.spm.spatial.realignunwarp.data(j).pmscan    = vdm;
            
    end             

    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5; 
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2; % 2edegree Bspline (default value)
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
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;    %  4edegree Bspline (default value)
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_3_Realign&Unwarp.mat'),'matlabbatch');
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);       
end   

function DoCoregister(w)

    % Get T1 structural file
    a = spm_select('FPList', w.T1Path, ['^' w.subName  '.*' w.anat_ref '.*\.nii$']); 
       
    % Get mean realigned and unwarped EPI
    meanUnwarp = spm_select('FPList', fullfile(w.funcPath, w.sessions{1}), ['^meanua' w.subName  '.*\.nii$']); 
   
    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(meanUnwarp);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(a);
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
   
	save(fullfile(w.subPath, 'SPM12_matlabbatch_4_Coregister.mat'),'matlabbatch');  
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
end

function DoSegment(w)

    % Select the T1 3D image coregisterede
    coregAnat = cellstr(spm_select('FPList', w.T1Path, ['^' w.subName  '.*' w.anat_ref '.*\.nii$'])); 
 

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
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];  %native
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = tmpSoftTissue;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];  %native
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
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];   %  Forward

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_5_Segment.mat'),'matlabbatch');     
end

function DoNormalise(w)
        
	% Get Field Deformation image
    forwardDeformation = spm_select('FPList', w.T1Path, ['^y_' w.subName  '.*' w.anat_ref '.*\.nii$']); 
    
            
    % Get coregistered structural image  
    coregAnat = spm_select('FPList', w.T1Path, ['^' w.subName  '.*' w.anat_ref '.*\.nii$']);  
    
    % Get Sliced EPI images of all runs
    EPI = {};
    % Loop on sessions
    for j=1:numel(w.sessions)            
        % Get EPI Realigned files without dummy files
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^ua' w.subName   '.*\.nii$']);        
        nScans = get_nii_frame(f);
        newF = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^ua' w.subName '.*\.nii$'], w.dummy+1:nScans);  
        f = cellstr(newF);
        EPI = vertcat(EPI, f);      
    end
    
    % Get c1  c2  and c3 
    c1 = spm_select('FPList', w.T1Path, ['^c1' w.subName  '.*\.nii$']); 
    c2 = spm_select('FPList', w.T1Path, ['^c2' w.subName   '.*\.nii$']); 
    c3 = spm_select('FPList', w.T1Path, ['^c3' w.subName   '.*\.nii$']); 
            
    clear matlabbatch;         
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {coregAnat};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(EPI);    
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    
    matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
    matlabbatch{3}.spm.spatial.normalise.write.subj.resample = cellstr(c1);    
    matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
    matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;   
    
    matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
    matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cellstr(c2);    
    matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
    matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;     
    
    matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = {forwardDeformation};
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(c3);    
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [w.thickness w.thickness w.thickness];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;        

    save(fullfile(w.subPath, 'SPM12_matlabbatch_6_Normalize.mat'),'matlabbatch');   
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  

end

function DoSmooth(w)

    clear matlabbatch;
    EPI = [];
    % Get Normalized EPI files of all sessions
    for j=1:numel(w.sessions)
    
        % Get EPI Realigned files without dummy files
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^wua' w.subName '.*\.nii$']);        
        nScans = get_nii_frame(f);
        newF = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^wua' w.subName '.*\.nii$'], w.dummy+1:nScans); 
        f = cellstr(newF);
        EPI = vertcat(EPI, f);           
    end
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(EPI);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [(w.thickness*2) (w.thickness*2) (w.thickness*2)];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    save(fullfile(w.subPath, 'SPM12_matlabbatch_7_Smooth.mat'),'matlabbatch'); 
        
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);     
end

function DoExplicitMask(w)

    %% Get normalized tissus (grey and white matter, CSF
  
    wc1 = spm_select('FPList', w.T1Path, ['^wc1' w.subName  '.*\.nii$']); 
    wc2 = spm_select('FPList', w.T1Path, ['^wc2' w.subName   '.*\.nii$']); 
    wc3 = spm_select('FPList', w.T1Path, ['^wc3' w.subName   '.*\.nii$']); 

    P = [wc1; wc2; wc3];  
    
    matlabbatch{1}.spm.util.imcalc.input = cellstr(P);
    matlabbatch{1}.spm.util.imcalc.output = fullfile(w.T1Path, 'explicitMask_wc1wc2wc3_0.3.nii');
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 +i3)>0.3';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;   
       
    save(fullfile(w.subPath, 'SPM12_matlabbatch_8_Mask.mat'),'matlabbatch'); 

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
end

function DoART_Outliers_regression(w)
  normEPI = {};
    for j=1:numel(w.sessions) 
        % Get EPI Realigned files without dummy files
        normEPI {j} = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^wua' w.subName '.*\.nii$']);        
       
    end
    CREx_art_batch(normEPI, w.subName, w.subPath);
end
