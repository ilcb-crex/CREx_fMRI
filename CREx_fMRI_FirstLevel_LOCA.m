function SPM12_FirstLevel_LOCA
    %%
    close all; clear; clc;  
    %%
    
    % Initialise SPM
    spm('defaults','fmri');  
    spm_jobman('initcfg');
    
    %%====================== Define a general struture 'w' =================
%     w.dataDir           = 'F:\IRM\Intermod2\Data';      % root directory
	w.dataDir           = 'D:\fMRI-data\Data';      % root directory

    w.subjects          = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK', '10TL', '11CL', '12SP',...
        '13CC', '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'};  % subject directory (parent=dataDir)

    w.subjects          = {'12SP'};  % subject directory (parent=dataDir)

    w.funcDir           =  'Functional';                                  	% functional directory (parent=subject)
    w.structDir         =  'Structural';                                    % structural directory (parent=subject)
    w.sessions          =  {'func05'};                                      % session directory (parent=functional)
    w.T1Dir             =  'anat01';                                        % T1 directory (parent=structural)
    w.fieldMapDir       =  'fieldmap01';                                    % fieldmap directory  (parent=structural)
    w.stimDir           =  'Stim';  
    
    % stimulation directory (parent=subject)
    w.dummy 			=  0;	  % number of dummy files
    w.TR                = 1.224;  % Repetition time (s)
    
    %====@@@@@====% Initialize WaitBar %====@@@@@====% 
     tic; s = clock; h = waitbar(0,{strrep(mfilename, '_', '-')...
        ' '...
        ['Loop on subjects - Subject ' w.subjects{1} ' in progress']...
        ' '...
        'will be updated after each subject...'});   
    %====@@@@@==========@@@@@=@@@@@==========@@@@@====%   
    
    % Loop on subjects
    for iS=1:numel(w.subjects)
           
        switch (w.subjects{iS})
            case {'01LA','02CE','pilote2'}
                w.AudioLabview      = {'Localizer'};
            otherwise
                w.AudioLabview      = {'LOCAVISUEL'};
        end    
        
        w.subName          = w.subjects{iS};

        w.subPath          = fullfile (w.dataDir,  'PREPROC', w.subjects{iS}); 
        w.funcPath         = fullfile (w.subPath, w.funcDir);
        w.stimPath         = fullfile (w.dataDir,  'PREPROC', 'Labview_Files',  w.subName, w.stimDir);
        w.structPath       = fullfile (w.subPath, w.structDir);       
        w.fieldmapPath     = fullfile (w.structPath, w.fieldMapDir);
        w.T1Path           = fullfile (w.structPath, w.T1Dir);          
       
        % Get parameters from slice onset matrix
        slice_times = load('slice_onsets_Intermod2_VI_1225ms.mat');
   
        w.slice_times   = slice_times.slice_onsets_Intermod2;
        w.refslice      = w.slice_times(size(w.slice_times,2)/2);       
        
        cd(w.subPath);      
        
        
        %%===== Do First Level
        bEvent = false;         % booleen for event design vs block design 
        bWithArt = false;       % booleen for model with vs without "Art Detection" regressors
        
        DoFirstLevel(w, bEvent, bWithArt);
        %%=====       
               
            
        %====@@@@@====% Update WaitBar & remaining time %====@@@@@====%
        if iS<numel(w.subjects)
            remainTime = (etime(clock,s)/iS) * (numel(w.subjects)-iS);
            endTime = datetime('now','Format','HH:mm')+seconds(remainTime);               
            waitbar(iS/numel(w.subjects),h,{strrep(mfilename, '_', '-') ' '...
                ['Loop on subjects - Subject ' w.subjects{iS+1} ' in progress']... 
                ['Remaining time < ',num2str( remainTime /60,'%.1f' ),' minutes' ]...
                ['End time ~ ' char(endTime)] });
        end
        %====@@@@@==========@@@@@=============@@@@@==========@@@@@====%
  
    end    
    waitbar(1,h,{strrep(mfilename, '_', '-') ' ' '========== Done !!! ==========' ' ' ['Took ' num2str(round(toc/60)) ' minutes'] });
end

function DoFirstLevel(w, bEvent, bWithArt)

    % Make new directory 'First level'
    if(bEvent)       
        if (bWithArt)
            firstDir = fullfile (w.dataDir, 'FIRST_LEVEL', 'LOCA', 'FirstLevel_Event_Art',  w.subName); 
        else
            firstDir = fullfile (w.dataDir, 'FIRST_LEVEL', 'LOCA', 'FirstLevel_Event',  w.subName); 
        end
    else
        if (bWithArt)
            firstDir = fullfile (w.dataDir, 'FIRST_LEVEL', 'LOCA', 'FirstLevel_Block_Art',  w.subName); 
        else
            firstDir = fullfile (w.dataDir, 'FIRST_LEVEL', 'LOCA', 'FirstLevel_Block',  w.subName); 
        end
    end
        
    if isdir (firstDir)
        delete ([firstDir '/*']);
    end
    mkdir(firstDir); 
           
    %==============================================================%
    %  fMRI model specification
    %==============================================================% 
      
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstDir}; 
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';    
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = w.TR;    
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = numel(w.slice_times)/3;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = numel(w.slice_times)/3/2; 
    
    for j=1:numel(w.sessions)
        %% Get the file of head movements (and prepare 'zeros' for regressors of no interest for contrasts)

        if (bWithArt)
            rpF = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^art_regression_outliers_and_movement_' '.*\.mat$']);  
            load(rpF); no_interest_reg{j} = zeros(1,size(R,2));
        else
            rpF = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), ['^rp' '.*\.txt$']);
            no_interest_reg{j} = zeros(1,size(load(rpF),2));
        end
        
        %% Get stimulation onsets 
        onset_file = {};

        if (bEvent)
            onset_file = cellstr(spm_select('FPList', w.stimPath, [w.AudioLabview{j}  '.*\event.mat$'])); 
        else
            onset_file = cellstr(spm_select('FPList', w.stimPath, [w.AudioLabview{j} '.*block.mat$'])); 
        end 

        %% Get EPI smoothed images         
        f = spm_select('FPList',  fullfile(w.funcPath, w.sessions{j}), [ '^swua'  w.subName  '.*\.nii$' ]);        
        nScans = get_nii_frame(f);
        EPI = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions{j}), ['^swua' w.subName  '.*\.nii$'], w.dummy+1:nScans); 
%%
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).scans = cellstr(EPI);  
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = onset_file;
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).regress = struct('name', {}, 'val', {});        
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg = {rpF}; 
        matlabbatch{1}.spm.stats.fmri_spec.sess(j).hpf = 128;       
        
    end   
    explicitMask = spm_select('FPList', w.T1Path, '^explicitMask_wc1wc2wc3_03.nii$'); 
    
%%	
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3; % Implicit mask threshold (default 0.8) 
    matlabbatch{1}.spm.stats.fmri_spec.mask = {explicitMask};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';       
          
    %==============================================================%
    %  Model Estimation
    %==============================================================%   

     matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(firstDir,'SPM.mat'));
     matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;    

    %==============================================================%
    %  Contrast manager
    %==============================================================%   
    
    % 4 Conditions 
    % C1 : FIXATION     
    % C2 : WORDS
    % C3 : NONWORDS
    % C4 : RESPONSE  
    nconds = 4;
    conds       = eye(nconds);
    Fixation    = conds(1,:);
    Words       = conds(2,:);
    Nonwords    = conds(3,:);
    Response    = conds(4,:);
 
	matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(firstDir,'SPM.mat'));   
 
    %% Contrasts T (betas)
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FIXATION';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [Fixation no_interest_reg{1}];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'WORDS';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [Words no_interest_reg{1}];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'NONWORDS';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [Nonwords no_interest_reg{1}];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'RESPONSE';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [Response no_interest_reg{1}];
     matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
     
    %% Contrasts T (comparisons)
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'WORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [(Words - Fixation)  no_interest_reg{1}]/2; % Division par 2 pas indispensable...
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'NONWORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [(Nonwords - Fixation)  no_interest_reg{1}]/2;
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'WORDS&NONWORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [((Words+Nonwords)/2-Fixation) no_interest_reg{1}]/2;
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'WORDS > NONWORDS';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [(Words - Nonwords)  no_interest_reg{1}]/2;
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    %% Contrasts F
 
    % Contrast F : Movement & outliers effect
    mat_size = max(cellfun(@(x) numel(x),no_interest_reg)); % Array size (rows) for "F-no_interest_reg" contrast
    matlabbatch{3}.spm.stats.con.consess{9}.fcon.name = 'F-no_interest_reg';
    matlabbatch{3}.spm.stats.con.consess{9}.fcon.weights = [zeros(mat_size,nconds)  eye(mat_size,numel(no_interest_reg{1}))];
    matlabbatch{3}.spm.stats.con.consess{9}.fcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{10}.fcon.name = 'F-TOUTES-CONDITIONS';
    matlabbatch{3}.spm.stats.con.consess{10}.fcon.weights = [conds zeros(nconds,numel(no_interest_reg{1}))];
    matlabbatch{3}.spm.stats.con.consess{10}.fcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete = 1;  
    
%%  
	if (bWithArt)
         if(bEvent)
            save(fullfile(firstDir, 'SPM12_matlabbatch_9_FirstLevel_Event_ART.mat'),'matlabbatch');
        else
            save(fullfile(firstDir, 'SPM12_matlabbatch_9_FirstLevel_Block_ART.mat'),'matlabbatch');
        end
    else
        if(bEvent)
            save(fullfile(firstDir, 'SPM12_matlabbatch_9_FirstLevel_Event.mat'),'matlabbatch');
        else
            save(fullfile(firstDir, 'SPM12_matlabbatch_9_FirstLevel_Block.mat'),'matlabbatch');
        end  
	end
%%         
	spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
    
end
