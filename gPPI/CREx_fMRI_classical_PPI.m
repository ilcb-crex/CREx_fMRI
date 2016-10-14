function CREx_fMRI_classical_PPI
%==========================================================================
% CLASSICAL PSYCHO-PHYSIOLOGICAL INTERACTION (PPI) 
%
% Script for Classical implementation of PPI using SPM12
% Author: Valérie Chanoine, Research Engineer at Brain and Language
% Co-authors from BLRI: Samuel Planton and Chotiga Pattadimalok
%==========================================================================


    %======================================================================
    % Initialise SPM
    %======================================================================
      
    spm('Defaults','fMRI');  
    spm_jobman('initcfg');  
    
    %======================================================================
    % Settings
    %======================================================================
        
    bEvent = true;
    
    w.dataPath           = 'F:\IRM\Intermod2\Data';      % data path  
	w.subjects          = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK',...
        '10TL', '11CL', '12SP', '13CC', '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'};  % subject directory (parent=dataDir)          
    w.funcDir           =  'Functional';     % functional directory (parent=subject)   
    w.structDir         =  'Structural';     % structural directory (parent=subject)
    w.T1Dir             =  'anat01';         % T1 directory (parent=structural)
    w.stimDir           =  'Stim';  
    w.TR                = 1.224;  % Repetition time (s)


    % Loop on subjects
	for iS=1:numel(w.subjects)
       
        %==================================================================
        % Subject Settings
        %==================================================================      
              
        w.subName          = w.subjects{iS};
        w.subPath          = fullfile (w.dataPath,  'PREPROC', w.subjects{iS}); 
        
        w.funcPath         = fullfile (w.subPath, w.funcDir);
        w.structPath       = fullfile (w.subPath, w.structDir);       
        w.T1Path           = fullfile (w.structPath, w.T1Dir);        
        w.stimPath         = fullfile (w.dataPath,  'Labview_Files',  w.subName, w.stimDir);  
        
        % Redefine session order for each subject      
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
        
        % Get parameters from slice onset matrix
        slice_times     = load('slice_onsets_Intermod2_VI_1225ms.mat');
        w.slice_times   = slice_times.slice_onsets_Intermod2;
          
   
        % Define directory 'GLM_E'
        if(bEvent)
            w.GLMDir = fullfile (w.dataPath, 'AUDIO', 'GLM_Event',  w.subName); 
        else
            w.GLMDir = fullfile (w.dataPath, 'AUDIO', 'GLM_Block',  w.subName); 
        end
        
        SpecifyGLM(w, bEvent);
         
        VOI1.name = 'vOT';
        VOI1.centre = [-47.5 -48.5 -14.5];
        VOI1.radius = 6;    
        VOI1.spm.contrast = 4; % COMPs_s + PERC_s_s   
        VOI1.ppi.name = [VOI1.name 'x(COMPs-PERCs)'];   
        
        VOI2.name = 'IFG_Tri';
        VOI2.centre = [-42 14 19];
        VOI2.radius = 6;  
        VOI2.spm.contrast = 4; % COMPs_s + PERC_s_s       
               
        DoClassicalPPI(w, bEvent, VOI1,VOI2);          
    end
 
end
     
function SpecifyGLM(w, bEvent)

        % Create new directory 'GLM_E'
        if isdir (w.GLMDir)
            delete ([w.GLMDir '/*']);
        end
        mkdir(w.GLMDir);

        %======================================================================
        % GLM SPECIFICATION, ESTIMATION, AND INFERENCE
        %======================================================================  
        SpecifyAndEstimateModel(w, bEvent);
        DefineContrats(w, bEvent);
end

function SpecifyAndEstimateModel(w, bEvent)
    %======================================================================
    % GLM SPECIFICATION, ESTIMATION, AND INFERENCE
    %======================================================================  
        
    clear matlabbatch;
    % Directory
    %----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.dir = {w.GLMDir}; 
    
    % Timing
    %----------------------------------------------------------------------    
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';    
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = w.TR;    
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = numel(w.slice_times)/3;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = numel(w.slice_times)/3/2;     
       
    % Sessions
    %----------------------------------------------------------------------    

     onset_file = {};
     EPI = {};
     rpF = {};
     
     for j=1:numel(w.sessions_First)
        % Get the file of head movements  
        rpF{j} = load( spm_select('FPList',  fullfile(w.funcPath, w.sessions_First{j}), ['^rp' '.*\.txt$']));
              

        % Get stimulation onsets 
        if (bEvent)           
            onset_file{j} = load(spm_select('FPList', w.stimPath, [w.AudioLabview{j}  '.*\event.mat$'])); 

        else
            %onset_file{j} = cellstr(spm_select('FPList', w.stimPath, [w.AudioLabview{j} '.*block.mat$'])); 
            onset_file{j} = load(spm_select('FPList', w.stimPath, [w.AudioLabview{j}  '.*block.mat'])); 
        end 

        % Concatenate EPI smoothed images                
        f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions_First{j}), ['^swua' w.subName  '.*\.nii$'], inf);    
        EPI = vertcat(EPI, cellstr(f));  
     end
     %%
       
    % Get onset file of an unique session
    %----------------------------------------------------------------------   
    onset_File = ConcatenateOnsetsInOneSession(w, onset_file)
    
    % Get headMvts file of an unique session
    %----------------------------------------------------------------------   
  
    R = [rpF{1} ; rpF{2}; rpF{3} ;rpF{4} ];  
    save(fullfile(w.GLMDir, 'HeadMvts.mat'), 'R');   
    HeadMvts_File = cellstr(spm_select('FPList', w.GLMDir, ['^HeadMvts' '*.\mat$']));
     
        
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = EPI;  
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = onset_File;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});        
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = HeadMvts_File; 

    % High-pass filter
    %------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;              
  
      
    % Explicit mask
    %---------------------------------------------------------------------- 
    
    explicitMask = spm_select('FPList', w.T1Path, '^explicitMask_wc1wc2wc3_0.3.nii$'); 
      
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3; % Implicit mask threshold (default 0.8) 
    matlabbatch{1}.spm.stats.fmri_spec.mask = {explicitMask};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';          
    
   
    
    % Block effect to model the different sessions
    % %--------------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'Session 1';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = kron([1 0 0 0]',ones(353,1));
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'Session 2';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = kron([0 1 0 0 ]',ones(353,1));
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Session 3';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = kron([0 0 1 0]',ones(353,1));
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'Session 4';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = kron([0 0 0 1]',ones(353,1));
     
    
        
    %======================================================================
    %  Model Estimation
    %======================================================================   

	matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat'));
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;   
    
    save(fullfile(w.GLMDir, 'SPM12_matlabbatch_1_classicalPPI_GLM_SpecMod.mat'),'matlabbatch');

	spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
end

function onset_path = ConcatenateOnsetsInOneSession(w, onset_file)
      
    %======================================================================
    % Concatenate names
    %======================================================================
    names = horzcat(onset_file{1}.names, onset_file{2}.names);
 
    %======================================================================
    % Concatenate onsets and durations
    %======================================================================
    NCOND               = numel(names)/2;
    total_duration      = 353 * w.TR;
    
    onsets_sess1_3      = {};
    onsets_sess2_4      = {} ;

    durations_sess1_3   = {};
    durations_sess2_4   = {} ; 

    for iC=1:NCOND
       onsets_sess1_3{iC} = [onset_file{1}.onsets{iC} ; (onset_file{3}.onsets{iC} + 2* total_duration)];
       onsets_sess2_4{iC} = [(onset_file{2}.onsets{iC} + total_duration) ; (onset_file{4}.onsets{iC} + 3 *total_duration)];

       durations_sess1_3{iC} = [onset_file{1}.durations{iC} ; onset_file{3}.durations{iC}];
       durations_sess2_4{iC} = [onset_file{2}.durations{iC}  ; onset_file{4}.durations{iC}];
    end
    
    onsets      = [onsets_sess1_3 onsets_sess2_4];
    durations   = [durations_sess1_3 durations_sess2_4];
 

    durations{5} = 0;
    durations{10} = 0;
    
    %======================================================================
    % Save the concatenation in one file
    %======================================================================  
    save(fullfile(w.GLMDir, 'onsets_OneSession.mat'),'names', 'onsets', 'durations');    
    onset_path = cellstr(spm_select('FPList', w.GLMDir, ['^onsets_OneSession' '*.\mat$']));    
  end


function DefineContrats(w, bEvent)
    %======================================================================
    % MODEL INFERENCE
    %======================================================================  
        
    clear matlabbatch;
    % Directory
    %----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat')); 

    % Inference
    %----------------------------------------------------------------------
    % CONDITIONS
    %%
    nconds          = 10;
    conds           = eye(nconds);
    no_Cond         = zeros(1,nconds);
    
    % 10
    REPOS_COMPR_SILENCE         = conds(1,:);
    REPOS_COMPR_BRUIT           = conds(2,:);
    TACHE_COMPR_SILENCE         = conds(3,:); 
    TACHE_COMPR_BRUIT           = conds(4,:);
    REPONSE_MANUELLE_COMPR    	= conds(5,:);
    REPOS_PERCE_SILENCE         = conds(6,:);
    REPOS_PERCE_BRUIT           = conds(7,:);
    TACHE_PERCE_SILENCE         = conds(8,:);
    TACHE_PERCE_BRUIT           = conds(9,:);
    REPONSE_MANUELLE_PERC     	= conds(10,:);  

   
    Sessions = zeros(1,4);
    HeadMvts = zeros(1,6);
  
    COMPs_s = TACHE_COMPR_SILENCE - REPOS_COMPR_SILENCE;
    PERCs_s = TACHE_PERCE_SILENCE - REPOS_PERCE_SILENCE;
    Tasks_s = COMPs_s + PERCs_s;
%%
    Fsession =  [eye(10,10) zeros(10,10)];
    %%
                                                         
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'F-TOUTES-CONDITIONS';
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = Fsession;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'COMPs_s';
	matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [COMPs_s Sessions HeadMvts];
	matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
 
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'PERCs_s';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [PERCs_s Sessions HeadMvts];
	matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none'; 
    

    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Tasks_s';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [Tasks_s Sessions HeadMvts];
	matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';   
  
    save(fullfile(w.GLMDir, 'SPM12_matlabbatch_2_classicalPPI_GLM_Contrasts.mat'),'matlabbatch');
   
    % Run Job
    %----------------------------------------------------------------------
	spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);    
  
end
  
function VOI_ExtractTimeSeries(w, VOI)
    %======================================================================
    % VOI: EXTRACTING TIME SERIES
    %======================================================================       
 
    clear matlabbatch;

    % Directory
    %----------------------------------------------------------------------
    matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat')); 
    matlabbatch{1}.spm.util.voi.adjust = 1;
    matlabbatch{1}.spm.util.voi.session = 1; %%%% 
    matlabbatch{1}.spm.util.voi.name = VOI.name;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
    matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = VOI.spm.contrast;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = VOI.centre;
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = VOI.radius;
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm = 1;
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';   
   
    % Run Job
    %----------------------------------------------------------------------
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);   
   
end

function DoClassicalPPI(w, bEvent, VOI1, VOI2)


    VOI_ExtractTimeSeries(w, VOI1);

    %======================================================================
    % Do Classical PSYCHO-PHYSIOLOGIC INTERACTION
    %======================================================================  
    clear matlabbatch


    % GENERATE PPI STRUCTURE
    %======================================================================
    matlabbatch{1}.spm.stats.ppi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat')); 
    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = cellstr(fullfile(w.GLMDir, ['VOI_' VOI1.name '_1.mat']));
    matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 1 -1; 3 1 1; 6 1 -1; 8 1 1]; % COMPs_s > PERCs_s
    matlabbatch{1}.spm.stats.ppi.name = VOI1.ppi.name;
    matlabbatch{1}.spm.stats.ppi.disp = 0;


    % OUTPUT DIRECTORY
    %==========================================================================

    w.PPIDir = fullfile (w.dataPath, 'AUDIO', 'PPI',  w.subName, 'PPI'); 
    if ~exist(w.PPIDir)
        mkdir(w.PPIDir);
    end
  
    
    % MODEL SPECIFICATION
    %======================================================================

    % Directory
    %--------------------------------------------------------------------------
    matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(w.PPIDir);

    
    % Timing
    %----------------------------------------------------------------------    
    matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';    
    matlabbatch{2}.spm.stats.fmri_spec.timing.RT = w.TR;   
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = numel(w.slice_times)/3;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = numel(w.slice_times)/3/2;  
    
    
    % Get the file of ppi structure 
    ppiF = fullfile( w.GLMDir, ['PPI_' VOI1.ppi.name '.mat']);  
    
    
    EPI = {};
	for iSess=1:numel(w.sessions_First) 
 
        % Get EPI smoothed images                
        f = spm_select('ExtFPList',  fullfile(w.funcPath, w.sessions_First{iSess}), ['^swua' w.subName  '.*\.nii$'], inf);    
        EPI = vertcat(EPI, cellstr(f));  
    end
    
    % Session
    %--------------------------------------------------------------------------
    matlabbatch{2}.spm.stats.fmri_spec.sess.scans = EPI;
    
    
   % Regressors
    %--------------------------------------------------------------------------
	% Create 'multi_block_regressors.mat'
    val = kron([1 0 0 0]',ones(353,1))
    val2 = kron([0 1 0 0]',ones(353,1));
    val3 = kron([0 0 1 0]',ones(353,1));
    R = [];
    R = [val val2 val3];
    
    save(fullfile( w.GLMDir,'multi_block_regressors.mat'), 'R');  
    
    matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = { ppiF;...
    fullfile(w.GLMDir,'multi_block_regressors.mat')};  
    
    
    % High-pass filter
    %--------------------------------------------------------------------------
    matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;


    % MODEL ESTIMATION
    %==========================================================================
    matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(w.PPIDir,'SPM.mat'));


    % INFERENCE
    %==========================================================================

    matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(w.PPIDir,'SPM.mat'));
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'PPI-Interaction';
    matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights =  [1 0 0 0 0 0 0];

    % RESULTS
    %==========================================================================
    matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(w.PPIDir,'SPM.mat'));
    matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{5}.spm.stats.results.conspec.thresh = 0.01;
    matlabbatch{5}.spm.stats.results.conspec.extent = 3;
    matlabbatch{5}.spm.stats.results.print = false;
        

    save(fullfile(w.GLMDir, ' SPM12_matlabbatch_3_classicalPPI_PPImodel.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
      

            
    % JUMP TO V5 AND OVERLAY ON A STRUCTURAL IMAGE
    %----------------------------------------------------------------------

    w.PPIDir = fullfile (w.dataPath, 'AUDIO', 'PPI',  w.subName, 'PPI'); 

	% Evaluated fields in xSPM (input)
    xSPM.swd = w.PPIDir; 
    xSPM.title = 'PPI-Interaction'; 
    xSPM.Ic = 1;     
    xSPM.n =1;
    xSPM.Im = [];   
    xSPM.pm = []; 
    xSPM.Ex = []; 
    xSPM.u = 0.01;
    xSPM.k = 3; 
    xSPM.thresDesc = 'none';

    [SPM,xSPM] = spm_getSPM(xSPM);
    xSPM.thresDesc = 'none';
    [hReg,xSPM,SPM] = spm_results_ui('Setup',xSPM);


    template = fullfile(spm('Dir'),'canonical', 'single_subj_T1.nii');

    spm_mip_ui('SetCoords',VOI2.centre);

    spm_sections(xSPM,findobj(spm_figure('FindWin','Interactive'),'Tag','hReg'),...
            template);
       

    % PSYCHO-PHYSIOLOGIC INTERACTION GRAPH
    %----------------------------------------------------------------------
    VOI_ExtractTimeSeries(w, VOI2);
  
    clear matlabbatch;
    
    VOI1_file = cellstr(fullfile(w.GLMDir, ['VOI_' VOI1.name '_1.mat']));
    VOI2_file = cellstr(fullfile(w.GLMDir, ['VOI_' VOI2.name '_1.mat']));
    
      
    % GENERATE PPI STRUCTURE: vOTxCOMP
    %----------------------------------------------------------------------  
    matlabbatch{1}.spm.stats.ppi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat')); 
    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = VOI1_file;
    matlabbatch{1}.spm.stats.ppi.type.ppi.u = [2 1 1];
    matlabbatch{1}.spm.stats.ppi.name = 'vOTxCOMP';
    matlabbatch{1}.spm.stats.ppi.disp = 0;

    % GENERATE PPI STRUCTURE: vOTxPERC
    %----------------------------------------------------------------------
    matlabbatch{2}.spm.stats.ppi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat'));
    matlabbatch{2}.spm.stats.ppi.type.ppi.voi = VOI1_file;
    matlabbatch{2}.spm.stats.ppi.type.ppi.u = [3 1 1];
    matlabbatch{2}.spm.stats.ppi.name = 'vOTxPERC';
    matlabbatch{2}.spm.stats.ppi.disp = 0;

    % GENERATE PPI STRUCTURE: IFGtxCOMP
    %----------------------------------------------------------------------  
    matlabbatch{3}.spm.stats.ppi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat'));
    matlabbatch{3}.spm.stats.ppi.type.ppi.voi = VOI2_file;
    matlabbatch{3}.spm.stats.ppi.type.ppi.u = [2 1 1];
    matlabbatch{3}.spm.stats.ppi.name = 'IFGtxCOMP';
    matlabbatch{3}.spm.stats.ppi.disp = 0;

    % GENERATE PPI STRUCTURE: IFGtxPERC
    %----------------------------------------------------------------------  
    matlabbatch{4}.spm.stats.ppi.spmmat = cellstr(fullfile(w.GLMDir,'SPM.mat'));
    matlabbatch{4}.spm.stats.ppi.type.ppi.voi = VOI2_file;
    matlabbatch{4}.spm.stats.ppi.type.ppi.u = [3 1 1];
    matlabbatch{4}.spm.stats.ppi.name = 'IFGtxPERC';
    matlabbatch{4}.spm.stats.ppi.disp = 0;

    spm_jobman('run',matlabbatch);

    % PLOT THE PPI INTERACTION VECTORS UNDER EACH EXPer CONDITION
    %----------------------------------------------------------------------  
    load('PPI_vOTxCOMP.mat');   PPI_vOTxCOMP = PPI;
    load('PPI_vOTxPERC.mat');   PPI_vOTxPERC  = PPI;
    load('PPI_IFGtxCOMP.mat');  PPI_IFGtxCOMP = PPI;
    load('PPI_IFGtxPERC.mat'); 	PPI_IFGtxPERC  = PPI;

    figure;
    plot(PPI_vOTxCOMP.ppi,PPI_IFGtxCOMP.ppi,'k.');
    hold on
    plot(PPI_vOTxPERC.ppi,PPI_IFGtxPERC.ppi,'r.');

    % BEST FIT LINES: COMP
    %----------------------------------------------------------------------
    x = PPI_vOTxCOMP.ppi(:);
    x = [x, ones(size(x))];
    y = PPI_IFGtxCOMP.ppi(:);
    B = x\y;
    y1 = B(1)*x(:,1)+B(2);
    plot(x(:,1),y1,'k-');

    % BEST FIT LINES: PERC
    %----------------------------------------------------------------------
    x = PPI_vOTxPERC.ppi(:);
    x = [x, ones(size(x))];
    y = PPI_IFGtxPERC.ppi(:);
    B = x\y;
    y1 = B(1)*x(:,1)+B(2);
    plot(x(:,1),y1,'r-');

    legend('COMP','PERC')
    xlabel('vOT activity')
    ylabel('IFGt response')
    title('Psychophysiologic Interaction')
%%
end
