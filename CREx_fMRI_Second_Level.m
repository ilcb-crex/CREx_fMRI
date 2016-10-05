function CREx_fMRI_Second_Level
    % Example of Second-Level implementation for SPM12
    % Author: Valérie Chanoine, Research Engineer at Brain and Language
    % Institute (http://www.blri.fr/)
    % Co-authors from BLRI: Samuel Planton and Chotiga Pattadimalok
    % Co-authors from fMRI platform:   Julien Sein, Jean-Luc Anton, Bruno Nazarian and Pascal Belin from fMRI
    % platform (fMRI Center, Timone Hospital, Marseille, France)
    % Date: Oct 10, 2016
    
    close all; clear all; clc;
      
    %% Initialise SPM
    spm('defaults','fmri');  
    spm_jobman('initcfg');   

    %======================================================================
    %  Define a general working strusture 'w'
    %======================================================================
    w.dataDir           = 'F:\IRM\Intermod2\Data\LOCA';    
    w.preprocDir        = 'PREPROC'; 
    w.firstDir          = 'FIRST_LEVEL';
    w.secondDir         = 'SECOND_LEVEL';
    
%     w.subjects          = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK',...
%         '10TL', '11CL', '12SP', '13CC', '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'};  % subject directory (parent=dataDir)
%      
    w.subjects = {'01LA', '02CE', '04AC', '05GV', '06EC', '08NB', '09BK','10TL', '11CL', '12SP',...
    '13CC', '14EP','17NR','18EF', '19AK', '20ET', '21LP','22SR', '23LG', '24IA', 'pilote2'};  % 21 subjects

     
        %==================================================================
        %  DO SECOND-LEVEL ANALYSIS
        %==================================================================  

        bEvent = false;
        bWithArt = false;
    
        DoCreateGroupMask(w);
        DoSecondLevel(w, bEvent, bWithArt); 
 
end         
 
function DoCreateGroupMask(w)

    %==================================================================
    %  Create binary mask for the subject group
    %==================================================================    
    outputDir = fullfile(w.dataDir, w.secondDir, 'ExplicitMask');

    if ~exist(outputDir,  'dir')
      mkdir (outputDir);        
    end

    %%
    fileMask = {};
    express=[]; 
    N = numel(w.subjects);
    for iS=1:N
   
        funcPath        = fullfile (w.dataDir, w.preprocDir, w.subjects{iS}, 'Structural', 'anat01');
        fileMask {iS}   = spm_select('FPList',  fullfile(funcPath), 'explicitMask_wc1wc2wc3_03.nii');
    
        % For Expression, add 'i1+iN' for avering the subject-specific masks  
        if (iS== N)
            express = [express 'i' num2str(iS)];
        else
            express = [express 'i' num2str(iS) '+'];
        end    
    end
    express=['((' express ')/' num2str(N) ')>.50'] ; % use  a threshold of 50%
    %%

    % Batch SPM
    clear matlabbatch;
    matlabbatch{1}.spm.util.imcalc.input =  fileMask';
    matlabbatch{1}.spm.util.imcalc.output = 'ExplicitMask.nii';
    matlabbatch{1}.spm.util.imcalc.outdir = {outputDir};
    matlabbatch{1}.spm.util.imcalc.expression = express;
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
end

function DoSecondLevel(w, bEvent,bWithArt)  

    %==============================================================
    %  2nd Level Analysis
    %==============================================================  

    w.ncond = 3;
    w.SecondLevelExplicitMask = fullfile(w.dataDir, w.secondDir, 'ExplicitMask','ExplicitMask.nii,1');
       
    if(bEvent)    
        w.FirstLevelDir  = fullfile(w.dataDir, w.firstDir, 'LOCA','FirstLevel_Event');
        w.SecondLevelDir = fullfile(w.dataDir, w.secondDir, 'LOCA','ANOVA_W_Event');
        if bWithArt
            w.FirstLevelDir  = fullfile(w.dataDir, w.firstDir, 'LOCA','FirstLevel_Event_Art');
            w.SecondLevelDir = fullfile(w.dataDir, w.secondDir, 'LOCA','ANOVA_W_Event_Art');
        end
    else
        w.FirstLevelDir  = fullfile(w.dataDir, w.firstDir, 'LOCA','FirstLevel_Block');
        w.SecondLevelDir = fullfile(w.dataDir, w.secondDir, 'LOCA','ANOVA_W_Block');
        if bWithArt
            w.FirstLevelDir  = fullfile(w.dataDir, w.firstDir, 'LOCA','FirstLevel_Block_Art');
            w.SecondLevelDir = fullfile(w.dataDir, w.secondDir, 'LOCA','ANOVA_W_Block_Art');
        end
    end
    
    %==============================================================
    %  fMRI model specification
    %==============================================================  
    
    clear matlabbatch;
    matlabbatch{1}.spm.stats.factorial_design.dir = {w.SecondLevelDir};

    for i=1:numel(w.subjects)
        %%
        matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(i).scans = {fullfile(w.FirstLevelDir,w.subjects{i},'con_0001.nii,1')
                                                                                  fullfile(w.FirstLevelDir,w.subjects{i},'con_0002.nii,1')
                                                                                  fullfile(w.FirstLevelDir,w.subjects{i},'con_0003.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(i).conds = [1 2 3];
    end

    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0; % Implicit Mask = 1
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {w.SecondLevelExplicitMask};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    %==============================================================
    %  Model Estimation
    %==============================================================   

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %==============================================================
    %  Contrast manager
    %==============================================================   

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

    % Contrasts T
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'WORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NONWORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'WORDS&NONWORDS > FIXATION';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 0.5 0.5];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'WORDS > NONWORDS';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{5}.fcon.name = 'F-all_conditions';
    matlabbatch{3}.spm.stats.con.consess{5}.fcon.weights =  eye(w.ncond)*(w.ncond-1)/w.ncond + (ones(w.ncond)-eye(w.ncond))*-1/w.ncond;
    matlabbatch{3}.spm.stats.con.consess{5}.fcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.delete = 1;

    %==============================================================
    %  Save & Run Batch
    %==============================================================   

    mkdir(char(w.SecondLevelDir))
    save([char(w.SecondLevelDir) '\SPM12_matlabbatch_10_2ndLevel_ANOVA_W.mat'],'matlabbatch');

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch)
end 