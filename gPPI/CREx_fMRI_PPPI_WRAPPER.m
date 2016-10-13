function CREx_fMRI_PPPI_WRAPPER
%==========================================================================
%   GPPI_WRAPPER
%
%   This script allows you to run Donald McLaren's gPPI toolbox to build 
%   and estimate psychophysiological interaction (PPI) models for multiple
%   subjects and for multiple seed volumes-of-interest (VOI). 
%   You can download the toolbox and find further information at:
%   http://www.brainmap.wisc.edu/PPI
%   Author: Valérie Chanoine, Research Engineer at Brain and Language
%   Co-authors from BLRI: Samuel Planton and Chotiga Pattadimalok
%==========================================================================

    clear all; 

    %======================================================================
    % Settings in a "working" structure w
    %======================================================================   
    
    % Paths
    %----------------------------------------------------------------------
    w.gPPIpath        = fullfile(spm('dir'), 'toolbox', 'PPPI'); % path for gPPI toolbox
    w.dataPath        = fullfile('F:\IRM\Intermod2\Data');       % data path
    addpath(w.gPPIpath)
    
    cd(fullfile(w.dataPath, 'AUDIO', 'GLM_Event'));
    
	% Subjects 
    %----------------------------------------------------------------------   
    % w.subjects     	= {'01LA', '02CE'};
    w.subjects     	= {'02CE'};
   
    % Regions
    %----------------------------------------------------------------------  
    w.regionfile = {'VOI_vOT_1.mat', 'VOI_IFG_Tri_1.mat'};
    w.regionname = {'vOT', 'IFG_Tri'};   
    
 	% Conditions of level 1 GLM
    %----------------------------------------------------------------------    
    w.conditions  	= {'COMPR-SILENCE' 'PERC-SILENCE'}; 
    w.weights    	= [1 -1];  % Weights: for traditional PPI, you must specify weight vector for each condition.
    w.PPIweights    = [1 0 0 zeros(1, 4)];      % numSessions = 4;          
    
  	% Define the method generalized (condition-specific) method 
    %----------------------------------------------------------------------     
    % 'trad' for traditional method, 'cond' for generalized (condition-specific) method 
    w.method = 'cond';

    
    % Loop on Subjects
    %---------------------------------------------------------------------- 
    for iSub=numel(w.subjects)
        w.GLMPath = fullfile(w.dataPath, 'AUDIO', 'GLM_Event', w.subjects{iSub})

        % Loop on regions
        %----------------------------------------------------------------------   
        for iVOI=1:numel(w.regionfile)

            %==============================================================          
            % Define the parameter structure P
            %==============================================================
            P.subject   = w.subjects{iSub};
            P.directory = w.GLMPath; %either the first-level SPM.mat directory, or if you are only estimating a PPI model, then the first-level PPI directory.  
            P.method  	= w.method;   % 'trad' for traditional method, 'cond' for generalized (condition-specific) method  
            
            P.VOI       = fullfile(w.GLMPath, w.regionfile{iVOI});
                %% OR
                % %             P.VOI.VOI = [maskpath filesep voiNAMES{v}];
                % %             P.VOI.masks = {'mask.img'}; % name of images to use to constrain definition of the seed region
                % %             P.VOI.thresh = .5;  %values to threshold the images above at (number of values must match number of masks)
                % %             P.VOI.min = 10;  % minimum number of voxels to accept a seed region as valid  
                % %             P.VOImin = 10;  % minimum number of voxels to accept a seed region as valid     
                
               
            P.Region    =   w.regionname{iVOI};     % name of output file(s)
            %P.contrast  =  'F-TOUTES-CONDITIONS';   % contrast to adjust  (f constrast number in the level 1 GLM)
            P.contrast  = 1;
            P.analysis  = 'psy';                    % 'psy' for psychophysiological interaction
            P.extract   = 'eig';                    % method of ROI extraction, eigenvariate ('eig') or mean ('mean')
            P.equalroi  = 1;                        % the ROIs must be the same size in all subjects; set to 0 to lift the restriction
            P.FLmask    = 0 ;                       % specifies that the ROI should be restricted using the
                                        %               mask.img from the 1L statistics. NOTE: default=0.        
            %P.VOI2
      
                
            if strcmp(w.method,'trad')
                P.Tasks = w.conditions;  % condition names of level 1 GLM
                P.Weights = w.weights;      % Weights: for traditional PPI, you must specify weight vector for each task.
                P.Contrasts(1).name = 'PPI';
                P.Contrasts(1).left = w.PPIweights;
                P.Contrasts(1).right = [];
            else
                P.Tasks = ['0' w.conditions]; % condition names of level 1 GLM
            end             
                
              
            P.Estimate  = 1;      % estimate the PPI design (1 = yes, 0 = no)
            P.CompContrasts = 0; % estimate any contrasts (1 = yes, Default = no)

            
            if strcmp(w.method,'cond') % Specify  contrasts oonly for 'Cond PPI'             
               for iC = numel(w.conditions)
                   P.Contrasts(iC).name = w.conditions{iC};
                   P.Contrasts(iC).left = w.conditions{iC};
                   P.Contrasts(iC).right = {'none'}; 
               end
            end 
            
            %==============================================================          
            % Run gPPI Toolbox
            %==============================================================           
            PPPI(P);
            
        end
    end
end



