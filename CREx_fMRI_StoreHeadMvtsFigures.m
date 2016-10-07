function CREx_fMRI_StoreHeadMvtsFigures
% Store RealignFigures.pdf file from subject folder 
% Author: Valérie Chanoine, Research Engineer at Brain and Language
% Research Institute (http://www.blri.fr/)
% Date: Oct 10, 2016

    clear all;

    subjects  = {'01LA', '02CE', '03NC', '04AC', '05GV', '06EC', '07MR', '08NB', '09BK',...
            '10TL', '11CL', '12SP', '13CC', '14EP','16SM','17NR','18EF', '19AK', '20ET', '21LP', '22SR', '23LG', '24IA', 'pilote2'}; 

    rootSrc     = 'F:\IRM\Intermod2\Data\PREPROC';     
    sessions    = {'func01', 'func02', 'func03', 'func04', 'func05'};

    bArtDetection =  true;
        
    for iSub=1:numel(subjects)

        if (bArtDetection)
            srcDir  = fullfile(rootSrc, subjects{iSub});
            src     = spm_select('FPList', srcDir, ['^' subjects{iSub} '_art_config.*\.png$']); 

            rootDest    = fullfile(rootSrc, 'HeadMvts_ART');
            if ~exist(rootDest, 'dir')
                mkdir(rootDest)
            end

            dest    = fullfile(rootDest, ['ART_Detection' subjects{iSub} '.png']);
            copyfile(src, dest);     

        else        

            for j=1:numel(sessions)
                srcDir  = fullfile(rootSrc, subjects{iSub}, 'Functional', sessions{j});
                src     = spm_select('FPList', srcDir, '^RealignFigures.*\.pdf$');  

                rootDest    = fullfile(rootSrc, 'HeadMvts_SPM');                
                if ~exist(rootDest, 'dir')
                     mkdir(rootDest)
                end

                dest    = fullfile(rootDest, ['HeadMvts_SPM_' subjects{iSub} '_' sessions{j} '.pdf']);
                copyfile(src, dest);
            end
        end
    end    
end