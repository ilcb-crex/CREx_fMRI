function CREx_fMRI_DoRealignFigure
    % Do Realign Figure from rp* file on SPM12
    % Author: Valérie Chanoine, Research Engineer at Brain and Language
    % Institute (http://www.blri.fr/)
    % Date: July 25, 2016

    clear all;
    set(0,'DefaultFigureColor',[1 1 1]); % Set color background : blank

    st.subjects    = {'01LA', '02CE'}; 
    

    st.rootDir     = 'F:\IRM\Intermod2\Data\AUDIO';   
    st.sessions     = {'func01', 'func02', 'func03', 'func04'};    

    
    for ix=1:numel(st.subjects)  
        for iS=1:numel(st.sessions)  
            path = fullfile( st.rootDir,  st.subjects{ix}, 'Functional', st.sessions{iS});
   
            rp = spm_load(spm_select('FPList', path, '^rp.*\.txt$'));
        
            figure;

            subplot(2,1,1);plot(rp(:,1:3));
            set(gca,'xlim',[0 size(rp,1)+1]);
            title ( [st.subjects{ix} ' translation']);
            legend('x translation', 'y translation', 'z translation');
            xlabel('image');
            ylabel('mm')

            subplot(2,1,2);plot(rp(:,4:6));
            set(gca,'xlim',[0 size(rp,1)+1]);
            title([st.subjects{ix} ' rotation']);
            legend('pitch', 'roll', 'yaw');
            xlabel('image');
            ylabel('degrees');


            namfig = fullfile(path, 'RealignFigures');
            saveas(gcf, namfig, 'pdf'); 
            close;
        end
	end
end