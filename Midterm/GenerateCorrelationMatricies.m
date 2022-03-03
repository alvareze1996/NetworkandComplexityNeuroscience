% Script to run seed-based connectivity between seeds from fMRI Prep Data

clear all;
close all;

%% Parameters to set:

TR = 1; %TR in seconds
fmri_t = 16;
fmri_t0 = 1;

root_folder = ['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RestingState\Data'];
%batch_folder = ['/Users/annakonova/Documents/Matlab_codes/RESTING_seed_conn_analysis/GLM_batches'];
addpath('C:\Users\ea472\Box\Alvarez_Emmanuel\05_Resources\MRIResources\Matlabscripts_MRI');
load(['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RestingState\Scripts\conn_parcellation_15ROI.mat'])

parcell_coords=conn_parcellation_15ROI.parcell_coords(:,1:3);
radius = 10;

smoothing_kernel = 6;
forcedo_smoothing = 0; %to re-smooth existing smoothed data, set to 1
tpcount = 695;
motion_parameter_number = 6;

%% bandpass filter 0.009-0.08 Hz
passband = [0.01 0.1]; 
pb = passband./((1/TR)/2);
[b,a] = butter(8,pb); % this is to generate a Butterworth filter with parameters a and b

%% Model specification:

SubID = dir([root_folder filesep 'sub-*']);
SubID = SubID(~startsWith({SubID.name}, '.'));

CorrMatricies = struct; 
conn_strength = [];
for sub=1:length(SubID)
    
    clear img_files_smoothed_cell
    
    image_folder = [root_folder filesep SubID(sub).name];
   
        %% construct seeds & WM & csf & motion
        %seeds
        load(['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RestingState\Scripts\SPM_file_to_load\SPM.mat']); 
        img_files_smoothed = spm_select('ExtFPList',image_folder,['^s6' SubID(sub).name '.*\.nii' '$'],1:tpcount);

        for i=1:size(img_files_smoothed,1)
            if ~isempty(strfind(img_files_smoothed(i,:),' '))
                img_files_smoothed_cell{i,1} = img_files_smoothed(i,1:length(img_files_smoothed(i,:))-length(strfind(img_files_smoothed(i,:),' ')));
            else
                img_files_smoothed_cell{i,1} = img_files_smoothed(i,:);
            end
        end
        
        
        conn_matrix=[];
        for seed=1:size(parcell_coords,1)
            ind_seed=[];
            for i=1:size(img_files_smoothed_cell,1)
                [temp] = AnnieROI(img_files_smoothed_cell{i},SPM,parcell_coords(seed,1:3),radius);
                ind_seed = [ind_seed; temp];
            end
            ind_seed = filtfilt(b,a,ind_seed); % band-pass filtering
            conn_matrix = [conn_matrix; ind_seed'];
        end
        conn_matrix = conn_matrix';
        
             
        
        %load nuisance variables

        
        motionfile = dir([image_folder filesep 'motion_rest.txt']);
        motion_regs = dlmread([image_folder filesep motionfile.name]);
        
        
        %CSF_seed
        CSF_seed = filtfilt(b,a,motion_regs(:,8));
        %CSF_seed = filtfilt(b,a,nuisance_regs.csf); % band-pass filtering
        
        %WM_seed
        WM_seed = filtfilt(b,a,motion_regs(:,7)); % band-pass filtering
        
        %motion_regs   
        
        motion_regs = motion_regs(:,1:6);
            
       %% compute partial connectivity strengths: Pearson's r
    [RHO,pval] = partialcorr(conn_matrix,conn_matrix,[WM_seed CSF_seed motion_regs]);
    RHO(logical(eye(size(RHO)))) = NaN; % remove identity line
    pval(logical(eye(size(RHO)))) = NaN; % remove identity line
    conn_matrix_strength_RHO(:,:,sub) = RHO;
    conn_matrix_strength_PVal(:,:,sub) = pval;
    
    %% convert r to Fisher's Z to approximate normal distribution
    z = 0.5.*(log(1+RHO) - log(1-RHO));
    conn_matrix_strength_Fishers_Z(:,:,sub) = z;
    
    
    %% save in matrix
    CorrMatricies(sub).Matrix = conn_matrix_strength_Fishers_Z;
    CorrMatricies(sub).pVal = conn_matrix_strength_PVal;
    
    end
   
    
