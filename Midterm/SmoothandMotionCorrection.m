%% Batch to move fMRI files and motion regressors from Box fMRIprep Format folder 
% This script works for all of the RiskAm subjects in via Box 
% Smoothed and Motion .tsv (6mm kernel) files in ea. Subjects' 'func' folder 
% Adapted from Alvarez Optimism Bias smoothing 

clear;
close all;

username = getenv('username');

% Path for folders 
SourceBox_data_path =  ['C:\Users\' username '\Box\RUBICMRI']; %source fmri folder RAW

Box_data_path = ['C:\Users\' username '\Box\Alvarez_Emmanuel\01_Projects\RestingState\Data\fmri']; %new folder copying smoothed files to

%% 
runNum = 4; %4 runs total w/ 30 trials each
smKernel = 6; %6mm smoothing kernel

% In case you need to clean files in the folder for future use!
% %Find the folder contents, i.e. Subjects
% d = dir(SourceBox_data_path);
% %Remove the files
% dFolders = d([d(:).isdir] == 1);
% %Get rid of extra folders
% dFolders = dFolders(~ismember({dFolders(:).name}, {'.','..' , 'logs'}));

% just look at folders, not .html files!
d = dir(SourceBox_data_path);
d = d(~startsWith({d.name}, '.'));
dFolders = d([d(:).isdir] == 1);

%% 

for s = 1:length(dFolders)
    wo_fs = dir([SourceBox_data_path filesep dFolders(s).name filesep 'sub-*']);
    w_fs = dir([SourceBox_data_path filesep dFolders(s).name filesep 'fmriprep']);
    if ~isempty(w_fs)
        close;
        subject_folder = [SourceBox_data_path filesep dFolders(s).name filesep 'fmriprep'];
        subject_ID = dir([subject_folder filesep 'sub-*']);
        subfolder = subject_ID.name;
        subject_ID = subject_ID([subject_ID(:).isdir] == 1);
        subject_ID = dir([subject_ID.folder filesep subject_ID.name]);
        subject_ID = dir([subject_ID(1).folder filesep 'ses-*']);
        subject_ID = subject_ID.name(5:end);
        
        %% Copy the motion correction and zipped Nii files obtained from Box fMRIprep locale into Box Optbias folder
        
        Box_fmriprep_func_path = [subject_folder filesep  subfolder filesep 'ses-' subject_ID filesep 'func']; %subject fmri data from source /getting files from
        cd(Box_fmriprep_func_path)
        
        Box_sub_func_path = [Box_data_path filesep dFolders(s).name filesep 'Processed_data']; %where we are moving files to
        if ~isfolder(Box_sub_func_path); mkdir(Box_sub_func_path); end
        
        
        for run = 1:runNum %4 total
            
            fMRI_zip = [subfolder '_ses-' subject_ID '_Rest_run-' num2str(run) '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'];
            motion_file = [subfolder '_ses-' subject_ID '_Rest_run-' num2str(run) '_desc-confounds_timeseries.tsv'];
            
            if exist([Box_fmriprep_func_path filesep fMRI_zip])~=0
                copyfile([Box_fmriprep_func_path filesep fMRI_zip], Box_sub_func_path) %Copy files to a separate folder for smoothing.
                copyfile([Box_fmriprep_func_path filesep motion_file], Box_sub_func_path)
                
                
                %% Unzip the downloaded files via gunzip(fMRI_zip_Nii) for the copied files
                cd(Box_sub_func_path)
                gunzip('*.gz')
                fMRI_file = extractBefore(fMRI_zip,'.gz');
                
                %% Smoothing of functional files (6mm kernel)
                
                spm('defaults','fmri');
                spm_jobman('initcfg');
                smooth = struct;
                
                %Load all volumes of fMRI file
                f4D_spm = spm_vol(fMRI_file);
                spm_size = size(f4D_spm);
                Nt = spm_size(1); %number of volumes
                fns={};
                for i = 1:Nt
                    fns{i} = [fMRI_file ',' num2str(i) ];
                end
                smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
                smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [smKernel smKernel smKernel];
                smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
                smooth.matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(smKernel) ];
                
                % Run smoothing with kernel 'sm' set at begining of this batch
                spm_jobman('run',smooth.matlabbatch);
                
                fMRI_s6_file = ['s' num2str(smKernel) fMRI_file];
                
                %% Get motion correction regressor file from fMRIprep files
                motcorrRun = tdfread(motion_file);
                motcorrRun = [motcorrRun.trans_x motcorrRun.trans_y motcorrRun.trans_z motcorrRun.rot_x motcorrRun.rot_y motcorrRun.rot_z];
                writematrix(motcorrRun, ['motionCorr_run' num2str(run) '.txt'],'Delimiter','tab')
                
                %% Move smoothed fMRI files and motion correction matrices into smoothed folder
                
                Subj_s6_path =[Box_sub_func_path filesep 'smooth_6mm'];
                mkdir(Subj_s6_path);
                movefile(fMRI_s6_file , Subj_s6_path)
                movefile(['motionCorr_run' num2str(run) '.txt'], Subj_s6_path)
                
                %% Housekeeping: delete zipped/unsmoothed fMRI and motion correction files
                delete(fMRI_zip)
                delete(fMRI_file)
                delete(motion_file)
            else
                display([[Box_fmriprep_func_path filesep fMRI_zip] ' does not exist']);
                continue;
            end
            
        end
        
    elseif ~isempty(wo_fs)
        close;
        subject_folder = [SourceBox_data_path filesep dFolders(s).name];
        subject_ID = dir([subject_folder filesep 'sub-*']);
            subfolder = subject_ID.name;
            subject_ID = subject_ID([subject_ID(:).isdir] == 1);
        subject_ID = dir([subject_ID.folder filesep subject_ID.name filesep 'ses-*'])
            subject_ID = subject_ID.name(5:end);
        
        %% Copy the motion correction and zipped Nii files obtained from Box fMRIprep locale into Box Optbias folder
        
        Box_fmriprep_func_path = [subject_folder filesep  subfolder filesep 'ses-' subject_ID filesep 'func']; %subject fmri data from source /getting files from
        cd(Box_fmriprep_func_path)
        
        Box_sub_func_path = [Box_data_path filesep dFolders(s).name filesep 'Processed_data']; %where we are moving files to
        if ~isfolder(Box_sub_func_path); mkdir(Box_sub_func_path); end
        
        
        for run = 1%:runNum % resting state run total
            
            fMRI_zip = [subfolder '_ses-' subject_ID '_task-OptBias_run-' num2str(run) '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'];
            motion_file = [subfolder '_ses-' subject_ID '_task-OptBias_run-' num2str(run) '_desc-confounds_timeseries.tsv'];
            
            if exist([Box_fmriprep_func_path filesep fMRI_zip])~=0
                copyfile([Box_fmriprep_func_path filesep fMRI_zip], Box_sub_func_path) %Copy files to a separate folder for smoothing.
                copyfile([Box_fmriprep_func_path filesep motion_file], Box_sub_func_path)
                
                
                %% Unzip the downloaded files via gunzip(fMRI_zip_Nii) for the copied files
                cd(Box_sub_func_path)
                gunzip('*.gz')
                fMRI_file = extractBefore(fMRI_zip,'.gz');
                
                %% Smoothing of functional files (6mm kernel)
                
                spm('defaults','fmri');
                spm_jobman('initcfg');
                smooth = struct;
                
                %Load all volumes of fMRI file
                f4D_spm = spm_vol(fMRI_file);
                spm_size = size(f4D_spm);
                Nt = spm_size(1); %number of volumes
                fns={};
                for i = 1:Nt
                    fns{i} = [fMRI_file ',' num2str(i) ];
                end
                smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
                smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [smKernel smKernel smKernel];
                smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
                smooth.matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(smKernel) ];
                
                % Run smoothing with kernel 'sm' set at begining of this batch
                spm_jobman('run',smooth.matlabbatch);
                
                fMRI_s6_file = ['s' num2str(smKernel) fMRI_file];
                
                %% Get motion correction regressor file from fMRIprep files
                motcorrRun = tdfread(motion_file);
                motcorrRun = [motcorrRun.trans_x motcorrRun.trans_y motcorrRun.trans_z motcorrRun.rot_x motcorrRun.rot_y motcorrRun.rot_z];
                writematrix(motcorrRun, ['motionCorr_run' num2str(run) '.txt'],'Delimiter','tab')
                
                %% Move smoothed fMRI files and motion correction matrices into smoothed folder
                
                Subj_s6_path =[Box_sub_func_path filesep 'smooth_6mm'];
                mkdir(Subj_s6_path);
                movefile(fMRI_s6_file , Subj_s6_path)
                movefile(['motionCorr_run' num2str(run) '.txt'], Subj_s6_path)
                
                %% Housekeeping: delete zipped/unsmoothed fMRI and motion correction files
                delete(fMRI_zip)
                delete(fMRI_file)
                delete(motion_file)
            else
                display([[Box_fmriprep_func_path filesep fMRI_zip] ' does not exist']);
                continue;
            end
            
        end       
        
    else 
        continue;
        display([dFolders(s).name  ' does not exist'])
        
    end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        

    
