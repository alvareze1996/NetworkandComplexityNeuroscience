%%
%%
close all; 
clear all; 
clc; 
%% Paths 

username = getenv('username'); % To switch between desktop and laptop

%%% Toolboxes paths
addpath('C:\workbench\bin_windows64\')
addpath('C:\Users\ea472\Documents\gifti-main\')
addpath('C:\Users\ea472\Documents\cifti-matlab-master')
cd('C:\Users\ea472\Documents\ColeAnticevicNetPartition-master\ColeAnticevicNetPartition-master')

%%% Project folder 

   CIFTI_Data = ['C:\Users\' username '\Box\Preprocessed files\fMRIPrep\w_freesurfer_w_cifti'];

%%% Subject list
    CIFTIFiles = dir(CIFTI_Data);
    CIFTIFiles = CIFTIFiles(~startsWith({CIFTIFiles.name}, '.'));   
    
%%
conn_strength = struct;
for sub=1:length(CIFTIFiles)
    SubjectName =  CIFTIFiles(sub).name;
    Session = ['ses-' extractAfter(CIFTIFiles(sub).name,8)];
    restingState =['sub-'  SubjectName '_' Session '_task-rest_space-fsLR_den-91k_bold.dtseries.nii'];
    ScanFile = [CIFTI_Data filesep SubjectName filesep 'fmriprep\sub-' SubjectName filesep Session filesep 'func' filesep restingState];
    SubID = extractBefore(CIFTIFiles(sub).name,9);
    %%
    if exist(ScanFile)
        
    %%  Cole-parcellation
        %Setting the parcel files to be the 718 parcels (cortical + subcortical)
        parcelCIFTIFile='CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_ReorderedByNetworks.dlabel.nii';
        parcelTSFilename=['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RestingState\Analyses\CIFTI\' SubjectName '_Atlas_CortSubcort.Parcels.LR.ptseries.nii'];

        %Set this to be your input fMRI data CIFTI file
        inputFile = ScanFile;
        
        %%% Parcellate file
        eval(['!wb_command -cifti-parcellate ' '"' inputFile '"' ' ' parcelCIFTIFile ' COLUMN ' parcelTSFilename ' -method MEAN'])

        %Load parcellated data (requires the ciftiopen function from the HCP website, FieldTrip)
        LR_dat = ciftiopen(parcelTSFilename,'wb_command');

        NUMPARCELS=718;
        tseriesMatSubj=LR_dat.cdata;

        %Loading other relevant files
        load('cortex_subcortex_community_order.mat');
        netorder=readtable('network_labelfile.txt','ReadVariableNames',false);
        netassignments=table2array(readtable('cortex_subcortex_parcel_network_assignments.txt','ReadVariableNames',false));

  %% Computing Pearson correlation-based functional connectivity 
        
        [FCmat, P] =corrcoef(tseriesMatSubj');
        FCmat_sorted=FCmat(indsort,indsort);

        %%% FDR corrected p-values
        [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(P);
        
        %%% Fisher's z conversion
        z = 0.5.*(log(1+FCmat_sorted) - log(1-FCmat_sorted));
        
      %% Save matrix
      
      conn_strength.coeff(:,:,sub) = FCmat_sorted;
        FCmat_sorted(logical(eye(size(FCmat_sorted)))) = NaN; % remove identity line
        conn_strength.coeff_noIDline(:,:,sub) = FCmat_sorted;
      conn_strength.Pval(:,:,sub) = P;
        P(logical(eye(size(P)))) = NaN; % remove identity line
        conn_strength.Pval_noIDline(:,:,sub) = P;
      conn_strength.z(:,:,sub)= z;
        z(logical(eye(size(z)))) = NaN; % remove identity line
        conn_strength.z_noIDline(:,:,sub) = z;
      conn_strength.adjP(:,:,sub) = adj_p;
      conn_strength.crit_P(:,:,sub) = crit_p;
      
        %% Save .mat correlation matrix
%         folder = ['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RestingState\Analyses\Correlations'];
% 
%         fileName = [folder filesep SubID];
% 
%         eval(sprintf('save %s.mat z',fileName))


    else 
        %% if no resting state - skip
        continue;
    end 
end 

  

 