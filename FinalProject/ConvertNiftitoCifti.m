%% Add workbench to path in case not in environment
addpath('C:\workbench\bin_windows64\')
%% Get template volumes

ciftiseperate ='sub-21760004122018_ses-122018_task-RISK_run-1_space-fsLR_den-91k_bold.dtseries.nii';
volumeCifti = 'ciftiVolume.nii';
volumeLabel = 'labels.dlabel.nii';

eval(['!wb_command -cifti-separate ' ciftiseperate ' COLUMN ' ' -volume-all ' volumeCifti ' -label ' volumeLabel])

%% Volume files to convert into CIFTI
conmap = 'spmT_0001.nii'; % Contrast map
volumespace = volumeCifti;
resampled = 'resampledvolumne.nii';
eval(['!wb_command  -volume-resample ' conmap ' ' volumespace ' ' 'TRILINEAR' ' ' resampled])

%% Get metric files from surface data

leftSurface = 'sub-21760004122018_ses-122018_hemi-L_midthickness.surf.gii';

rightSurface ='sub-21760004122018_ses-122018_hemi-R_midthickness.surf.gii';

    leftMetric = 'left_metric.shape.gii';

    rightMetric = 'right_metric.shape.gii';

 eval(['!wb_command -volume-to-surface-mapping ' resampled ' ' leftSurface ' ' leftMetric ' -trilinear'])
 eval(['!wb_command -volume-to-surface-mapping ' resampled ' ' rightSurface ' ' rightMetric ' -trilinear'])

%% Surface spheres 

%%% Input spheres files from free surfer 
L_sphere = 'lh.sphere.reg';
R_sphere = 'rh.sphere.reg';

    %%% GIFTI spheres
    save(gifti(L_sphere), 'lh.surf.gii')
    save(gifti(R_sphere), 'rh.surf.gii')
   
    L_sphere = 'lh.surf.gii';
    R_sphere = 'rh.surf.gii';

%%% Failed creating an arbitrary sphere with size 
    %  eval(['!wb_command -surface-create-sphere 534431 ' L_sphere])
    %  eval(['!wb_command -surface-flip-lr ' L_sphere ' ' R_sphere])
    %  eval(['!wb_command -set-structure ' R_sphere ' CORTEX_RIGHT'])
    %  eval(['!wb_command -set-structure ' L_sphere ' CORTEX_LEFT'])


%% Convert NIFTI into cifti

%%% output cifti files
output='output.dscalar.nii';

%%% function to convert nifti to cifti
eval(['!wb_command -cifti-create-dense-scalar -volume ' resampled ' ' volumeLabel ' -left-metric ' leftMetric ' -right-metric ' rightMetric ' ' output])

%%% Created a cifti! However, is missing information for parcellation
%cifti_read(output)
%% Attemtped to resample as a way to help with parcellation

sampleCifti= 'CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_ReorderedByNetworks.dscalar.nii';
LsphereOut =  'S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
RsphereOut= 'S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
ciftiOUT = 'ciftiOutput.dscalar.nii';

eval(['!wb_command -cifti-resample ' output ' COLUMN ' sampleCifti ' COLUMN ' ' ADAP_BARY_AREA ' ' TRILINEAR ' ciftiOUT...
    ' -left-spheres ' L_sphere ' '    LsphereOut  ' -left-area-surfs ' leftSurface ' ' LsphereOut...
    ' -right-spheres ' R_sphere ' '   RsphereOut ' -right-area-surfs ' rightSurface ' '   RsphereOut])
    %%% ERROR: left current sphere doesn't match input cifti 

%% Additional gerry rigging

    %%% Create a template to fit size of con map
    x = cifti_read('sub-21760004122018_ses-122018_task-rest_space-fsLR_den-91k_bold.dtseries.nii')
    x.cdata = mean(x.cdata,2);
    x.diminfo{1, 2}.length  = 1;
    ciftisave(x,'testtemplate.dtseries.nii')
    ciftiTemplate = 'testtemplate.dtseries.nii';

    %%% Take conmap again 
    conmap = 'spmT_0001.nii';
    %%% output name
    outputCifti = 'conmap.dtseries.nii';

    eval(['!wb_command -cifti-convert -from-nifti ' conmap ' ' ciftiTemplate ' ' outputCifti])



