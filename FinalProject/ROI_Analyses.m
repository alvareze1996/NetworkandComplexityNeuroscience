%% RISK AMBIGUITY ROI extraction 
clear; close all;



%% Paths
project = 'RiskAmbiguity';
ADNlab_path = 'C:\Users\ea472\Box\';
ADNlab_project_path = [ADNlab_path filesep project];
ADNlab_data_path = [ADNlab_project_path filesep '01_a_Data'];
ADNlab_behavior_path = [ADNlab_data_path filesep 'Behavior'];
ADNlab_fMRI_path = [ADNlab_path filesep 'RiskGLM'];
ADNlab_clinical_path = [ADNlab_data_path filesep 'Clinical'];
ADNlab_MASTERclinical_path = [ADNlab_path filesep '09_Master_Clinical_Data'];
ADNlab_script_path = [ADNlab_project_path filesep '01_b_Analysis_and_Scripts'];


%% Load model-free and model-based data
load([ADNlab_behavior_path filesep 'Choices_allsub.mat']); %stick to group average parameters for figures
load([ADNlab_behavior_path filesep 'group_modelfit_parameters_fmincon.mat']); %all subjects with behavior
modeling_table=table(...
    group_modelfit_parameters.sub, ...
    group_modelfit_parameters.alpha_softmax, ...
    group_modelfit_parameters.beta_softmax, ...
    group_modelfit_parameters.gamma_softmax, ...
    group_modelfit_parameters.AIC_softmax, ...
    group_modelfit_parameters.BIC_softmax, ...
    group_modelfit_parameters.LL_softmax, ...
    group_modelfit_parameters.r2_softmax, ...
    'VariableNames',{'SubID','alpha_softmax','beta_softmax','gamma_softmax','AIC_softmax','BIC_softmax','r2_softmax','LL_softmax'});
RiskAmb_MFMBdata = join(Choices_allsub,modeling_table);

%% Load clinical data and demographics

TableDataPath = ['C:\Users\ea472\Box\RiskAmbiguity\01_a_Data'];
load([TableDataPath filesep 'Clinical' filesep 'RiskAmb_table1.mat'])

%% List of regions of interest
ROI_list={'vmPFC','VS_L','VS_R','midbrain','insula_L','insula_R','amygdala_L','amygdala_R','dACC','dlPFC_L','dlPFC_R','lOFC_L','lOFC_R','IFG_L','IFG_R','PCC','vmPFC_riskyVSambig','leftIFG_forManny'};%,'thalamus_L','thalamus_R'
MNI_allROI = zeros(length(ROI_list),4); MNI_allROI(MNI_allROI==0)=NaN;

%% Load template SPM for AnnieROI.m
%load([ADNlab_script_path filesep 'Template_SPM.mat']); %template SPM.mat file
load([ADNlab_fMRI_path filesep 'Second_level_GLM' filesep 'RiskAmb_Decision_AmbigMF_LotteryAmbig_DX' filesep 'SPM.mat'])

%% Select maps for ROI extraction (contrast maps)
RiskAmb_contrasts={
    'RiskAmb_Decision_LotteryAmount', ...
    'RiskAmb_Decision_LotteryProba', ...
    };
ROImaps_path = [ADNlab_fMRI_path filesep 'Contrast_maps']; % directory where you have all subjects' contrast maps

%% Extract ROIs 

for contrast_nb = 1:length(RiskAmb_contrasts)
    cd(ROImaps_path)
    con_name = char(RiskAmb_contrasts(contrast_nb)); %name of the 1st level GLM contrasts to study maps for.
    dsubmaps = dir([ROImaps_path filesep '*_' con_name '*.nii']); %list of all maps for a given contrast (several maps per subject because multiple levels)

    ROI_extracted_mean_allsub=[];
    for sub=1:length(dsubmaps)
        ROI_extracted_mean_con=[];
        imagecon=dsubmaps(sub).name;
        subject=str2num(extractBefore(imagecon,'_RiskAmb'));
        
               
        for ROI=1:length(ROI_list)
            ROI_name = char(ROI_list(ROI));
            switch ROI
                %%% value regions
                case 1 %vmPFC
                    MNI = [-1 46 -7];
                case 2 %ventral striatum left
                    MNI = [-10 10 -4];
                case 3 %ventral striatum right
                    MNI = [10 10 -4];
                case 4 %midbrain
                    MNI = [0 -20 -20];
                    
                %%% risk/error regions
                case 5 %insula left
                    MNI = [-38 16 -7];
                case 6 %insula right
                    MNI = [38 16 -7];
                case 7 %amygdala left
                    MNI = [-21 -2 -22];
                case 8 %amygdala right
                    MNI = [21 -2 -22];
                case 9 %dACC
                    MNI = [0 14 27];
                case 10 %dlPFC left
                    MNI = [-36 -4 60];
                case 11 %dlPFC right
                    MNI = [36 -4 60];
                case 12 %lOFC left
                    MNI = [-35 45 -10];
                case 13 %lOFC right
                    MNI = [35 45 -10];
%                 case 14 %thalamus left
%                     MNI = [-15 -15 6];
%                 case 15 %thalamus right
%                     MNI = [15 -15 6];
                case 14 %IFG left (Tali Sharot)
                    MNI = [-36 33 8];
                case 15 %IFG right (Tali Sharot)
                    MNI = [36 33 8];
                case 16 %PCC (Bartra NIMG2013)
                    MNI = [-4 -30 36];
                case 17 %vmPFC RiskAmb peak SV risk vs ambig
                    MNI = [-4 42 -13];
                case 18 %left IFG for Manny's NRSA grant Dec2021
                    MNI = [-48.5 17.7 18.3];
            end
            radius=5; %can be 5 or 10mm
            MNI_allROI(ROI,:) = [MNI radius];
            [extracted_mean_con] = AnnieROI(imagecon,SPM,MNI,radius);
            ROI_extracted_mean_con=[ROI_extracted_mean_con , extracted_mean_con]; %all ROI for 1 subject and 1 contrast
        end
%         ROI_extracted_mean_sub=[subject , sub_dx , ROI_extracted_mean_con, sub_sex , sub_age , sub_alpha , sub_beta, sub_gamma, sub_AIC, sub_BIC, sub_r2, sub_LL , sub_druguse_pastweek, sub_druguse_nextweek,sub_craving_pastweek,sub_craving_rightnow ,sub_craving_globalNow]; %include subID in first column
        ROI_extracted_mean_sub=[subject , sub_dx , ROI_extracted_mean_con, sub_sex , sub_age , sub_alpha , sub_beta, sub_gamma, sub_AIC, sub_BIC, sub_r2, sub_LL, sub_craving_pastweek, sub_druguse_pastweek]; %include subID in first column
        ROI_extracted_mean_allsub=[ROI_extracted_mean_allsub ; ROI_extracted_mean_sub]; %all ROI for all subjects and 1 contrast
    end
    
end




