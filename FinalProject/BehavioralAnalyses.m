%% Behavioral figures/statistics 
username = getenv('username'); % To switch between desktop and laptop

%%% Clinical Data Path
RiskAmbDataPath=['C:\Users\ea472\Box\Alvarez_Emmanuel\01_Projects\RiskAmbig'];
load([RiskAmbDataPath filesep 'Data_RiskAmb.mat'])
TableDataPath = ['C:\Users\' username '\Box\RiskAmbiguity\01_a_Data'];

%% Demographics 
load([TableDataPath filesep 'Clinical' filesep 'RiskAmb_table1.mat'])
Sublist=unique(Data_RiskAmb.SubID);
Dx = RiskAmb_table.Group;
    CON_index = strcmp(Dx,'CONTROL');
    OUD_index = strcmp(Dx,'OUD');

%% GLMEs 

%%% Model- free
    %%% How do p, A and v of the lottery influence the choice made by the participants between the safe and lottery options, and do we observe group differences between CTR and OUD?
    glme = fitglme(Data_RiskAmb,'choice_w_ref ~ 1 +  prob*val*Group +  ambig*val*Group +(1|SubID)','Distribution','Binomial','Link','logit','FitMethod','Laplace'); %logistic regression
    disp(glme)


%%% Model- based 
    %%% How do subjective values drive choices, and do we observe group difference between CTR and OUD?
        %%% Overall 
        glme = fitglme(Data_RiskAmb,'choice_w_ref ~ 1 +  SV*Group +(1|SubID)','Distribution','Binomial','Link','logit','FitMethod','Laplace'); %logistic regression
        disp(glme)


        %%% Trial type interaction (BETTER FIT FOR MODEL-BASED ANALYSES)
        glme = fitglme(Data_RiskAmb,'choice_w_ref ~ 1 +  SV*Group*trialtype +(1|SubID)','Distribution','Binomial','Link','logit','FitMethod','Laplace'); %logistic regression
        disp(glme)
        
               %%% Broken by trials --- DO NOT USE 
              glme = fitglme(Data_RiskAmb(strcmp(Data_RiskAmb.trialtype,'Ambig'),:),'choice_w_ref ~ 1 +  SV*Group +(1|SubID)','Distribution','Binomial','Link','logit','FitMethod','Laplace'); %logistic regression
              disp(glme)
              
              glme = fitglme(Data_RiskAmb(strcmp(Data_RiskAmb.trialtype,'Risk'),:),'choice_w_ref ~ 1 +  SV*Group +(1|SubID)','Distribution','Binomial','Link','logit','FitMethod','Laplace'); %logistic regression
              disp(glme)
 
 %% Figures 
     oudcolor = [1.0000 0.6250 0.4766];
     controlcolor = [0.4375 0.5000 0.5625];  
 
 %%% Risk Tolerance 
 alpha=[];beta=[];
 for sub=1:length(Sublist)
     alpha(sub,1) = nanmean(Data_RiskAmb.alpha_softmax(Data_RiskAmb.SubID==Sublist(sub)));
     beta(sub,1) = nanmean(Data_RiskAmb.beta_softmax(Data_RiskAmb.SubID==Sublist(sub)));
 end 
 
   %%% Known risk          
      figure;
      hold on
      %%% Averages
        data =  [nanmean(alpha(CON_index)); nanmean(alpha(OUD_index))];

        b = bar(data,'FaceColor','flat');
          b.CData(1,:) = controlcolor;
          b.CData(2,:) = oudcolor;
          
         %%% Hold scatter
         ydata1 = alpha(CON_index);    
         ydata2 = alpha(OUD_index);
         
           
    xdata1= repmat([1],size(ydata1,1), 1);
    xdata2= repmat([2],size(ydata2,1), 1);
    
    scatter(xdata1(:), ydata1(:),25, 'o', 'jitter','on', 'jitterAmount', 0.01,'MarkerFaceColor',controlcolor,'MarkerEdgeColor','k' )
    scatter(xdata2(:), ydata2(:),25, 'o', 'jitter','on', 'jitterAmount', 0.01,'MarkerFaceColor',oudcolor,'MarkerEdgeColor','k' )
  
    %%% House keeping
        xlim([0 3])
        xticks([1 2])
        set(gca,'xticklabels',{'Controls', 'OUD'},'FontWeight','bold')
        %xticklabels({'Controls', 'OUD'},'FontWeight','bold')
        xtickangle(45)
        ylabel({'Risk tolerance (\alpha)'})
        %ylim([-3 3])
        set(gca,'tickdir','out')
        axis square
 
    %%% Unknown risk 
    
    
         figure;
      hold on
      %%% Averages
        data =  [nanmean(beta(CON_index)); nanmean(beta(OUD_index))];

        b = bar(data,'FaceColor','flat');
          b.CData(1,:) = controlcolor;
          b.CData(2,:) = oudcolor;
          
         %%% Hold scatter
         ydata1 = beta(CON_index);    
         ydata2 = beta(OUD_index);
         
           
    xdata1= repmat([1],size(ydata1,1), 1);
    xdata2= repmat([2],size(ydata2,1), 1);
    
    scatter(xdata1(:), ydata1(:),25, 'o', 'jitter','on', 'jitterAmount', 0.01,'MarkerFaceColor',controlcolor,'MarkerEdgeColor','k' )
    scatter(xdata2(:), ydata2(:),25, 'o', 'jitter','on', 'jitterAmount', 0.01,'MarkerFaceColor',oudcolor,'MarkerEdgeColor','k' )
  
    %%% House keeping
        xlim([0 3])
        xticks([1 2])
        set(gca,'xticklabels',{'Controls', 'OUD'},'FontWeight','bold')
        xtickangle(45)
        ylabel({'Ambiguity tolerance (\beta)'})
        %ylim([-3 3])
        set(gca,'tickdir','out')
        axis square
 
