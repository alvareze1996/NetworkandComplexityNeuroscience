%% Seperate correlation matricies for 50 subjects

sublist = [conn_parcellation_15ROI.subjects];
sublist =cellstr(sublist);

sublist = append('s',sublist);


for i=1:size(conn_parcellation_15ROI.conn_matrix_strength_Fishers_Z,3)
    equation = [sublist{i} '=conn_parcellation_15ROI.conn_matrix_strength_Fishers_Z  (:,:,' num2str(i) ')'];
    eval(equation);
    eval([sublist{i} '(isnan('  sublist{i} '))=1'])
    eval(['save ' sublist{i} '.mat' ' ' sublist{i}]);
end 
 


    