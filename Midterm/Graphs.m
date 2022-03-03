%% Subset of subjects (n=14)

deg=[]; % collect degrees of nodes that are above threshold
figure;
    for i=1:14
        nodes = CorrMatricies(15).pVal(:,:,i)<=0.05;
        G = graph(nodes,'lower') ;
           subplot(2,7,i)
           plot(G)
           temp = degree(G);
           deg = [deg temp];    
           
    end

MostDegrees = mean(deg,2);


