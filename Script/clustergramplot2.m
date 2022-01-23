function h2  = clustergramplot2( drug, teststats, results_ablation_testAUC, results_mean_AUCs_test, DrugOrder )

% fix the labels

combi2 = nchoosek( 1:6, 2 );
combi3 = nchoosek( 1:6, 3 );

combinations = cell( 41, 1 );

for j = 1:6
   
    combinations{ j, 1 } = j;
    
end

for j = 1:size( combi2, 1 )
   
    combinations{ j+6, 1 } = combi2( j, : );
    
end

for j = 1:size( combi3, 1 )
   
    combinations{ j+6+size( combi2, 1 ), 1 } = combi3( j, : );
    
end

prexvalues = { 'Mut, ', 'CNV, ', 'Meth, ', 'Tis, ', 'PWAct, ', 'Gex, ' };
xvalues = {};
index = 1;

for k = 1:size( combinations, 1 )
    
    temp = [ prexvalues{ combinations{ k, 1 } } ];
    xvalues{ index, 1 } = temp( 1:length( temp ) - 2 ); %#ok<*AGROW>
    index = index + 1;
    
end

yvalues = { 'Neural network', 'LASSO-regularised lin. regr.', 'E. net-regularised lin. regr.', 'Ridge-regularised lin. regr.','LASSO-regularised log. regr.', ...
'E. net-regularised log. regr.', 'Ridge-regularised log. regr.', 'LASSO-regularised SVM', 'Ridge-regularised SVM', 'Bagged DT ensemble', 'Boosted DT ensemble', ...
'Random forest' };

% calculate the index 

index = nansum( results_ablation_testAUC{ drug } ) > 0;

% first clustergram with absolute values

% temp = results_ablation_testAUC{ drug };
% temp = temp( :, index );
% temp = temp( [ 1:9,11:13 ], : );
% 
% h = clustergram( temp, 'Symmetric', 'false', 'Colormap', 'jet', 'ColumnLabels', xvalues( index ), 'RowLabels', yvalues, 'ColumnLabelsRotate', 45, 'Annotate', true  );
% 
% addTitle( h, [ 'Ablationstudies for ', DrugOrder{ drug, 1 }, ': absolute mean ROC-AUCs in testing for reduced models' ] );


% second clustergram with normalized values

rat = results_mean_AUCs_test( drug, 7:18);

temp = results_ablation_testAUC{ drug }';
temp = temp( index, : );
temp = temp( :, [ 1:9,11:13 ] );

temp = temp./rat;


%h2 = clustergram( temp', 'Symmetric', 'false', 'Colormap', 'jet', 'ColumnLabels', xvalues( index ), 'RowLabels', yvalues, 'ColumnLabelsRotate', 45, 'Annotate', true  );
h2 = clustergram( temp, 'Symmetric', 'false', 'Colormap', 'turbo', 'ColumnLabels', ...
    yvalues, 'RowLabels', xvalues( index ), 'ColumnLabelsRotate', 35, ...
    'RowLabelsRotate', 35, 'Annotate', true, 'AnnotColor', 'w'  ); % 'AnnotColor', [ 0.6, 0.6, 0.6 ]

%addTitle( h2, [ 'Ablationstudies for ', DrugOrder{ drug, 1 }, ': relative mean ROC-AUCs in testing for reduced models' ] );


end