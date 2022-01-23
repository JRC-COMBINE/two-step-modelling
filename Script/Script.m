% Script for the Thesis

%% Run the two-step modelling framework on every drug in the GDSC data set

for drug = 1:265

    [ models1, models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part ] = twostepmodel( drug );

    parsave( strcat( "res_", num2str( drug ) ), models1, models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part );
    
end

%% Sort the results into arrays

% AUCs as measures of predictive performance

results_mean_AUCs_tr = zeros( 265, 18 );
results_mean_AUCs_test = zeros( 265, 18 );
results_std_AUCs_tr = zeros( 265, 18 );
results_std_AUCs_test = zeros( 265, 18 );

for count = 1:265

    str = strcat( 'res_', num2str( count ), '.mat' );
    
    load( str, 'AUCs' );
    
    results_mean_AUCs_tr( count, : ) = nanmean( cell2mat( AUCs( :, 1 ) ) ); %#ok<*NANMEAN> 
    results_mean_AUCs_test( count, : ) = nanmean( cell2mat( AUCs( :, 2 ) ) );
    
    results_std_AUCs_tr( count, : ) = nanstd( cell2mat( AUCs( :, 1 ) ) ); %#ok<*NANSTD> 
    results_std_AUCs_test( count, : ) = nanstd( cell2mat( AUCs( :, 2 ) ) );

end


% Test AUCs as measures of the predictive performance of ablation models

results_ablation_testAUC = cell( 1, 265 );

for count = 1:265
    
    resultmat = NaN( 13, 41 );

    str = strcat( 'res_', num2str( count ), '.mat' ) %#ok<*NOPTS>
    
    load( str, 'ablstudies_results' );
    
    temp = ablstudies_results( :, 3 );
    
    for l = 1:13 % run through all algorithms
    
        tempmatrix = NaN( 10, 41 );
        
        if l == 10 % Naïve Bayes does not have a ROC-AUC
            
            for k = 1:10 % run through all folds
                
                temptemp = temp{ k, 1 }{ l, 1 };
                tempmatrix( k, : ) = temptemp( 1, : );
                
            end
            
        else
        
            for k = 1:10 % run through all folds
                
                temptemp = temp{ k, 1 }{ l, 1 };
                tempmatrix( k, : ) = temptemp( 2, : );
                
            end
        
        end
        
        ind = sum( isnan( tempmatrix ) ) + sum( tempmatrix == 0 ) <= 5;
    
        resultmat( l, ind ) = nanmean( tempmatrix( :, ind ) );
        
    end
    
    results_ablation_testAUC{ 1, count } = resultmat;

end


%% Chapter 3 - Data

load( 'SomaticMutation.mat', 'SomaticMutation' );
load( 'CNV.mat', 'CNV' );
load( 'Hypermethylation.mat', 'Hypermethylation' );
load( 'TissueType.mat', 'TissueType' );
load( 'GeneExpression.mat', 'GeneExpression' );
load( 'PathwayActivation.mat', 'PathwayActivation' );
load( 'Response_AUC.mat', 'Response_AUC' );

load( 'TissueOrder.mat', 'TissueOrder' );
load( 'Celllines_info.mat', 'Celllinesinfo' );
load( 'DrugOrder.mat', 'DrugOrder' );

% Figure 3.1a

colours = [ 117, 146, 213; 189, 136, 213; 223, 197, 213; 188, 221, 213; 
            36, 145, 213; 126, 87, 126; 227, 125, 126; 32, 86, 126; 
            0, 114, 189; 203, 218, 250; 220, 218, 35; 229, 218, 171; 
            214, 218, 210; 138, 174, 79; 244, 174, 36; 80, 197, 247;
            229, 113, 46; 64, 64, 64; 45, 113, 219; 25, 74, 11;   
            236, 74, 14; 208, 74, 106; 163, 23, 24; 91, 23, 99;  
            248, 236, 63; 248, 206, 176; 216, 238, 148; 9, 42, 79;  
            109, 11, 34; 125, 55, 24 ];
colours = colours./255;

figure;

r = treemap( tissues, 1.6, 1 );
plotRectangles( r, TissueOrder( :, 2 ), colours )

% Figure 3.1b

cancertypes = unique( Celllinesinfo( 2:969, 10 ));

cancers = zeros( length( cancertypes ), 1 );

for k = 1:length( cancertypes )
    
    cancers( k, 1 ) = sum( strcmp( Celllinesinfo( 2:969, 10 ), cancertypes{ k, 1 } ) );
    
end

% 14 cell lines cannot be classified ('UNABLE TO CLASSIFY') and 171 cell lines have not been matched to a TCGA label; 
% therefore, they are to be excluded from the overall count

cancers = cancers( 2:31 );
cancertypes = cancertypes( 2:31 );

figure;

r = treemap( cancers, 1.6, 1 );
plotRectangles( r, cancertypes, colours )

% Figure 3.2

drugtypes = unique( DrugOrder( :, 4 ) ); 

drugnumbers = zeros( 21, 1 );

for k = 1:21
    
    drugnumbers( k, 1 ) = sum( cell2mat( DrugOrder( :, 5 ) ) == k );
    
end

colours2 = [ 117, 146, 213; 189, 136, 213; 223, 197, 213; 188, 221, 213; 
             36, 145, 213; 126, 87, 126; 227, 125, 126; 32, 86, 126; 
             28, 81, 17; 203, 218, 250; 220, 218, 35; 138, 174, 79;  
             244, 174, 36; 154, 174, 96; 229, 113, 46; 45, 113, 219; 
             217, 83, 25; 208, 74, 106; 91, 23, 99; 248, 236, 63; 
             109, 11, 34]./255;

figure;

r = treemap( drugnumbers, 1.6, 1 );
plotRectangles( r, drugtypes, colours2 )

% Figure 3.3a

mutsums = nansum( SomaticMutation ); %#ok<*NANSUM> 
cnvsums = nansum( abs( CNV ) );
methsums = nansum( Hypermethylation );

figure;
    
boxchart( [mutsums', cnvsums', methsums' ], 'JitterOutliers', 'on', 'MarkerStyle', '.' )

box on;
grid on;
ylim( [ -100, 700 ] )
yticks( 0:100:700 )

xticklabels( { 'Somatic mutations', 'CNVs', 'Hypermethylations' } )
ylabel( 'Number of cell lines featuring a molecular event' )

set( gca, 'FontSize', 14 )

% Figure 3.3b

temp = sum( ~isnan( Response_AUC ) );

figure;
histogram( temp, 40 )

hold on;

plot( [ 500, 500 ], [ 0, 90 ], 'k-.' )
plot( [ 600, 600 ], [ 0, 90 ], 'k-.' )

text( 490, 80.5, 'Sample size \leq 500 \leftarrow', 'FontSize', 14, 'HorizontalAlignment', 'right'  )
text( 610, 80.5, '\rightarrow Sample size \geq 600', 'FontSize', 14, 'HorizontalAlignment', 'left'  )

hold off

xticks( 300:50:950 )
xlim( [ 300, 950 ] )

grid on; box on;

set( gca, 'FontSize', 14 )
xlabel( 'Number of samples per drug' )

% Figure 3.4

numsamples = sum( ~isnan( Response_AUC) );
numperclass = cell( 21, 1 );

for k = 1:21
    
    ind =  cell2mat( DrugOrder( :, 5 ) ) == k ;
    numperclass{ k, 1 } = numsamples( 1, ind );
    
end

allperclassa = zeros( 21, 1 );
lowperclassa = zeros( 21, 1 );

for k = 1:21
    
    allperclassa( k, 1 ) = length( numperclass{ k, 1 } );
    lowperclassa( k, 1 ) = sum( numperclass{ k, 1 } <= 500 );
    
end

figure;

b1 = bar( 1:21, allperclassa, 0.9, 'FaceColor', [0 0.4470 0.7410] );

hold on

b2 = bar( 1:21, lowperclassa, 0.6, 'FaceColor', [0.8500 0.3250 0.0980] );

hold off

xtips1 = b2( 1 ).XEndPoints;
ytips1 = b2( 1 ).YEndPoints;
labels1 = string( b2( 1 ).YData );

text( xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Color', 'white', 'FontSize', 14 )

xtips1 = b1( 1 ).XEndPoints;
ytips1 = b1( 1 ).YEndPoints;
labels1 = string( b1( 1 ).YData );

text( xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'Color', 'black', 'FontSize', 14 )

xticks( 1:21 )
xticklabels( drugtypes )
xtickangle( 35)

set( gca, 'FontSize', 14 )

grid on;
box on;

ylabel( 'Number of drug compounds' )

legend( { 'Sample size \geq 600', 'Sample size \leq 500' }, 'FontSize', 14 )
legend( 'boxoff' )


%% Chapter 5

% Figure 5.1

temp = zeros( 1, 6 );

temp( 1, 1 ) = sum( nanmax( results_mean_AUCs_test, [], 2 ) < 0.5   ); %#ok<*NANMAX> 
temp( 1, 2 ) = sum( ( 0.5 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.6 )  );
temp( 1, 3 ) = sum( ( 0.6 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.7 )  );
temp( 1, 4 ) = sum( ( 0.7 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.8 )  );
temp( 1, 5 ) = sum( ( 0.8 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.9 )  );
temp( 1, 6 ) = sum( 0.9 <= nanmax( results_mean_AUCs_test, [], 2 ) );

labels={ '24', '85', '109', '43', '3' };

figure; 
pie( temp( 2:6 ), labels )

legend( { '0.5 \leq ROC-AUC_t_e_s_t < 0.6', '0.6 \leq ROC-AUC_t_e_s_t < 0.7', ...
    '0.7 \leq ROC-AUC_t_e_s_t < 0.8', '0.8 \leq ROC-AUC_t_e_s_t < 0.9', ...
    '0.9 \leq ROC-AUC_t_e_s_t' } );
legend( 'boxoff' )

set( gca, 'FontSize', 14 )

DrugOrder( nanmax( results_mean_AUCs_test, [], 2 ) < 0.5, 1:4 )
DrugOrder( nanmax( results_mean_AUCs_test, [], 2 ) > 0.9, 1:4 )

% Figure 5.2

algorithms = { 'Neural network', 'LASSO-regularised lin.regr.',	'E.net-regularised lin.regr.',	'Ridge-regularised lin.regr.',	'LASSO-regularised log.regr.', ...
    'E.net-regularised log.regr.', 'Ridge-regularised log.regr.', 'Ridge-regularised SVM', 'LASSO-regularised SVM', 'Bagged DT ensemble', 'Boosted DT ensemble', 'Random forest' };

figure;
b1 = distributionPlot( results_mean_AUCs_tr, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] );

hold on

b2 = distributionPlot( results_mean_AUCs_test, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

boxplot( results_mean_AUCs_tr, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:17.95 );
boxplot( results_mean_AUCs_test, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:18.05 );

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( 'Mean ROC-AUCs per model' )
ylim( [ 0, 1 ] )
xlim( [ 0, 19 ] )

xticklabels( [ 'Mutation-based model', 'CNV-based model', 'Hypermethylation-based model', 'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', algorithms  ] )

set( gca, 'FontSize', 14 )

legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Training', 'Testing' }, 'Location', 'southeast', 'FontSize', 14 );
legend( 'boxoff' )


[ val, ind ] = sort( results_mean_AUCs_test( :, 1 ), 'descend', 'MissingPlacement', 'last' );

DrugOrder( ind( 1:16 ), : )

sum( cell2mat( DrugOrder( :, 5 ) ) == 4 )
hygecdf( 11, 265, 19, 16, 'upper' )

val( 1:16 )
DrugOrder( ind(1:8), : )

sum( strcmp( DrugOrder( :, 3 ), 'MEK1, MEK2' ) )
hygecdf( 6, 265, 7, 8, 'upper' )

DrugOrder( ind( 1:16 ), : )
sum( strcmp( DrugOrder( :, 3 ), 'BRAF' ) )
hygecdf( 4, 265, 5, 16, 'upper' )

% Figure 5.3

high = sum( ~isnan( Response_AUC ), 1 ) >= 600;
low = sum( ~isnan( Response_AUC ), 1 ) < 500;

overf = results_mean_AUCs_tr - results_mean_AUCs_test;

figure;
b1 = distributionPlot( overf( high, : ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] );

hold on

b2 = distributionPlot( overf( low, : ), 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

b3 = boxplot( overf( high, : ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:17.95 );
b4 = boxplot( overf( low, : ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:18.05 );

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( 'Model-wise difference between training and test ROC-AUC' )
ylim( [ -0.2, 0.6 ] )
xlim( [ 0, 19 ] )

xticklabels( [ 'Mutation-based model', 'CNV-based model', 'Hypermethylation-based model', 'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', algorithms  ] )

set( gca, 'FontSize', 14 )

legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Sample size \geq 600', 'Sample size \leq 500' }, 'Location', 'southeast', 'FontSize', 14 );
legend( 'boxoff' )

% Figure C.1

indices = cell( 18, 5 );

for k = 1:18
    
    indices{ k, 1 } = find( results_mean_AUCs_tr( :, k ) < 0.6 );
    indices{ k, 2 } = find( 0.6 <= results_mean_AUCs_tr( :, k ) & results_mean_AUCs_tr( :, k ) < 0.7 );
    indices{ k, 3 } = find( 0.7 <= results_mean_AUCs_tr( :, k ) & results_mean_AUCs_tr( :, k ) < 0.8 );
    indices{ k, 4 } = find( 0.8 <= results_mean_AUCs_tr( :, k ) & results_mean_AUCs_tr( :, k ) < 0.9 );
    indices{ k, 5 } = find( 0.9 <= results_mean_AUCs_tr( :, k ) );
    
    
end

overfits = cell( 18, 5 );

for k = 1:18
    
    overfits{ k, 1 } = overf( indices{ k, 1 }, k );
    overfits{ k, 2 } = overf( indices{ k, 2 }, k );
    overfits{ k, 3 } = overf( indices{ k, 3 }, k );
    overfits{ k, 4 } = overf( indices{ k, 4 }, k );
    overfits{ k, 5 } = overf( indices{ k, 5 }, k );
    
    
end

figure;
t = tiledlayout( 5, 1 );
t.TileSpacing = 'compact';
t.Padding = 'compact';

colours = { [0 0.4470 0.7410], [0.4660, 0.6740, 0.1880], [0.9290, 0.6940, 0.1250], [0.8500 0.3250 0.0980], [0.6350 0.0780 0.1840] };
b = zeros( 1, 5 );

for k = 1:5
    
    oftemp = [];
    g = [];
    
    for l = 1:18
        
        if isempty( overfits{ l, k } )
            
            oftemp = [ oftemp; NaN ]; %#ok<*AGROW> 
            g = [ g; l ];
            
        else
        
            oftemp = [ oftemp; overfits{ l, k } ];
            g = [ g; l * ones( size( overfits{ l, k } ) ) ];
            
        end
        
    end
    
    ax = nexttile;
    
    btemp = distributionPlot( overfits( :, k ), 'color', colours{ k }, 'histOpt', 1, 'showMM', 0 )
    b( 1, k ) = btemp{ 1, 1 }( 14, 1 );
    
    hold on;
    
    boxplot( oftemp, g, 'PlotStyle','traditional', 'Widths', 0.1, 'BoxStyle', 'filled', 'Colors', colours{ k }, 'Symbol', 'k.', 'MedianStyle', 'target', 'Positions', 1:1:18 )
    
    hold off;
    
    grid on;
    box on;
    
    set( gca, 'FontSize', 14 )
    
    ylim( [ -0.15, 0.55 ] )
    xticklabels( [] )
    xlim( [ 0, 19 ] )
    
end

ylabel( t, 'Model-wise difference between training and test ROC-AUC', 'FontSize', 14 )
xtickangle( 25 )
xticklabels( [ 'Mutation-based model', 'CNV-based model', 'Hypermethylation-based model', ...
    'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', algorithms  ] )

lgd = legend( ax, b, { '0.5 \leq ROC-AUC_t_r < 0.6', '0.6 \leq ROC-AUC_t_r < 0.7', '0.7 \leq ROC-AUC_t_r < 0.8', '0.8 \leq ROC-AUC_t_r < 0.9', ...
    '0.9 \leq ROC-AUC_t_r' }, 'Orientation', 'Horizontal', 'NumColumns', 5 );
lgd.Layout.Tile = 'North';

% Figure 5.4

bestfirst = nanmax( results_mean_AUCs_test( :, 1:6 ), [], 2 );
bestsecond = nanmax( results_mean_AUCs_test( :, 7:18 ), [], 2 );

improvement = bestsecond - bestfirst;

sum( improvement > 0 )

figure;
distributionPlot( improvement, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on
    
boxplot( improvement, 'PlotStyle','traditional', 'Widths', 0.3, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
    
hold off

grid on;
box on;

ylim( [ -0.05, 0.125 ] )
yticks( -0.05:0.025:0.125 )
xticks( [] )

ylabel( 'Difference between best test ROC-AUC of 1^{st}- and 2^{nd}-step models' )
set( gca, 'FontSize', 12 )

% Figure 5.5

% Analyse the improvement for targetted and cytotoxic drugs separately
% Cytotoxic drugs: DNA replication and Cytoskeleton

index01 = cell2mat( DrugOrder( :, 5 ) ) == 17;
index02 = cell2mat( DrugOrder( :, 5 ) ) == 2;
index03 = cell2mat( DrugOrder( :, 5 ) ) == 20;

index2 = find( index01 == 1 | index02 == 1 ); % cytotoxic
index1 = find( index01 == 0 & index02 == 0 & index03 == 0 ); % targeted
index3 = find( index03 == 1 ); % other

bestfirst = nanmax( results_mean_AUCs_test( :, 1:6 ), [], 2 );
bestsecond = nanmax( results_mean_AUCs_test( :, 7:18 ), [], 2 );

improvement1 = bestsecond( index1 ) - bestfirst( index1 );
improvement2 = bestsecond( index2 ) - bestfirst( index2 );
improvement3 = bestsecond( index3 ) - bestfirst( index3 );

improvementcell = { improvement1, improvement2, improvement3 };
improvementindex = zeros( 265, 1 );
improvementindex( index1 ) = 1;
improvementindex( index2 ) = 2;
improvementindex( index3 ) = 3;

[ ~, p1 ] = ttest2( improvement1, improvement2, 'Vartype', 'unequal' ) 
[ ~, p2 ] = ttest2( improvement1, improvement3, 'Vartype', 'unequal' ) 
[ ~, p3 ] = ttest2( improvement2, improvement3, 'Vartype', 'unequal' ) 


figure;

distributionPlot( improvementcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on
    
boxplot( improvement, improvementindex, 'PlotStyle','traditional', 'Widths', 0.3, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

% p-value between groups 1 and 2

plot( [ 1, 2 ], [ 0.125, 0.125 ], 'k', 'HandleVisibility', 'off' );
plot( [ 1, 1 ], [ 0.125, 0.12 ], 'k', 'HandleVisibility', 'off' );
plot( [ 2, 2 ], [ 0.125, 0.12 ], 'k', 'HandleVisibility', 'off' );

text( 1.5, 0.13, [ 'p = ', num2str( p1 ) ], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )

% p-value between groups 2 and 3

plot( [ 2, 3 ], [ 0.1, 0.1 ], 'k', 'HandleVisibility', 'off' );
plot( [ 2, 2 ], [ 0.1, 0.095 ], 'k', 'HandleVisibility', 'off' );
plot( [ 3, 3 ], [ 0.1, 0.095 ], 'k', 'HandleVisibility', 'off' );

text( 2.5, 0.105, [ 'p = ', num2str( p3 ) ], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )

% p-value between groups 1 and 3

plot( [ 1, 3 ], [ -0.07, -0.07 ], 'k', 'HandleVisibility', 'off' );
plot( [ 1, 1 ], [ -0.07, -0.065 ], 'k', 'HandleVisibility', 'off' );
plot( [ 3, 3 ], [ -0.07, -0.065 ], 'k', 'HandleVisibility', 'off' );

text( 2, -0.085, [ 'p = ', num2str( p2 ) ], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )

hold off

grid on;
box on;

yticks( -0.075:0.025:0.125 )
ylim( [ -0.1, 0.15 ] )
xticks( 1:3 )
xticklabels( { 'Targeted drugs', 'Cytotoxic drugs', 'Other' } )
xtickangle( 0 )

ylabel( 'Difference in the best test ROC-AUC between 1^{st}- and 2^{nd}-step models' )
set( gca, 'FontSize', 14 )

% Figure 5.6

[ bestval, bestalg ] = nanmax( results_mean_AUCs_test, [], 2 );

counts = zeros( 18, 1 );

for k = 1:18
    
    counts( k, 1 ) = sum( bestalg == k );
    
end

figure; 
b = bar( counts, 'FaceAlpha', 0.7 );

ylim( [ 0, 55 ] )
xlim( [ 0, 19 ] )

xtips1 = b( 1 ).XEndPoints;
ytips1 = b( 1 ).YEndPoints;
labels1 = string( b( 1 ).YData );
text( xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12 )

xticks( 1:18 )
xtickangle( 25 )
xticklabels( [ 'Mutation-based model', 'CNV-based model', 'Hypermethylation-based model', ...
    'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', algorithms  ] )

ylabel( 'Cases where an algorithm produces the best-performing model' )

grid on;
box on;

set( gca, 'FontSize', 12 )

% Figure 5.7

bestperf = cell( 1, 18 );

for k = 1:18
    
    bestperf{ 1, k } = bestval( bestalg == k );
    
end

figure;
distributionPlot( bestperf, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( bestval, bestalg, 'PlotStyle','traditional', 'Widths', 0.1, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 3:18 )
    
hold off

box on;
grid on;

xlim( [ 0, 19] )

xticks( 1:18 )
xtickangle( 25 )
xticklabels( [ 'Mutation-based model', 'CNV-based model', 'Hypermethylation-based model', ...
    'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', algorithms  ] )

yticks( 0.45:0.05:1 )
yticklabels( { '', '0.5', '', '0.6', '', '0.7', '', '0.8', '', '0.9', '', '1' } )
ylim( [ 0.45, 1 ] )

ylabel( 'Mean test ROC-AUCs of the best-performing models' )

set( gca, 'FontSize', 14 )

% Figure 5.8

y = results_mean_AUCs_test( :, 7:18 );

g1 = cell( 265, 12 ); % algorithms
g2 = cell( 265, 12 ); % drug classes
g3 = cell( 265, 12 ); % sample numbers

drugs = {'ABL', 'DNArepl', 'EGFR', 'ERK MAPK', 'GenInt', 'IGFR', 'JNK p38', 'PI3K', 'RTK', 'TOR', 'WNT', 'ApopReg', 'CellCyc', 'ChrHistAcet', ...
    'ChrHistMeth', 'ChrOther', 'CytSkel', 'Met', 'Mit', 'Other', 'p53' };

algs = { 'NN', 'LLinR', 'ELinR', 'RLinR','LLogR', 'ELogR', 'RLogR', 'LSVM', 'RSVM', 'BagDT', 'BoostDT', 'RF' };

nums = { 'high', 'low' };

Y = sum( ~isnan( Response_AUC ) ) >= 600;

for l = 1:265
    
    for k = 1:12
        
        g1{ l, k } = algs{ k };
        
    end
    
    
    for k = 1:12
        
        g2{ l, k } = drugs{ DrugOrder{ l, 5 } };
        
    end
    
    if Y( l )
        
        for k = 1:12
        
            g3{ l, k } = nums{ 1 };
        
        end
        
    else
        
        for k = 1:12
        
            g3{ l, k } = nums{ 2 };
        
        end
        
    end
    
end

% reformat

y = reshape( y, 1, 12*265 );
g1 = reshape( g1, 1, 12*265 );
g2 = reshape( g2, 1, 12*265 );
g3 = reshape( g3, 1, 12*265 );


temp = algs;

g1num = zeros( size( g1 ) );

for j = 1:size( temp, 2 )
    
    g1num( 1, strcmp( algs{ j }, g1 ) ) = j;
    
end


temp = drugs;

g2num = zeros( size( g2 ) );

for j = 1:size( temp, 2 )
    
    g2num( 1, strcmp( temp{ j }, g2 ) ) = j;
    
end


temp = unique( g3 );

g3num = zeros( size( g3 ) );

for j = 1:size( temp, 2 )
    
    g3num( 1, strcmp( temp{ j }, g3 ) ) = j;
    
end


tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y( g2num == k )';

end

figure;

btemp = distributionPlot( tempcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( y, g2num, 'PlotStyle', 'compact', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

hold off

xlim( [ 0, 22 ] )
xticks( 1:21 )
xticklabels( drugtypes )
xtickangle( 35 )
ylabel( 'Mean test ROC-AUC' )
ylim( [ 0, 1 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure 5.9 

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

ablation_drops = cell( 265, 6 );

for drug = 1:265 % run through the drugs
    
    index = nansum( results_ablation_testAUC{ drug } ) > 0;
    temp = results_ablation_testAUC{ drug };
    
    for k = 1:6 % run through the inputs
        
        inputs = index( 1:6 );
        
        if inputs( k )
            
            removed = k;
            compare = 0; % for the original, full model
            
            for l = 7:41
                
                temp2 = combinations{ l, 1 };
                
                if ~isempty( intersect( temp2, k) )
                    
                    removed = [ removed, l ]; 
                    temp2 = setdiff( temp2, k );
                    
                    for m = 1:41
                        
                        if isequal( combinations{ m, 1 }, temp2 )
                            
                            compare = [ compare, m ]; 
                            
                        end
                        
                    end
                    
                end
                
            end
            
            tempdrop = NaN( 12, length( compare ) );
            
            tempdrop( :, 1 ) = temp( [ 1:9, 11:13 ], removed( 1 ) )./results_mean_AUCs_test( drug, 7:18 )';
            
            for p = 2:length( compare ) % run through the ablation models of interest
                
                if index( p )
                    
                    tempdrop( :, p ) = temp( [ 1:9, 11:13 ], removed( p ) )./temp( [ 1:9, 11:13 ], compare( p ) );
                end
                
            end
            
            ablation_drops{ drug, k } = tempdrop; 
            
        end
        
    end

end


datatypes = { ' somatic mutation', 'CNV', 'hypermethylation', 'tissue-descriptor', 'pathway activation', 'gene expression' };

figure;
t = tiledlayout( 6, 1 );
t.TileSpacing = 'compact';
t.Padding = 'compact';

for k = 1:6

    tempmat = [];
    s = [];
    sm = 0;
    
    for l = 1:21
        
        drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
        
        tempdrop = ablation_drops( drugs, k );
        
        temptempmat = NaN( length( drugs ), 12, 16 );
        
        for m = 1:length( drugs )
            
            if ~isempty( tempdrop{ m, 1 } )
                
                temptempmat( m, :, : ) = tempdrop{ m, 1 };
                
            end
            
        end
        
        tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ]; %#ok<*AGROW>
        s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
        
        if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
           
            sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
            
        end
        
    end
    
    nexttile;
    
    temptempmat = NaN( sm, 21 );
    
    for l = 1:21
        
        temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
        
    end
        
    b = distributionPlot( temptempmat, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

    hold on
    
    boxplot( temptempmat, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'outline', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
    
    hold off
    
    grid on;
    box on;
    
    xlim( [ 0, 22 ] )
    xticks( 1:21 )
    
    ylim( [ 0.5, 1.5 ] )
    yticks( 0.5:0.25:1.5 )
   
    if k >= 6
    
        xticklabels( drugtypes )
        xtickangle( 35 )
    
    else
        
        xticklabels( {} )
    
    end
    
    set( gca, 'FontSize', 12 )
    
    title( [  datatypes{ k }, ' data' ] )
    
end

ylabel( t, 'Ratio of mean test ROC-AUCs between models utilising particular data types and those disregarding them' )
title( t, 'Loss of predictive performance upon removing: ', 'FontWeight', 'bold', 'FontSize', 14 )

% Figures C.2-C.4

ablation_drops2 = cell( 265, 15 );

for drug = 1:265 % run through the drugs
    
    index = nansum( results_ablation_testAUC{ drug } ) > 0;
    temp = results_ablation_testAUC{ drug };
    
    for k = 1:15 % run through the pairs of inputs
        
        inputs = index( 1:6 );
        
        if sum( inputs( combi2( k, : ) ) ) == 2
            
            removed = k + 6;
            compare = 0; % for the original, full model
            
            for l = 22:41
                
                temp2 = combinations{ l, 1 };
                
                if length( intersect( temp2, combi2( k, : ) ) ) == 2
                    
                    removed = [ removed, l ]; 
                    temp2 = setdiff( temp2, combi2( k, : ) );
                    
                    for m = 1:41
                        
                        if isequal( combinations{ m, 1 }, temp2 )
                            
                            compare = [ compare, m ]; 
                            
                        end
                        
                    end
                    
                end
                
            end
            
            tempdrop = NaN( 12, length( compare ) );
            
            tempdrop( :, 1 ) = temp( [ 1:9, 11:13 ], removed( 1 ) )./results_mean_AUCs_test( drug, 7:18 )';
            
            for p = 2:length( compare ) % run through the ablation models of interest
                
                if index( p )
                    
                    tempdrop( :, p ) = temp( [ 1:9, 11:13 ], removed( p ) )./temp( [ 1:9, 11:13 ], compare( p ) );
                    
                end
                
            end
            
            ablation_drops2{ drug, k } = tempdrop; 
            
        end
        
    end

end

for n = 1:3
    
    figure;
    t = tiledlayout( 5, 1 );
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    for k = 5 * ( n - 1 ) + 1:5 * n
        
        tempmat = [];
        s = [];
        sm = 0;
        
        for l = 1:21
            
            drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
            
            tempdrop = ablation_drops2( drugs, k );
            
            temptempmat = NaN( length( drugs ), 12, 5 );
            
            for m = 1:length( drugs )
                
                if ~isempty( tempdrop{ m, 1 } )
                    
                    temptempmat( m, :, : ) = tempdrop{ m, 1 };
                    
                end
                
            end
            
            tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ];
            s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
            
            if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
                
                sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
                
            end
            
        end
        
        nexttile;
        
        temptempmat = NaN( sm, 21 );
        
        for l = 1:21
            
            temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
            
        end
        
        b = distributionPlot( temptempmat, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )
        
        hold on
        
        boxplot( temptempmat, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'outline', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
        
        hold off
        
        grid on;
        box on;
        
        xlim( [ 0, 22 ] )
        xticks( 1:21 )
        
        if ~mod( k, 5 )
            
            xticklabels( drugtypes )
            xtickangle( 35 )
            
        else
            
            xticklabels( {} )
            
        end
        
        ylim( [ 0.5, 1.5 ] )
        yticks( 0.5:0.25:1.5 )
        
        set( gca, 'FontSize', 12 )
        
        temptemp = combi2( k, : );
        title( [ datatypes{ temptemp( 1 ) }, ' and ', datatypes{ temptemp( 2 ) }, ' data' ] )
        
    end
    
    ylabel( t, 'Ratio of mean test ROC-AUCs between models utilising particular data types and those disregarding them' )
    title( t, 'Loss of predictive performance upon removing: ', 'FontWeight', 'bold', 'FontSize', 14 )
    
end

% Figures C.5-C.8

ablation_drops3 = cell( 265, 20 );

for drug = 1:265 % run through the drugs
    
    index = nansum( results_ablation_testAUC{ drug } ) > 0;
    temp = results_ablation_testAUC{ drug };
    
    for k = 1:20 % run through the combinations of 3 inputs
        
        inputs = index( 1:6 );
        
        if sum( inputs( combi3( k, : ) ) ) == 3

            ablation_drops3{ drug, k } = temp( [ 1:9, 11:13 ], k + 21 )./results_mean_AUCs_test( drug, 7:18 )'; 
            
        end
        
    end

end

for n = 1:4
    
    figure;
    t = tiledlayout( 5, 1 );
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    for k = 5 * ( n - 1 ) + 1:5 * n
        
        tempmat = [];
        s = [];
        sm = 0;
        
        for l = 1:21
            
            drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
            
            tempdrop = ablation_drops3( drugs, k );
            
            temptempmat = NaN( length( drugs ), 12 );
            
            for m = 1:length( drugs )
                
                if ~isempty( tempdrop{ m, 1 } )
                    
                    temptempmat( m, : ) = tempdrop{ m, 1 };
                    
                end
                
            end
            
            tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ];
            s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
            
                        
            if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
                
                sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
                
            end
            
        end
        
        nexttile;
        
        temptempmat = NaN( sm, 21 );
        
        for l = 1:21
            
            temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
            
        end
        
        b = distributionPlot( temptempmat, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )
        
        hold on
        
        boxplot( temptempmat, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'outline', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
        
        hold off
        
        grid on;
        box on;
        
        xlim( [ 0, 22 ] )
        xticks( 1:21 )
        
        if ~mod( k, 5 )
            
            xticklabels( drugtypes )
            xtickangle( 35 )
            
        else
            
            xticklabels( {} )
            
        end
        
        ylim( [ 0.5, 1.5 ] )
        yticks( 0.5:0.25:1.5 )
        
        set( gca, 'FontSize', 12 )
        
        temptemp = combi3( k, : );
        title( [ datatypes{ temptemp( 1 ) }, ', ', datatypes{ temptemp( 2 ) }, ' and ', datatypes{ temptemp( 3 ) }, ' data' ] )
        
    end
    
    ylabel( t, 'Ratio of mean test ROC-AUCs between models utilising particular data types and those disregarding them' )
    title( t, 'Loss of predictive performance upon removing: ', 'FontWeight', 'bold', 'FontSize', 14 )
    
end

% Figure 5.10

stds = NaN( 41, 21 );
means = NaN( 41, 21 );
medians = NaN( 41, 21 ); 
q95 = NaN( 41, 21 );
q05 = NaN( 41, 21 );
q75 = NaN( 41, 21 );
q25 = NaN( 41, 21 );

for k = 1:6
    
    tempmat = [];
    s = [];
    sm = 0;
    
    for l = 1:21
        
        drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
        
        tempdrop = ablation_drops( drugs, k );
        
        temptempmat = NaN( length( drugs ), 12, 16 );
        
        for m = 1:length( drugs )
            
            if ~isempty( tempdrop{ m, 1 } )
                
                temptempmat( m, :, : ) = tempdrop{ m, 1 };
                
            end
            
        end
        
        tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ]; %#ok<*AGROW>
        s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
        
        if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
            
            sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
            
        end
        
    end
    
    temptempmat = NaN( sm, 21 );
    
    for l = 1:21
        
        temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
        
    end
    
   stds( k, : ) = nanstd( temptempmat ); %#ok<*NANSTD> 
   means( k, : ) = nanmean( temptempmat );
   medians( k, : ) = nanmedian( temptempmat );
   q95( k, : ) = prctile( temptempmat, 95 );
   q05( k, : ) = prctile( temptempmat, 5 );
   q75( k, : ) = quantile( temptempmat, 0.75 );
   q25( k, : ) = quantile( temptempmat, 0.25 );
    
end

for k = 1:15
    
    tempmat = [];
    s = [];
    sm = 0;
    
    for l = 1:21
        
        drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
        
        tempdrop = ablation_drops2( drugs, k );
        
        temptempmat = NaN( length( drugs ), 12, 5 );
        
        for m = 1:length( drugs )
            
            if ~isempty( tempdrop{ m, 1 } )
                
                temptempmat( m, :, : ) = tempdrop{ m, 1 };
                
            end
            
        end
        
        tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ];
        s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
        
        if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
            
            sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
            
        end
        
    end
    
    temptempmat = NaN( sm, 21 );
    
    for l = 1:21
        
        temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
        
    end
    
    stds( k + 6, : ) = nanstd( temptempmat );
    means( k + 6, : ) = nanmean( temptempmat );
    medians( k + 6, : ) = nanmedian( temptempmat );
    q95( k + 6, : ) = prctile( temptempmat, 95 );
    q05( k + 6, : ) = prctile( temptempmat, 5 );
   q75( k + 6, : ) = quantile( temptempmat, 0.75 );
   q25( k + 6, : ) = quantile( temptempmat, 0.25 );
    
end

for k = 1:20
        
    tempmat = [];
    s = [];
    sm = 0;
    
    for l = 1:21
        
        drugs = find( cell2mat( DrugOrder( :, 5 ) ) == l );
        
        tempdrop = ablation_drops3( drugs, k );
        
        temptempmat = NaN( length( drugs ), 12 );
        
        for m = 1:length( drugs )
            
            if ~isempty( tempdrop{ m, 1 } )
                
                temptempmat( m, : ) = tempdrop{ m, 1 };
                
            end
            
        end
        
        tempmat = [ tempmat; temptempmat( ~isnan( temptempmat ) ) ];
        s = [ s; l * ones( size( temptempmat( ~isnan( temptempmat ) ) ) ) ];
        
        if size( temptempmat( ~isnan( temptempmat ) ), 1 ) > sm
            
            sm = size( temptempmat( ~isnan( temptempmat ) ), 1 );
            
        end
        
    end
    
    temptempmat = NaN( sm, 21 );
    
    for l = 1:21
        
        temptempmat( 1:sum( s == l ),l ) = tempmat( s == l, 1 );
        
    end
    
    stds( k + 21, : ) = nanstd( temptempmat );
    means( k + 21, : ) = nanmean( temptempmat );
    medians( k + 21, : ) = nanmedian( temptempmat ); %#ok<*NANMEDIAN> 
    q95( k + 21, : ) = prctile( temptempmat, 95 );
    q05( k + 21, : ) = prctile( temptempmat, 5 );
    q75( k + 21, : ) = quantile( temptempmat, 0.75 );
    q25( k + 21, : ) = quantile( temptempmat, 0.25 );
    
end

prexvalues = { 'Mut, ', 'CNV, ', 'Meth, ', 'Tis, ', 'PWAct, ', 'Gex, ' };
xvalues = {};
index = 1;

for k = 1:size( combinations, 1 )
    
    temp = [ prexvalues{ combinations{ k, 1 } } ];
    xvalues{ index, 1 } = temp( 1:length( temp ) - 2 ); %#ok<*SAGROW,*AGROW>
    index = index + 1;
    
end

figure;

h = heatmap( drugtypes, xvalues, 100*(1 - medians), 'Colormap',jet, 'MissingDataLabel', 'Missing', 'MissingDataColor', [ 1 1 1 ], ...
    'ColorbarVisible','on', 'CellLabelFormat', '%0.2f', 'GridVisible', 'off' );

set( gca, 'FontSize', 11 ) 
set( struct( h ).NodeChildren( 3 ), 'YTickLabelRotation', 40 );
set( struct( h ).NodeChildren( 3 ), 'XTickLabelRotation', 40 );

nanmin( nanmin( 100*(1 - medians) ) ) %#ok<*NANMIN> 
nanmax( nanmax( 100*(1 - medians) ) )

% Figure C.9

hm = NaN( 41, 21 );
hm( q75 - q25 <= 0.05 ) = ( 1 - medians( q75 - q25 <= 0.05 ) ) * 100; 

figure;

h = heatmap( drugtypes, xvalues, hm, 'Colormap',jet, 'MissingDataLabel', 'Heterogeneous', 'MissingDataColor', [ 1 1 1 ], ...
    'ColorbarVisible','on', 'CellLabelFormat', '%0.2f', 'GridVisible', 'on' );

set( gca, 'FontSize', 11 ) 
set( struct( h ).NodeChildren( 3 ), 'YTickLabelRotation', 40 );
set( struct( h ).NodeChildren( 3 ), 'XTickLabelRotation', 40 );

% Figure C.10

hm = NaN( 41, 21 );
hm( sign( q95 - 1 ) == sign( q05 - 1 ) ) = ( 1 - medians( sign( q95 - 1 ) == sign( q05 - 1 ) ) ) * 100; 

figure;

h = heatmap( drugtypes, xvalues, hm, 'Colormap',jet, 'MissingDataLabel', 'Not sign.', 'MissingDataColor', [ 1 1 1 ], ...
    'ColorbarVisible','on', 'CellLabelFormat', '%0.2f', 'GridVisible', 'on' );

set( gca, 'FontSize', 11 ) 
set( struct( h ).NodeChildren( 3 ), 'YTickLabelRotation', 40 );
set( struct( h ).NodeChildren( 3 ), 'XTickLabelRotation', 40 );

% Figure 5.11 and Table B.4

hm = NaN( 41, 21 );
hm( sign( q95 - 1 ) == sign( q05 - 1 ) & q75 - q25 <= 0.05 & ( 1 -medians ) >= 0.05 ) = ...
    ( 1 - medians( sign( q95 - 1 ) == sign( q05 - 1 ) & q75 - q25 <= 0.05 & ( 1 - medians ) >= 0.05 ) ) * 100; 

indexr = sum( ~isnan(hm), 2 ) >= 1;
indexc = sum( ~isnan(hm) ) >= 1;

hm = hm( indexr, indexc );

figure;

h = heatmap( drugtypes( indexc ), xvalues( indexr ), hm, 'Colormap',jet, 'MissingDataLabel', 'Not significant', 'MissingDataColor', [ 1 1 1 ], ...
    'ColorbarVisible','on', 'CellLabelFormat', '%0.2f', 'GridVisible', 'on' );

set( gca, 'FontSize', 12 ) 
set( struct( h ).NodeChildren( 3 ), 'YTickLabelRotation', 40 );
set( struct( h ).NodeChildren( 3 ), 'XTickLabelRotation', 40 );

% Figure 5.12 and Table B.5

load( 'PathwayOrder.mat' )
load( 'GeneOrder.mat' )

PCOrder = { 'PC_1', 'PC_2', 'PC_3', 'PC_4', 'PC_5', 'PC_6', 'PC_7' };

drugpatterns = cell( 265, 6 );

for k = 1:265 % exclude the class 'Others'
        
        str = [ 'res_', num2str( k ), '.mat' ];
        load( str, 'var_sets' )
        load( str, 'models1' )
        
        tempmutations = [];
        tempcnvs = [];
        tempmethylations = [];
        temptissues = [];
        temppathways = [];
        tempPCs = [];
        
        for l = 1:10 % run through the folds
            
            tempmutations = [ tempmutations, var_sets{ l, 5 } ];
            tempcnvs = [ tempcnvs, var_sets{ l, 6 }' ];
            tempmethylations = [ tempmethylations, var_sets{ l, 7 }' ];
            temptissues = [ temptissues, var_sets{ l, 4 } ];
            
            temppathways = [ temppathways, find( models1{ l, 5 }.Coefficients.pValue <= 0.05/11 )' ];
            tempPCs = [ tempPCs, find( models1{ l, 6 }.Coefficients.pValue <= 0.05/7 )' ];
            
        end
        
        % test if the mutations are stable for the drug in question
        % if so, add them
        
        u1 = unique( tempmutations );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( tempmutations == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 1 } = u1( logical( index ) );
    
        % test if the CNVs are stable for the drug in question
        % if so, add them
        
        u1 = unique( tempcnvs );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( tempcnvs == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 2 } = u1( logical( index ) );
        
        % test if the methylations are stable for the drug in question
        % if so, add them
        
        u1 = unique( tempmethylations );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( tempmethylations == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 3 } = u1( logical( index ) );
        
        % test if the tissues are stable for the drug in question
        % if so, add them
        
        u1 = unique( temptissues );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( temptissues == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 4 } = u1( logical( index ) );
        
        % test if the pathways are stable for the drug in question
        % if so, add them
        
        u1 = unique( temppathways );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( temppathways == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 5 } = u1( logical( index ) );
        
        % test if the PCs are stable for the drug in question
        % if so, add them
        
        u1 = unique( tempPCs );
        index = zeros( size( u1 ) );
        
        for m = 1:length( u1 ) 
            
            index( m ) = sum( tempPCs == u1( m ) ) >= 7;
            
        end
        
        drugpatterns{ k, 6 } = u1( logical( index ) );
        
        k        
        
end

figure;
t = tiledlayout( 5, 4 ); 
t.TileSpacing = 'compact';

indices = [ 1:19, 21 ];
features = cell( 20, 1 );

for j = 1:20
    
    drugpatterns_temp = drugpatterns( cell2mat( DrugOrder( :, 5 ) ) == indices( j ), : );
    n = size( drugpatterns_temp, 1 );
    
    % count how often a feature occurs (in %)

    freqs_temp = cell( 2, 6 );
    
    for k = 1:6
        
        list = [];
        
        for l = 1:n
            
            list = [ list, drugpatterns_temp{ l, k } ];
            
        end
        
        ulist = unique( list );
        flist = zeros( size( ulist ) );
        
        for l = 1:length( ulist )
            
            flist( l ) = sum( list == ulist( l ) );
            
        end
        
        freqs_temp{ 1, k } = ulist;
        freqs_temp{ 2, k } = 100 * flist./n;
    
    end
    
    freqs_temp{ 3, 1 } = GeneOrder( freqs_temp{ 1, 1 } );
    freqs_temp{ 3, 2 } = GeneOrder( freqs_temp{ 1, 2 } );
    freqs_temp{ 3, 3 } = GeneOrder( freqs_temp{ 1, 3 } );
    
    freqs_temp{ 3, 4 } = TissueOrder( freqs_temp{ 1, 4 }, 2 );
    
    freqs_temp{ 3, 5 } = PathwayOrder( freqs_temp{ 1, 5 } );
    
    freqs_temp{ 3, 6 } = PCOrder( freqs_temp{ 1, 6 } );
    
    % sort out features that are just found once in order to not overcrowd the
    % plot
    
    for k = 1:6
        
        temp = find( freqs_temp{ 2, k } ~= 100 * 1/n );
        temp2 = freqs_temp{ 2, k };
        temp3 = freqs_temp{ 3, k }
        
        freqs_temp{ 4, k } = temp2( temp2 ~= 100 * 1/n );
        freqs_temp{ 5, k } = temp3( temp ); %#ok<FNDSB> 
        
    end
    
    % reordering
    
    for k = 1:6
        
        if size( freqs_temp{ 2, k }, 1 ) ~= 1
            
            freqs_temp{ 2, k } = freqs_temp{ 2, k }';
            
        end
        
    end

    features{ j, 1 } = freqs_temp;

    nexttile;
    stepplot2( freqs_temp, j );
    
    xticks( 0:20:100 )
    
    title( drugtypes{ indices( j ), 1 } )
    
end

xlabel( t, 'Percentage of drugs per class' )
ylabel( t, 'Number of significant predictive features' )

% Figure 5.13

drugpatterns_temp = drugpatterns( cell2mat( DrugOrder( :, 5 ) ) == indices( 19 ), : );
n = size( drugpatterns_temp, 1 );
    
% count how often a feature occurs (in %)
    
freqs_temp = cell( 2, 6 );
    
for k = 1:6

    list = [];

    for l = 1:n

        list = [ list, drugpatterns_temp{ l, k } ];

    end

    ulist = unique( list );
    flist = zeros( size( ulist ) );

    for l = 1:length( ulist )

        flist( l ) = sum( list == ulist( l ) );

    end

    freqs_temp{ 1, k } = ulist;
    freqs_temp{ 2, k } = 100 * flist./n;

end
    
freqs_temp{ 3, 1 } = GeneOrder( freqs_temp{ 1, 1 } );
freqs_temp{ 3, 2 } = GeneOrder( freqs_temp{ 1, 2 } );
freqs_temp{ 3, 3 } = GeneOrder( freqs_temp{ 1, 3 } );

freqs_temp{ 3, 4 } = TissueOrder( freqs_temp{ 1, 4 }, 2 );

freqs_temp{ 3, 5 } = PathwayOrder( freqs_temp{ 1, 5 } );

freqs_temp{ 3, 6 } = PCOrder( freqs_temp{ 1, 6 } );

% sort out features that are just found once in order to not overcrowd the
% plot

for k = 1:6

    temp = find( freqs_temp{ 2, k } ~= 100 * 1/n );
    temp2 = freqs_temp{ 2, k };
    temp3 = freqs_temp{ 3, k };

    freqs_temp{ 4, k } = temp2( temp2 ~= 100 * 1/n );
    freqs_temp{ 5, k } = temp3( temp ); %#ok<FNDSB> 

end
    
% reordering
    
for k = 1:6

    if size( freqs_temp{ 2, k }, 1 ) ~= 1

        freqs_temp{ 2, k } = freqs_temp{ 2, k }';

    end

end

freqs_19 = freqs_temp;
drugpatterns_19 = drugpatterns_temp;
 
[ ~, t ] = featureplot( freqs_19 )  

% Single features of interest

% H&L tissue and CNV in CDKN2A 

ind = find( cell2mat( DrugOrder( :, 5 ) ) == 19 );

tis = TissueType( :, 10 ) == 1;
cnv = abs( CNV( :, 8605 ) ) == 1;

meant = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meant( k, 1 ) = nanmean( Response_AUC( tis, ind( k ) ) );
    meant( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, tis ), ind( k ) ) );
    meant( k, 3 ) = sign( meant( k, 1 ) - meant( k, 2 ) );
    meant( k, 4 ) = nanmedian( Response_AUC( tis, ind( k ) ) );
    meant( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, tis ), ind( k ) ) );
    meant( k, 6 ) = sign( meant( k, 4 ) - meant( k, 5 ) );
    
end

meanc = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanc( k, 1 ) = nanmean( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 3 ) = sign( meanc( k, 1 ) - meanc( k, 2 ) );
    meanc( k, 4 ) = nanmedian( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 6 ) = sign( meanc( k, 4 ) - meanc( k, 5 ) );
    
end

% Hypermethylations

m3 =  strcmp( GeneOrder, 'ADHFE1' ) == 1 ;
m4 =  strcmp( GeneOrder, 'PABPC3' ) == 1 ;
m5 =  strcmp( GeneOrder, 'RRS1' ) == 1 ; % 3 and 5 are correlated

m3 = abs( Hypermethylation( :, m3 ) ) == 1; %#ok<*NASGU> % -1 in seven of them, 1,3,4 not
m4 = abs( Hypermethylation( :, m4 ) ) == 1; % -1 everywhere
m5 = abs( Hypermethylation( :, m5 ) ) == 1;

meanm = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanm( k, 1 ) = nanmean( Response_AUC( m4, ind( k ) ) );
    meanm( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, m4 ), ind( k ) ) );
    meanm( k, 3 ) = sign( meanm( k, 1 ) - meanm( k, 2 ) );
    meanm( k, 4 ) = nanmedian( Response_AUC( m4, ind( k ) ) );
    meanm( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, m4 ), ind( k ) ) );
    meanm( k, 6 ) = sign( meanm( k, 4 ) - meanm( k, 5 ) );
    
end

% Figure 5.14

drugpatterns_14 = drugpatterns( cell2mat( DrugOrder( :, 5 ) ) == 14, : );

% count how often a feature occurs (in %)

freqs_14 = cell( 2, 6 );

for k = 1:6
    
    list = [];
    
    for l = 1:10
        
        list = [ list, drugpatterns_14{ l, k } ];
        
    end
    
    ulist = unique( list );
    flist = zeros( size( ulist ) );
    
    for l = 1:length( ulist )
        
        flist( l ) = sum( list == ulist( l ) );
        
    end
    
    freqs_14{ 1, k } = ulist;
    freqs_14{ 2, k } = 100 * flist./10;
    
end

freqs_14{ 3, 1 } = GeneOrder( freqs_14{ 1, 1 } );
freqs_14{ 3, 2 } = GeneOrder( freqs_14{ 1, 2 } );
freqs_14{ 3, 3 } = GeneOrder( freqs_14{ 1, 3 } );

freqs_14{ 3, 4 } = TissueOrder( freqs_14{ 1, 4 }, 2 );

freqs_14{ 3, 5 } = PathwayOrder( freqs_14{ 1, 5 } );

freqs_14{ 3, 6 } = PCOrder( freqs_14{ 1, 6 } );

% sort out features that are just found once in order to not overcrowd the
% plot

for k = 1:6
  
    temp = find( freqs_14{ 2, k } ~= 100 * 1/10 );  
    temp2 = freqs_14{ 2, k };
    temp3 = freqs_14{ 3, k }
    
    freqs_14{ 4, k } = temp2( temp2 ~= 100 * 1/10 );
    freqs_14{ 5, k } = temp3( temp ); %#ok<FNDSB> 
    
end

featureplot( freqs_14 )

% H&L tissue and CNV in CDKN2A 

ind = find( cell2mat( DrugOrder( :, 5 ) ) == 14 );

tis = TissueType( :, 10 ) == 1;
cnv = abs( CNV( :, 8605 ) ) == 1;

meant = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meant( k, 1 ) = nanmean( Response_AUC( tis, ind( k ) ) );
    meant( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, tis ), ind( k ) ) );
    meant( k, 3 ) = sign( meant( k, 1 ) - meant( k, 2 ) );
    meant( k, 4 ) = nanmedian( Response_AUC( tis, ind( k ) ) );
    meant( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, tis ), ind( k ) ) );
    meant( k, 6 ) = sign( meant( k, 4 ) - meant( k, 5 ) );
    
end

meanc = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanc( k, 1 ) = nanmean( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 3 ) = sign( meanc( k, 1 ) - meanc( k, 2 ) );
    meanc( k, 4 ) = nanmedian( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 6 ) = sign( meanc( k, 4 ) - meanc( k, 5 ) );
    
end

% Hypermethylations

m1 =  strcmp( GeneOrder, 'TOX2' ) == 1 ;
m2 =  strcmp( GeneOrder, 'PIK3R1' ) == 1 ;
m3 =  strcmp( GeneOrder, 'ADHFE1' ) == 1 ;
m4 =  strcmp( GeneOrder, 'PABPC3' ) == 1 ;
m5 =  strcmp( GeneOrder, 'RRS1' ) == 1 ;
m6 =  strcmp( GeneOrder, 'EXD3' ) == 1 ;
m7 =  strcmp( GeneOrder, 'RBP1' ) == 1 ;
m8 =  strcmp( GeneOrder, 'MT1E' ) == 1 ;
m9 =  strcmp( GeneOrder, 'EEF1D' ) == 1 ;
m10 =  strcmp( GeneOrder, 'NAPRT1' ) == 1 ;

% Check if combinations are perfectly correlated and therefore redundant

m1 = abs( Hypermethylation( :, m1 ) ) == 1;
m2 = abs( Hypermethylation( :, m2 ) ) == 1;
m3 = abs( Hypermethylation( :, m3 ) ) == 1;
m4 = abs( Hypermethylation( :, m4 ) ) == 1;
m5 = abs( Hypermethylation( :, m5 ) ) == 1;
m6 = abs( Hypermethylation( :, m6 ) ) == 1;
m7 = abs( Hypermethylation( :, m7 ) ) == 1;
m8 = abs( Hypermethylation( :, m8 ) ) == 1;
m9 = abs( Hypermethylation( :, m9 ) ) == 1;
m10 = abs( Hypermethylation( :, m10 ) ) == 1;

ms = [ m1, m2, m3, m4, m5, m6, m7, m8, m9, m10 ];


eq = zeros( 10, 10 );

for k = 1:10
    
    for j = 1:10
        
        eq( j, k ) = isequal( ms( :, j ), ms( :, k ) );
        
    end
    
end

% m3 and m5 are perfectly correlated
% m9 and m10 are perfectly correlated

meanm = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanm( k, 1 ) = nanmean( Response_AUC( m9, ind( k ) ) );
    meanm( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, m9 ), ind( k ) ) );
    meanm( k, 3 ) = sign( meanm( k, 1 ) - meanm( k, 2 ) );
    meanm( k, 4 ) = nanmedian( Response_AUC( m9, ind( k ) ) );
    meanm( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, m9 ), ind( k ) ) );
    meanm( k, 6 ) = sign( meanm( k, 4 ) - meanm( k, 5 ) );
    
end

% -1,-1,-1,-1,-1,-1,-1,-1 

% Somatic mutation in MYC

mut = SomaticMutation( :, 7040 ) == 1;

meant = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meant( k, 1 ) = nanmean( Response_AUC( mut, ind( k ) ) );
    meant( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, mut ), ind( k ) ) );
    meant( k, 3 ) = sign( meant( k, 1 ) - meant( k, 2 ) );
    meant( k, 4 ) = nanmedian( Response_AUC( mut, ind( k ) ) );
    meant( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, mut ), ind( k ) ) );
    meant( k, 6 ) = sign( meant( k, 4 ) - meant( k, 5 ) );
    
end

% Figure 5.15


drugpatterns_3 = drugpatterns( cell2mat( DrugOrder( :, 5 ) ) == 3, : );

% count how often a feature occurs (in %)

freqs_3 = cell( 2, 6 );

for k = 1:6
    
    list = [];
    
    for l = 1:9
        
        list = [ list, drugpatterns_3{ l, k } ];
        
    end
    
    ulist = unique( list );
    flist = zeros( size( ulist ) );
    
    for l = 1:length( ulist )
        
        flist( l ) = sum( list == ulist( l ) );
        
    end
    
    freqs_3{ 1, k } = ulist;
    freqs_3{ 2, k } = 100 * flist./9;
    
end

freqs_3{ 3, 1 } = GeneOrder( freqs_3{ 1, 1 } );
freqs_3{ 3, 2 } = GeneOrder( freqs_3{ 1, 2 } );
freqs_3{ 3, 3 } = GeneOrder( freqs_3{ 1, 3 } );

freqs_3{ 3, 4 } = TissueOrder( freqs_3{ 1, 4 }, 2 );

freqs_3{ 3, 5 } = PathwayOrder( freqs_3{ 1, 5 } );

PCOrder = { 'PC_1', 'PC_2', 'PC_3', 'PC_4', 'PC_5', 'PC_6', 'PC_7' }

freqs_3{ 3, 6 } = PCOrder( freqs_3{ 1, 6 } );

% sort out features that are just found once in order to not overcrowd the
% plot

for k = 1:6
  
    temp = find( freqs_3{ 2, k } ~= 100 * 1/9 );  
    temp2 = freqs_3{ 2, k };
    temp3 = freqs_3{ 3, k }
    
    freqs_3{ 4, k } = temp2( temp2 ~= 100 * 1/9 );
    freqs_3{ 5, k } = temp3( temp ); %#ok<FNDSB> 
    
end


h = featureplot( freqs_3 );

% Single features of interest

ind = find( cell2mat( DrugOrder( :, 5 ) ) == 3 );

tis1 = TissueType( :, 10 ) == 1; % H&L tissue
tis2 = TissueType( :, 15 ) == 1; % Oesophagus
cnv = abs( CNV( :, 7782 ) ) == 1; % CN amplifications in ErbB2
meth = Hypermethylation( :, 9250 ) == 1; % Hypermethylation in THY1


meant = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meant( k, 1 ) = nanmean( Response_AUC( tis1, ind( k ) ) );
    meant( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, tis1 ), ind( k ) ) );
    meant( k, 3 ) = sign( meant( k, 1 ) - meant( k, 2 ) );
    meant( k, 4 ) = nanmedian( Response_AUC( tis1, ind( k ) ) );
    meant( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, tis1 ), ind( k ) ) );
    meant( k, 6 ) = sign( meant( k, 4 ) - meant( k, 5 ) );
    
end % +1 mostly

meant = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meant( k, 1 ) = nanmean( Response_AUC( tis2, ind( k ) ) );
    meant( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, tis2 ), ind( k ) ) );
    meant( k, 3 ) = sign( meant( k, 1 ) - meant( k, 2 ) );
    meant( k, 4 ) = nanmedian( Response_AUC( tis2, ind( k ) ) );
    meant( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, tis2 ), ind( k ) ) );
    meant( k, 6 ) = sign( meant( k, 4 ) - meant( k, 5 ) );
    
end % -1 mostly

meanc = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanc( k, 1 ) = nanmean( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 3 ) = sign( meanc( k, 1 ) - meanc( k, 2 ) );
    meanc( k, 4 ) = nanmedian( Response_AUC( cnv, ind( k ) ) );
    meanc( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, cnv ), ind( k ) ) );
    meanc( k, 6 ) = sign( meanc( k, 4 ) - meanc( k, 5 ) );
    
end % -1 six times

meanm = zeros( length( ind ), 6 );

for k = 1:length( ind )
    
    meanm( k, 1 ) = nanmean( Response_AUC( meth, ind( k ) ) );
    meanm( k, 2 ) = nanmean( Response_AUC( setdiff( 1:265, meth ), ind( k ) ) );
    meanm( k, 3 ) = sign( meanm( k, 1 ) - meanm( k, 2 ) );
    meanm( k, 4 ) = nanmedian( Response_AUC( meth, ind( k ) ) );
    meanm( k, 5 ) = nanmedian( Response_AUC( setdiff( 1:265, meth ), ind( k ) ) );
    meanm( k, 6 ) = sign( meanm( k, 4 ) - meanm( k, 5 ) );
    
end % -1

% Case study 1: Lapatinib

load( 'res_40.mat' )

% Figure 5.16

figure;

btemp = distributionPlot( Response_AUC( :, 40 ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on
    
boxplot( Response_AUC( :, 40 ), 'PlotStyle','traditional', 'Widths', 0.3, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
    
hold off

grid on;
box on;

ylim( [ 0, 1.1 ] )
yticks( 0.1:0.1:1 )
xticks( [] )

set( gca, 'FontSize', 14 )
ylabel( 'Cellular responsiveness to Lapatinib (AUC)', 'FontSize', 16 )

% Figure 5.17

index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;

temp = cell2mat( AUCs( :, 2 ) );
temp = temp( :, 1:6 );
temp = reshape( temp, 60, 1 );
xvals = [ ones( 10, 1 ); 2* ones( 10, 1 ); 3 * ones( 10, 1 ); 4 * ones( 10, 1 ); 5 * ones( 10, 1 ); 6 * ones( 10, 1 ); ]

figure;

b = boxchart( xvals, temp, 'BoxFaceColor', [0 0.4470 0.7410], 'MarkerColor', [0 0.4470 0.7410] )

ylim( [ 0, 1 ] )

xlabels = { '', 'Somatic mutation', 'CNV', 'Hypermethylation', 'Tissue descriptors', 'Pathway activation', 'Gene expression', '' };

xticklabels( xlabels )
xtickangle( 35 )

title( ' ' )

ylabel( 'Distribution of test ROC-AUCs' )

set( gca, 'FontSize', 14 )

box on;
grid on;

% Figure 5.18

[ ~, ~, meancoeffs_n40 ] = heatmapplot( 40, coeffs, DrugOrder );

set( 0, 'ShowHiddenHandles', 'on' )
allhnds = get( 0, 'Children' );
ht = findall( allhnds, 'Tag', 'HeatMapAxes' );
set( ht, 'FontSize', 10 )

nanmean( meancoeffs_n40') %#ok<*UDIM> 
nanmedian( meancoeffs_n40')

% Figure 5.19

h2 = clustergramplot2( 40, teststats, results_ablation_testAUC, results_mean_AUCs_test, DrugOrder );

% Figure 5.20

[ completesets40, importance40 ] = varplot4( 40, Response_AUC( :, 40 ), data );
[ ~, ~, completesetsext_40 ]  = varplot6vert( completesets40, 40, data );

% Figure 5.21a

find( strcmp( GeneOrder, 'SPOCK3' ) == 1 ) % 15455
find( strcmp( GeneOrder, 'MICALL1' ) == 1 ) % 2125
find( strcmp( GeneOrder, 'CORO6' ) == 1 ) % 11493
find( strcmp( GeneOrder, 'FAM171A1' ) == 1 ) % 8666
find( strcmp( GeneOrder, 'ACADM' ) == 1 ) % 4481

temp1 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); % only SPOCK3, 4
temp2 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 1 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); % only MICALL1, 2
temp3 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); % only CORO6, 3
temp4 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 1 & SomaticMutation( :, 4481 ) == 0, 40 ); % only FAM171A1, 8
temp5 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 1, 40 ); % only ACADM, 3

% Combinations of two

temp6 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 1 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); %  1
temp7 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); %  1
temp8 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 0 ...
        & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 1, 40 ); % 1
temp9 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 1 ...
        & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); % 1
temp10 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
         & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 1, 40 ); % 1
temp11 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
         & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 1 & SomaticMutation( :, 4481 ) == 1, 40 ); % 1

% Combinations of three

temp12 = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 1 ...
         & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 1 & SomaticMutation( :, 4481 ) == 1, 40 ); % 1
temp13 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 0 ...
         & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 1 & SomaticMutation( :, 4481 ) == 0, 40 ); % 1
temp14 = Response_AUC( SomaticMutation( :, 15455 ) == 1 & SomaticMutation( :, 2125 ) == 0 ...
         & SomaticMutation( :, 11493 ) == 1 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 1, 40 ); % 1


% No combinations of four or five are found in the data

tempnone = Response_AUC( SomaticMutation( :, 15455 ) == 0 & SomaticMutation( :, 2125 ) == 0 ...
           & SomaticMutation( :, 11493 ) == 0 & SomaticMutation( :, 8666 ) == 0 & SomaticMutation( :, 4481 ) == 0, 40 ); % 358

tempdata = [ temp11; temp3; temp9; temp2; tempnone; temp8; temp1; temp10; temp4; temp14; temp5; temp7; temp12; temp6; temp13 ];
g = [ ones( size( temp11 ) ); 2 * ones( size( temp3 ) ); 3 * ones( size( temp9 ) ); 4 * ones( size( temp2 ) ); 5 * ones( size( tempnone ) ); ...
      6 * ones( size( temp8 ) ); 7 * ones( size( temp1 ) ); 8 * ones( size( temp10 ) ); 9 * ones( size( temp4 ) ); 10 * ones( size( temp14 ) ); ...
      11 * ones( size( temp5 ) ); 12 * ones( size( temp7 ) ); 13 * ones( size( temp12 ) ); 14 * ones( size( temp6 ) ); ...
      15 * ones( size( temp13 ) ) ];

figure;

b = distributionPlot( { temp11; temp3; temp9; temp2; tempnone; temp8; temp1; temp10; temp4; temp14; temp5; temp7; temp12; temp6; temp13 }, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( tempdata, g, 'PlotStyle','traditional', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

grid on;  box on;

xticks( 0:16 )
xticklabels( { '', 'FAM171A1 & ACADM', 'CORO6', 'MICALL1 & CORO6', 'MICALL1', 'None', 'SPOCK3 & ACADM', ...
             'SPOCK3', 'CORO6 & ACADM', 'FAM171A1', 'SPOCK3 & CORO6 & ACADM', 'ACADM', 'SPOCK3 & CORO6', ...
             'MICALL1 & FAM171A1 & ACADM', 'MICALL1 & SPOCK3', 'SPOCK3 & CORO6 & FAM171A1', '' } );
xtickangle( 25 )
xlim( [ 0, 16 ] )

ylabel( 'Cellular sensitivity to Lapatinib (AUC)' )
ylim( [ 0, 1 ] )

set( gca, 'FontSize', 14 )

% Figure 5.21b

temp1 = Response_AUC( abs( CNV( :, 7782 ) ) == 1 & abs( CNV( :, 15369 ) ) == 0 & abs( CNV( :, 10859 ) ) == 0, 40 ); % only ERBB2
temp2 = Response_AUC( abs( CNV( :, 15369 ) ) == 1 & abs( CNV( :, 7782 ) ) == 0 & abs( CNV( :, 10859 ) ) == 0, 40 ); % only ZNF
temp3 = Response_AUC( abs( CNV( :, 10859 ) ) == 1 & abs( CNV( :, 7782 ) ) == 0 & abs( CNV( :, 15369 ) ) == 0, 40 ); % only RAD21

% There is no partial set with 2/3 combinations that has a measured response to Lapatinib

tempset = Response_AUC( abs( CNV( :, 15369 ) ) == 1 & abs( CNV( :, 7782 ) ) == 1 & abs( CNV( :, 10859 ) ) == 1, 40 ); % all
tempnone = Response_AUC( abs( CNV( :, 15369 ) ) == 0 & abs( CNV( :, 7782 ) ) == 0 & abs( CNV( :, 10859 ) ) == 0, 40 );

tempdata = [ tempnone; temp2; temp3; temp1; tempset ];
g = [ ones( size( tempnone ) ); 2 * ones( size( temp2 ) ); 3 * ones( size( temp3 ) ); 4 * ones( size( temp1 ) ); 5 * ones( size( tempset ) ) ];

figure;

distributionPlot( { tempnone; temp2; temp3; temp1; tempset }, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( tempdata, g, 'PlotStyle','traditional', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

grid on;  box on;

xticklabels( {'None', 'RAD21', 'ERBB2', 'ZNF292/SYNCRIP', 'All three' } )
xlabel( '' )

ylabel( 'Cellular sensitivity to Lapatinib (AUC)' )
ylim( [ 0, 1.1 ] )

set( gca, 'FontSize', 14 )

% Figure 5.21c

temp1 = Response_AUC( TissueType( :, 6 ) == 1, 40 ); % Breast, 13
temp2 = Response_AUC( TissueType( :, 10 ) == 1, 40 ); % Hematopoietic and lymphoid cells, 141
temp3 = Response_AUC( TissueType( :, 11 ) == 1, 40 ); % Kidney, 9
temp4 = Response_AUC( TissueType( :, 15 ) == 1, 40 ); % Oesophagus, 9
temp5 = Response_AUC( TissueType( :, 29 ) == 1, 40 ); % Urinary tract, 2
others = Response_AUC( TissueType( :, 6 ) == 0 & TissueType( :, 10 ) == 0 & ...
    TissueType( :, 11 ) == 0 & TissueType( :, 15 ) == 0 & TissueType( :, 29 ) == 0, 40 ); % All others, 213

tempdata = [ temp2; others; temp1; temp3; temp4; temp5 ];
g = [ ones( size( temp2 ) ); 2 * ones( size( others ) ); 3 * ones( size( temp1 ) ); ...
      4 * ones( size( temp3 ) ); 5 * ones( size( temp4 ) ); 6 * ones( size( temp5 ) ) ];

figure;

b = distributionPlot( { temp2; others; temp1; temp3; temp4; temp5 }, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( tempdata, g, 'PlotStyle','traditional', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

grid on;  box on;

xticklabels( { 'H&L tissue', 'Others', 'Breast', 'Kidney', 'Oesophagus', 'Urinary tract' } )
xlabel( '' )

ylabel( 'Cellular sensitivity to Lapatinib (AUC)' )
ylim( [ 0, 1.1 ] )

set( gca, 'FontSize', 14 )

% Continuous-valued features

% Pathways

means_pw = zeros( 10, 11 );
sign_pw = zeros( 10, 11 );

for k = 1:10
   
    temp = models1{ k, 5 }.Coefficients.pValue;
    sign_pw( k, temp <= 0.05/11  ) = 1;
    
    means_pw( k, : ) = temp;
    
end

means_pw = nanmedian( means_pw );
PathwayOrder( sum( sign_pw ) > 5 )
means_pw( sum( sign_pw ) > 5 )

% NFkB, 6/10 times

% PCs

means_pc = zeros( 10, 7 );
sign_pc = zeros( 10, 7 );

for k = 1:10
   
    temp = models1{ k, 6 }.Coefficients.pValue;
    sign_pc( k, temp <= 0.05/7  ) = 1;
    
    means_pc( k, : ) = temp;
    
end

means_pc = nanmean( means_pc );
means_pc( sum( sign_pc ) > 5 )

% Case study 2: Refametinib

load( 'res_263.mat' )

nanmin( Response_AUC( :, 263 ) )
nanmax( Response_AUC( :, 263 ) )

% Figure 5.22

figure;

btemp = distributionPlot( Response_AUC( :, 263 ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on
    
boxplot( Response_AUC( :, 263 ), 'PlotStyle','traditional', 'Widths', 0.3, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
    
hold off

grid on;
box on;

ylim( [ 0, 1.1 ] )
yticks( 0.1:0.1:1 )
xticks( [] )

set( gca, 'FontSize', 14 )
ylabel( 'Cellular responsiveness to Refametinib (AUC)', 'FontSize', 16 )

% Figure 5.23

index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;

temp = cell2mat( AUCs( :, 2 ) );
temp = temp( :, 1:6 );
temp = reshape( temp, 60, 1 );
xvals = [ ones( 10, 1 ); 2* ones( 10, 1 ); 3 * ones( 10, 1 ); 4 * ones( 10, 1 ); 5 * ones( 10, 1 ); 6 * ones( 10, 1 ); ]

figure;

b = boxchart( xvals, temp, 'BoxFaceColor', [0 0.4470 0.7410], 'MarkerColor', [0 0.4470 0.7410] )

ylim( [ 0, 1 ] )

xlabels = { '', 'Somatic mutation', 'CNV', 'Hypermethylation', 'Tissue descriptors', 'Pathway activation', 'Gene expression', '' };

xticklabels( xlabels )
xtickangle( 35 )

title( ' ' )

ylabel( 'Distribution of test ROC-AUCs' )

set( gca, 'FontSize', 14 )

box on;
grid on;

results_mean_AUCs_test( 263, : )
nanmax( results_mean_AUCs_test( 263, 7:18 ) ) - nanmax( results_mean_AUCs_test( 263, 1:6 ) )

% Figure 5.24

[ ~, ~, meancoeffs_n263 ] = heatmapplot( 263, coeffs, DrugOrder );

set( 0, 'ShowHiddenHandles', 'on' )
allhnds = get( 0, 'Children' );
ht = findall( allhnds, 'Tag', 'HeatMapAxes' );
set( ht, 'FontSize', 10 )

nanmean( meancoeffs_n263')
nanmedian( meancoeffs_n263')

% Figure 5.25

h2 = clustergramplot2( 263, teststats, results_ablation_testAUC, results_mean_AUCs_test, DrugOrder );

set( 0, 'ShowHiddenHandles', 'on' )
allhnds = get( 0, 'Children' );
ht = findall( allhnds, 'Tag', 'HeatMapAxes' );
set( ht, 'FontSize', 10 )

% Figure 5.26

[ completesets263, importance263 ] = varplot4( 263, Response_AUC( :, 263 ), data );
completesetsext_263  = varplot6vert( completesets263, 263, data );

% Figure 5.27a

temp1 = Response_AUC( TissueType( :, 5 ) == 1, 263 ); % Bone, 38
temp2 = Response_AUC( TissueType( :, 6 ) == 1, 263 ); % Breast, 43
temp3 = Response_AUC( TissueType( :, 10 ) == 1, 263 ); % Hematopoietic & lymphoid cells, 148
temp4 = Response_AUC( TissueType( :, 12 ) == 1, 263 ); % Large intestine, 46
temp5 = Response_AUC( TissueType( :, 14 ) == 1, 263 ); % Lung, 154
temp6 = Response_AUC( TissueType( :, 22 ) == 1, 263 ); % Skin, 50
temp7 = Response_AUC( TissueType( :, 27 ) == 1, 263 ); % Thyroid, 13

others = Response_AUC( TissueType( :, 5 ) == 0 & TissueType( :, 6 ) == 0 & ...
         TissueType( :, 10 ) == 0 & TissueType( :, 12 ) == 0 & TissueType( :, 14 ) == 0 & ...
         TissueType( :, 22 ) == 0 & TissueType( :, 27 ) == 0, 263 ); % all others, 338

tempdata = [ temp3; temp1; temp2; temp5; others; temp4; temp7; temp6 ];
g = [ ones( size( temp3 ) ); 2 * ones( size( temp1 ) ); 3 * ones( size( temp2 ) ); ...
      4 * ones( size( temp5 ) ); 5 * ones( size( others ) ); ...
      6 * ones( size( temp4 ) ); 7 * ones( size( temp7 ) ); 8 * ones( size( temp6 ) ) ];

figure;

b = distributionPlot( { temp3; temp1; temp2; temp5; others; temp4; temp7; temp6 }, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( tempdata, g, 'PlotStyle','traditional', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

grid on;  box on;

xticklabels( { 'H&L tissue', 'Bone', 'Breast', 'Lung', 'Others', 'Large intestine', 'Thyroid', 'Skin' } )
xlabel( '' )
xtickangle( 25 )
xlim( [ 0, 9 ] )

ylabel( 'Cellular sensitivity to Refametinib (AUC)' )
ylim( [ 0, 1.05 ] )

set( gca, 'FontSize', 14 )

% Figure 5.27b

temp1 = Response_AUC( SomaticMutation( :, 9638 ) == 1 & SomaticMutation( :, 6460 ) == 0, 263 ); % only BRAF, 96
temp2 = Response_AUC( SomaticMutation( :, 9638 ) == 0 & SomaticMutation( :, 6460 ) == 1, 263 ); % only KRAS, 116

tempset = Response_AUC( SomaticMutation( :, 9638 ) == 1 & SomaticMutation( :, 6460 ) == 1, 263 ); % both, 6
tempnone = Response_AUC( SomaticMutation( :, 9638 ) == 0 & SomaticMutation( :, 6460 ) == 0, 263 );

tempdata = [ tempnone; temp2; temp1; tempset ];
g = [ ones( size( tempnone ) ); 2 * ones( size( temp2 ) ); 3 * ones( size( temp1 ) ); 4 * ones( size( tempset ) ) ];

figure;

b = distributionPlot( { tempnone; temp2; temp1; tempset }, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( tempdata, g, 'PlotStyle','traditional', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

grid on;  box on;

xticks( 0:5 )
xticklabels( { '', 'None', 'KRAS', 'BRAF', 'BRAF & KRAS' } )
xlim( [ 0, 5 ] )

ylabel( 'Cellular sensitivity to Refametinib (AUC)' )
ylim( [ 0, 1.1 ] )

set( gca, 'FontSize', 14 )

% Figure 5.28a

means_pw = zeros( 10, 11 );
sign_pw = zeros( 10, 11 );

for k = 1:10
   
    temp = models1{ k, 5 }.Coefficients.pValue;
    sign_pw( k, temp <= 0.05/11  ) = 1;
    
    means_pw( k, : ) = temp;
    
end

means_pw = nanmean( means_pw );
PathwayOrder( sum( sign_pw ) > 5 )
means_pw( sum( sign_pw ) > 5 )

% {'PI3K'}, p = 1.0074e-04, found to be significant 10/10 times    

figure;
plot( PathwayActivation( :, 6 ), Response_AUC( :, 263 ), 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410] )
grid on

xlabel( 'Activation score for the PI3K pathway' )
ylabel( 'Sensitivity to Refametinib (AUC)' )
set( gca, 'FontSize', 14 )

index = ~isnan( PathwayActivation( :, 6 ).*Response_AUC( :, 263 ) );
corrcoef( PathwayActivation( index, 6 ), Response_AUC( index, 263 ) ) 

% Figure 5.28b

[ ~, score,~ ] = pca( GeneExpression );
scores = score( :, 1:2 );

means_pc = zeros( 10, 7 );
sign_pc = zeros( 10, 7 );

for k = 1:10
   
    temp = models1{ k, 6 }.Coefficients.pValue;
    sign_pc( k, temp <= 0.05/7  ) = 1;
    
    means_pc( k, : ) = temp;
    
end

means_pc = nanmean( means_pc );
means_pc( sum( sign_pc ) > 5 )

% PC 1, found to be significant 10/10 times
% PC 2, found to be significant 10/10 times

figure; 
scatter( scores( :, 1 ), scores( :, 2 ), 28, Response_AUC( :, 263 ), 'filled' )
grid on

xlabel( 'PC 1' )
ylabel( 'PC 2' )
set( gca, 'FontSize', 14 )

colormap parula

box on;
grid on;

% Figure 5.29a

figure; 
scatter( scores( :, 1 ), scores( :, 2 ), 28, PathwayActivation( :, 6 ), 'o', 'filled' )

grid on; box on;

xlabel( 'PC 1' )
ylabel( 'PC 2' )
set( gca, 'FontSize', 14 )

% Figure 5.29b

temp5 = TissueType( :, 14 ) == 1; % Lung, 154
index1 = strcmp( Celllinesinfo( :, 8 ), 'lung_NSCLC' ); %109
index2 = strcmp( Celllinesinfo( :, 8 ), 'lung_SCLC' ); % 61

index1 = index1( 2:969 );
index2 = index2( 2:969 );

nsclc = temp5 .* index1;
sclc = temp5 .* index2;

temp1 = TissueType( :, 5 ) == 1; % Bone, 38
temp2 = TissueType( :, 6 ) == 1; % Breast, 43
temp3 = TissueType( :, 10 ) == 1; % Hematopoietic, 148
temp4 = TissueType( :, 12 ) == 1; % Large intestine, 46
temp6 = TissueType( :, 22 ) == 1; % Skin, 50
temp7 = TissueType( :, 27 ) == 1; % Thyroid, 13

nanmedian( Response_AUC( temp1, 263 ) ) % 3
nanmedian( Response_AUC( temp2, 263 ) ) % 4
nanmedian( Response_AUC( temp3, 263 ) ) % 2
nanmedian( Response_AUC( temp4, 263 ) ) % 6
nanmedian( Response_AUC( temp6, 263 ) ) % 8
nanmedian( Response_AUC( temp7, 263 ) ) % 7
nanmedian( Response_AUC( logical( nsclc ), 263 ) ) % 5
nanmedian( Response_AUC( logical( sclc ), 263 ) ) % 1

% Ordered wrt to median AUC: SCLC, hematopoietic, bone, breast, NSCLC,
% large intestine, thyroid, skin
 
colours = parula( 8 );
 
figure; 
scatter( scores( temp3, 1 ), scores( temp3, 2 ), 32, colours( 1, : ), '<', 'filled' )

hold on

scatter( scores( logical( sclc ), 1 ), scores( logical( sclc ), 2 ), 32, colours( 2, : ), 'o', 'filled' )
scatter( scores( temp1, 1 ), scores( temp1, 2 ), 32, colours( 3, : ), 'd', 'filled' )
scatter( scores( temp2, 1 ), scores( temp2, 2 ), 50, colours( 4, : ), 'p', 'filled' )
scatter( scores( logical( nsclc ), 1 ), scores( logical( nsclc ), 2 ), 32, colours( 5, : ), 'v', 'filled' )
scatter( scores( temp4, 1 ), scores( temp4, 2 ), 50, colours( 6, : ), 'h', 'filled' )
scatter( scores( temp7, 1 ), scores( temp7, 2 ), 32, colours( 7, : ), 's', 'filled' )
scatter( scores( temp6, 1 ), scores( temp6, 2 ), 32, colours( 8, : ), '^',  'filled' )

grid on; box on;

xlabel( 'PC 1' )
ylabel( 'PC 2' )
set( gca, 'FontSize', 14 )

lgd = legend( { 'H&L tissue','SCLC' 'Bone', 'Breast', 'NSCLC', 'L. intestine', 'Thyroid', 'Skin' }, 'Location', 'North' )
lgd.NumColumns = 2;


%% Chapter 6

% Run the one-step models on every drug in the GDSC data set

for drug = 1:2

    [ models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part ] = onestepmodel( drug );

    parsave_one( strcat( "res_", num2str( drug ) ), models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part );
    
end

% Sort the results into arrays: Use AUCs as measures of predictive performance

mean_AUCs_tr = zeros( 265, 12 );
mean_AUCs_test = zeros( 265, 12 );

for count = 1:265

    str = strcat( 'res_', num2str( count ), '.mat' )
    
    load( str, 'AUCs' );
    
    mean_AUCs_tr( count, : ) = nanmean( cell2mat( AUCs( :, 1 ) ) );
    mean_AUCs_test( count, : ) = nanmean( cell2mat( AUCs( :, 2 ) ) );
    

end

% Figure 6.1 

temp = zeros( 1, 6 );

temp( 1, 1 ) = sum( nanmax( results_mean_AUCs_test, [], 2 ) < 0.5   ); 
temp( 1, 2 ) = sum( ( 0.5 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.6 )  );
temp( 1, 3 ) = sum( ( 0.6 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.7 )  );
temp( 1, 4 ) = sum( ( 0.7 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.8 )  );
temp( 1, 5 ) = sum( ( 0.8 <= nanmax( results_mean_AUCs_test, [], 2 ) ) & ( nanmax( results_mean_AUCs_test, [], 2 ) < 0.9 )  );
temp( 1, 6 ) = sum( 0.9 <= nanmax( results_mean_AUCs_test, [], 2 ) );

labels = { '24', '85', '109', '43', '3' };

figure; 
pie( temp( 2:6 ), labels )

legend( { '0.5 \leq ROC-AUC_t_e_s_t < 0.6', '0.6 \leq ROC-AUC_t_e_s_t < 0.7', ...
          '0.7 \leq ROC-AUC_t_e_s_t < 0.8', '0.8 \leq ROC-AUC_t_e_s_t < 0.9', ...
          '0.9 \leq ROC-AUC_t_e_s_t' } );
legend( 'boxoff' )

set( gca, 'FontSize', 14 )
labels = { '24', '85', '109', '43', '3' };

% Figure 6.2

temp = zeros( 1, 6 );

temp( 1, 1 ) = sum( nanmax( mean_AUCs_test, [], 2 ) < 0.5   ); 
temp( 1, 2 ) = sum( ( 0.5 <= nanmax( mean_AUCs_test, [], 2 ) ) & ( nanmax( mean_AUCs_test, [], 2 ) < 0.6 )  );
temp( 1, 3 ) = sum( ( 0.6 <= nanmax( mean_AUCs_test, [], 2 ) ) & ( nanmax( mean_AUCs_test, [], 2 ) < 0.7 )  );
temp( 1, 4 ) = sum( ( 0.7 <= nanmax( mean_AUCs_test, [], 2 ) ) & ( nanmax( mean_AUCs_test, [], 2 ) < 0.8 )  );
temp( 1, 5 ) = sum( ( 0.8 <= nanmax( mean_AUCs_test, [], 2 ) ) & ( nanmax( mean_AUCs_test, [], 2 ) < 0.9 )  );
temp( 1, 6 ) = sum( 0.9 <= nanmax( mean_AUCs_test, [], 2 ) );

labels = { '6', '71', '122', '55', '11' };

figure; 
p = pie( temp( 2:6 ), labels )

lgd = legend( { '0.5 \leq ROC-AUC_t_e_s_t < 0.6', '0.6 \leq ROC-AUC_t_e_s_t < 0.7', ...
                '0.7 \leq ROC-AUC_t_e_s_t < 0.8', '0.8 \leq ROC-AUC_t_e_s_t < 0.9', ...
                '0.9 \leq ROC-AUC_t_e_s_t' } ); 
legend( 'boxoff' )

set( gca, 'FontSize', 14 )
labels = { '6', '71', '122', '55', '11' };

% Enrichment

temp = nanmax( mean_AUCs_test' ); 
[ val, ind ] = nanmin( temp ) % % 0.5425, 193
DrugOrder( ind ) % VX-702

[ val, ind ] = nanmax( temp ) % 0.9340, 263
DrugOrder( ind ) % RDEA119_2

sum( temp >= 0.9 ) % 11


% Models with test ROC-AUCs over 0.9

find( nanmax( mean_AUCs_test, [], 2 ) >= 0.9 )
DrugOrder( ans, : ) %#ok<NOANS> 

hygecdf( 5, 265, 19, 11, 'upper' ) * 5 % ERK MAPK: 2.2119e-05
hygecdf( 3, 265, 9, 11, 'upper' ) * 5  % EGFR: 1.8564e-04
hygecdf( 1, 265, 9, 11, 'upper' ) * 5 % ABL: 0.0482, not significant in multiple testing
hygecdf( 1, 265, 3, 11, 'upper' ) * 5  % p53: 0.0046

% Models with test ROC-AUCs under 0.6

find( nanmax( mean_AUCs_test, [], 2 ) <= 0.6 )
DrugOrder( ans, : ) %#ok<NOANS> 

hygecdf( 1, 265, 9, 6, 'upper' ) * 6 % ABL: 0.0144
hygecdf( 1, 265, 6, 6, 'upper' ) * 6  % JNK: 0.0062
hygecdf( 1, 265, 7, 6, 'upper' ) * 6  % TOR: 0.0086
hygecdf( 1, 265, 11, 6, 'upper' ) * 6  % cell cycle: 0.0215
hygecdf( 1, 265, 74, 6, 'upper' ) * 6  % other: 0.5356


numsamples = sum( ~isnan( Response_AUC) );

hygecdf( length( intersect( find( nanmax( mean_AUCs_test, [], 2 ) <= 0.6 ), find( numsamples <= 500 ) ) ), 265, 48, 6, 'upper' ) % not significant
hygecdf( length( intersect( find( nanmax( mean_AUCs_test, [], 2 ) <= 0.6 ), find( numsamples > 600 ) ) ), 265, 217, 6, 'upper' ) % not significant

hygecdf( length( intersect( find( nanmax( mean_AUCs_test, [], 2 ) > 0.9 ), find( numsamples <= 500 ) ) ), 265, 48, 11, 'upper' ) % not significant
hygecdf( length( intersect( find( nanmax( mean_AUCs_test, [], 2 ) > 0.9 ), find( numsamples > 600 ) ) ), 265, 217, 11, 'upper' ) % not significant

bestres_one = nanmax( mean_AUCs_test, [], 2 );
bestres_two = nanmax( results_mean_AUCs_test( :, 7:18 ), [], 2 );

% Figure 6.3

improvement_two_one = bestres_two' - bestres_one';

figure;
btemp = distributionPlot( improvement_two_one', 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on
    
boxplot( improvement_two_one, 'PlotStyle','traditional', 'Widths', 0.3, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )
    
hold off

grid on;
box on;

ylim( [ -0.25, 0.05 ] )
yticks( -0.25:0.05:0.05 )
yticklabels( -0.25:0.05:0.05 )
xticks( [] )

set( gca, 'FontSize', 14 )
ylabel( '\Delta_{perf} between two-step and one-step models', 'FontSize', 16 )

% Additional information

% 20 negative outliers
% maximum: 0.036738
% minimum: -0.19507
% 25th percentile: -0.040309
% 75th percentile: -0.0081791
% lower adjacent: -0.088298

sum( improvement_two_one >= 0 ) % 37
nanmedian( improvement_two_one ) % -0.0225
iqr( improvement_two_one ) % 0.0321

% Enrichment for sample size

numsamples = sum( ~isnan( Response_AUC) );

temp = find( bestres_two - bestres_one <= -0.0882984 );
temp2 = find( improvement_two_one >= 0 ); 

hygecdf( length( intersect( temp, find( numsamples <= 500 ) ) ), 265, 48, length( temp ), 'upper' ) % 4.3868e-15
hygecdf( length( intersect( temp, find( numsamples > 600 ) ) ), 265, 217, length( temp ), 'upper' ) % 1.0000, no enrichment

% Drugs for which one-step models drastically outperform two-step models are enriched for small sample sizes

hygecdf( length( intersect( temp2, find( numsamples <= 500 ) ) ), 265, 48, length( temp2 ), 'upper' ) % 0.8451, no enrichment
hygecdf( length( intersect( temp2, find( numsamples > 600 ) ) ), 265, 217, length( temp2 ), 'upper' ) % 0.0629, no enrichment

% Enrichment for drug classes

d = unique( cell2mat( DrugOrder( temp, 5 ) ) );

for k = 1:length( d )
    
    dtemp = d( k );
    hygecdf( sum( cell2mat( DrugOrder( temp, 5 ) ) == d( k ) ), 265, sum( cell2mat( DrugOrder( :, 5 ) ) == d( k ) ), length( temp ), 'upper' ) 
    hygecdf( sum( cell2mat( DrugOrder( temp, 5 ) ) == d( k ) ), 265, sum( cell2mat( DrugOrder( :, 5 ) ) == d( k ) ), length( temp ), 'upper' ) * length( d )
    drugtypes( d( k ) )
    
end

% Only ABL signalling is significantly enriched here with a p-value of p = 6.2132e-06

% Figure 6.4

figure;
b1 = distributionPlot( mean_AUCs_tr, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] )

hold on

b2 = distributionPlot( mean_AUCs_test, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] )

b3 = boxplot( mean_AUCs_tr, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:11.95 )
b4 = boxplot( mean_AUCs_test, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:12.05 )

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( 'Mean ROC-AUCs per model' )
ylim( [ 0, 1 ] )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 13 )

lgd = legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Training', 'Testing' }, 'Location', 'southeast' );
legend( 'boxoff' )

% Figure C.11

ov1 = mean_AUCs_tr - mean_AUCs_test;
ov2 = results_mean_AUCs_tr - results_mean_AUCs_test;

indices = cell( 12, 5 );

for k = 1:12
    
    indices{ k, 1 } = find( mean_AUCs_tr( :, k ) < 0.6 );
    indices{ k, 2 } = find( 0.6 <= mean_AUCs_tr( :, k ) & mean_AUCs_tr( :, k ) < 0.7 );
    indices{ k, 3 } = find( 0.7 <= mean_AUCs_tr( :, k ) & mean_AUCs_tr( :, k ) < 0.8 );
    indices{ k, 4 } = find( 0.8 <= mean_AUCs_tr( :, k ) & mean_AUCs_tr( :, k ) < 0.9 );
    indices{ k, 5 } = find( 0.9 <= mean_AUCs_tr( :, k ) );
    
    
end

overfits1 = cell( 12, 5 );

for k = 1:12
    
    overfits1{ k, 1 } = ov1( indices{ k, 1 }, k );
    overfits1{ k, 2 } = ov1( indices{ k, 2 }, k );
    overfits1{ k, 3 } = ov1( indices{ k, 3 }, k );
    overfits1{ k, 4 } = ov1( indices{ k, 4 }, k );
    overfits1{ k, 5 } = ov1( indices{ k, 5 }, k );
    
    
end

figure;
t = tiledlayout( 5, 1 );
t.TileSpacing = 'compact';
t.Padding = 'compact';

colours = { [0 0.4470 0.7410], [0.4660, 0.6740, 0.1880], [0.9290, 0.6940, 0.1250], [0.8500 0.3250 0.0980], [0.6350 0.0780 0.1840] };
b = zeros( 1, 5 );

for k = 1:5
    
    oftemp = [];
    g = [];
    
    for l = 1:12
        
        if isempty( overfits1{ l, k } )
            
            oftemp = [ oftemp; NaN ];
            g = [ g; l ];
            
        else
        
            oftemp = [ oftemp; overfits1{ l, k } ];
            g = [ g; l * ones( size( overfits1{ l, k } ) ) ];
            
        end
        
    end
    
    ax = nexttile;
    
    btemp = distributionPlot( overfits1( :, k ), 'color', colours{ k }, 'histOpt', 1, 'showMM', 0 )
    b( 1, k ) = btemp{ 1, 1 }( 7, 1 );
    
    hold on;
    
    boxplot( oftemp, g, 'PlotStyle','traditional', 'Widths', 0.1, 'BoxStyle', 'filled', 'Colors', colours{ k }, 'Symbol', 'k.', 'MedianStyle', 'target', 'Positions', 1:1:12 )
    
    hold off;
    
    grid on;
    box on;
    
    set( gca, 'FontSize', 14 )
    
    ylim( [ -0.2, 0.6 ] )
    yticks( -0.2:0.2:0.6 )
    xticklabels( [] )
    xlim( [ 0, 13 ] )
    
end

yticks( -0.2:0.2:0.6 )
ylabel( t, 'Model-wise difference between training and test ROC-AUC', 'FontSize', 14 )
xtickangle( 25 )
xticklabels( algorithms )

lgd = legend( ax, b, { '0.5 \leq ROC-AUC_t_r < 0.6', '0.6 \leq ROC-AUC_t_r < 0.7', '0.7 \leq ROC-AUC_t_r < 0.8', '0.8 \leq ROC-AUC_t_r < 0.9', ...
              '0.9 \leq ROC-AUC_t_r' }, 'Orientation', 'Horizontal', 'NumColumns', 5 );
lgd.Layout.Tile = 'North';

% Figure 6.5

improvement_test = results_mean_AUCs_test( :, 7:18 ) - mean_AUCs_test;
improvement_train = results_mean_AUCs_tr( :, 7:18 ) - mean_AUCs_tr;

figure;
b1 = distributionPlot( improvement_train, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] )

hold on

b2 = distributionPlot( improvement_test, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] )

b3 = boxplot( improvement_train, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:11.95 )
b4 = boxplot( improvement_test, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:12.05 )

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( '\Delta_{perf} per model' )
ylim( [ -0.3, 0.3 ] )
yticks( -0.3:0.1:0.3 )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 13 )

lgd = legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Training', 'Testing' }, 'Location', 'southeast' );
legend( 'boxoff' )

% Figure C.12

figure;

b1 = distributionPlot( ov2( high, 7:18 ) - ov1( high, : ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] )

hold on

b2 = distributionPlot( ov2( low, 7:18 ) - ov1( low, : ), 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] )

b3 = boxplot( ov2( high, 7:18 ) - ov1( high, : ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:11.95 )
b4 = boxplot( ov2( low, 7:18 ) - ov1( low, : ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:12.05 )

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( '\Delta_{over} between two-step and one-step models' )
ylim( [ -0.25, 0.25 ] )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 13 )

lgd = legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Sample size \geq 600', 'Sample size \leq 500' }, 'Location', 'southeast' );
legend( 'boxoff' )

% Figure 6.6a

[ bestval, bestalg ] = nanmax( results_mean_AUCs_test, [], 2 );
[ bestval1, bestalg1 ] = nanmax( mean_AUCs_test, [], 2 );

counts = zeros( 12, 1 );
counts2 = zeros( 12, 1 );

for k = 1:12
    
    counts( k, 1 ) = sum( bestalg1 == k );
    counts2( k, 1 ) = sum( bestalg == k + 6 );
    
end

figure; 
b = bar( [counts2, counts ], 1, 'FaceAlpha', 0.8 );

ylim( [ 0, 100 ] )

xtips1 = b( 1 ).XEndPoints;
ytips1 = b( 1 ).YEndPoints;
labels1 = string( b( 1 ).YData );
text( xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )

xtips2 = b( 2 ).XEndPoints;
ytips2 = b( 2 ).YEndPoints;
labels2 = string( b( 2 ).YData );
text( xtips2, ytips2, labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )

xlim( [ 0, 13 ] )
xticks( 1:12 )
xtickangle( 25 )
xticklabels( algorithms )

ylabel( 'Cases of algorithms producing the best-performing models' )

grid on;
box on;

set( gca, 'FontSize', 14 )

lgd = legend( { 'Two-step models', 'One-step models' }, 'Location', 'northwest', 'FontSize', 14 );
legend( 'boxoff' )

% Figure 6.6b

bestperf1 = cell( 1, 12 );
bestperf = cell( 1, 12 );

for k = 1:12
    
    bestperf1{ 1, k } = bestval1( bestalg1 == k );
    bestperf{ 1, k } = bestval( bestalg == k + 6 );
    
end

figure;
btemp = distributionPlot( bestperf, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] )

hold on
    
btemp2 = distributionPlot( bestperf1, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] )
b3 = boxplot( bestval( bestalg >= 7 ), bestalg( bestalg >= 7 ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:11.95 )
b4 = boxplot( bestval1, bestalg1, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:12.05 )    
hold off

box on;
grid on;

xlim( [ 0, 13] )

xticks( 1:12 )
xtickangle( 25 )
xticklabels( algorithms )

yticks( 0.45:0.05:1 )
yticklabels( { '', '0.5', '', '0.6', '', '0.7', '', '0.8', '', '0.9', '', '1' } )
ylim( [ 0.45, 1 ] )

ylabel( 'Mean test ROC-AUCs of the best-performing models' )

set( gca, 'FontSize', 14 )

lgd = legend( [ btemp{ 1, 1 }( 1, 1 ), btemp2{ 1, 1 }( 1, 1 ) ], { 'Two-step models', 'One-step models' }, 'Location', 'southwest', 'FontSize', 14 );
legend( 'boxoff' )

% Figure 6.7 - see Figure 5.8

% Figure C.13

figure;

t = tiledlayout( 4, 3 );

t.TileSpacing = 'compact';

for k = 7:18
   
    nexttile;
    histogram( results_mean_AUCs_test( :, k ) )
    xlabel( algorithms{ k - 6 } )
    
    ylim( [ 0, 60 ] )
    yticks( 0:20:60 )
    
    xlim( [ 0.4, 1 ] )
    xticks( 0.4:0.1:1 )
    
    grid on;
    box on;
    
end

title( t, 'Distribution of predictive performance across distinct algorithms' ) 

% Table B.6

for k = 7:18
    
    nanvar( results_mean_AUCs_test( :, k ) ) %#ok<*NANVAR> 
    
end

nanvar( results_mean_AUCs_test( :, 11 ) )/nanvar( results_mean_AUCs_test( :, 16 ) ) 

range( nanmedian( results_mean_AUCs_test( :, 7:18 )  ) ) 


temp1 = results_mean_AUCs_test( Y == 0, 7:18 );
temp1 = reshape( temp1, 1, size( temp1, 1 ) * 12 );

temp2 = results_mean_AUCs_test( Y == 1, 7:18 );
temp2 = reshape( temp2, 1, size( temp2, 1 ) * 12 );

nanvar( temp1 )
nanvar( temp2 )

% Figure C.14

figure;

t3 = tiledlayout( 1, 2 );

t3.TileSpacing = 'compact';

nexttile;

temp = results_mean_AUCs_test( Y == 1, 7:18 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \geq 600' )

yticks( 0:25:300 )
yticklabels( { '0', '', '50', '', '100', '', '150', '', '200', '', '250', '', '300' } )

xlim( [ 0.3, 1 ] )
xticks( 0.3:0.1:1 )

grid on;
box on;

nexttile;

temp = results_mean_AUCs_test( Y == 0, 7:18 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \leq 500' )

xlim( [ 0.3, 1 ] )
xticks( 0.3:0.1:1 )

grid on;
box on;

title( t3, 'Distribution of predictive performances between drugs with large and small sample sizes' ) 

% Table B.7

for k = 1:21
    
    k
    temp = results_mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == k, 7:18 );
    nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )
    
end

temp1 = results_mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 18, 7:18 ); % smallest
temp1 = nanvar( reshape( temp1, 1, size( temp1, 1 ) * 12 ) );

temp2 = results_mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 4, 7:18 ); % largest
temp2 = nanvar( reshape( temp2, 1, size( temp2, 1 ) * 12 ) );

temp2/temp1

% Figure C.15

figure;

t2 = tiledlayout( 7, 3 );

t2.TileSpacing = 'compact';

for k = 1:21
   
    nexttile;
    
    temp = results_mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == k, 7:18 );
    histogram( reshape( temp, 1, size( temp, 1 ) * 12 ), 10 )
    xlabel( drugtypes{ k } )
    
    temp = yticks;
    
    if temp( end ) == 50
       
        yticks( 0:25:50 )
        
    end
    
    xlim( [ 0.4, 1 ] )
    xticks( 0.4:0.1:1 )
    
    grid on;
    box on;
    
end

title( t2, 'Distribution of predictive performance across distinct drug classes' ) 

% Figure C.16

figure;

t4 = tiledlayout( 4, 3 ); 

t4.TileSpacing = 'compact';

for k = 7:18
   
    ax = nexttile;
    
    h1 = histogram( results_mean_AUCs_test( Y==1, k ) - nanmean( results_mean_AUCs_test( Y==1, k ) ) ); % large sample size
    h1.Normalization = 'probability';
    h1.BinWidth = 0.05;
    
    hold on;
    
    h2 = histogram( results_mean_AUCs_test( Y==0, k ) - nanmean( results_mean_AUCs_test( Y==0, k ) ), 'FaceColor', [0.6350, 0.0780, 0.1840] ); % small sample size
    h2.Normalization = 'probability';
    h2.BinWidth = 0.05;
    
    hold off
    
    xlabel( algorithms{ k - 6 } )
    xlim( [ -0.4, 0.4] )
    xticks( -0.4:0.1:0.4 )
    
    grid on;
    box on;
    
    if k == 18
        
        lgd = legend( ax, [ h1, h2 ], { 'Sample size \geq 600', 'Sample size \leq 500' }, 'Orientation', 'Horizontal', ...
            'NumColumns', 2, 'FontSize', 11 );
        lgd.Layout.Tile = 'North';
        
    end
    
end

title( t4, 'Distribution of predictive performances across distinct algorithms for large and small sample sizes' ) 

% Table B.8

for k = 7:18
   
    k
    nanvar( results_mean_AUCs_test( Y==0, k ) )
    nanvar( results_mean_AUCs_test( Y==1, k ) )
    
end

nanvar( results_mean_AUCs_test( Y==0, 7 ) )/nanvar( results_mean_AUCs_test( Y==0, 16 ) )

% Table 6.1

[ ~, tbl1, stats1, ~ ] = anovan( y, { g1 g3 }, 'model', 'interaction', 'varnames', { 'Alg','Num' } );

% Table 6.2

[ ~, tbl2, stats2, ~ ] = anovan( y, { g1 g3 }, 'model', 'linear', 'varnames', { 'Alg','Num' } );

% Figure 6.8

figure;
[ c, m, h, gnames ] = multcompare( stats2, 'Dimension', 1 )

barcoords = [ 0.660898514936552,0.68544502963;0.63850202370969,0.663048538407713; ...
              0.645594813374003,0.670141328072026; 0.664850462570252,0.689396977268275; ...
              0.66502897853191,0.689575493229933; 0.664406110583397,0.68895262528142; ...
              0.664219233212693,0.688765747910716; 0.663867725020073,0.688414239718096; ...
              0.664414115280844,0.6889606299788671; 0.664364101791083,0.688910616489106; ...
              0.664364829026275,0.688911343724298; 0.663957662951625,0.688504177649648 ];

figure;
[ c2, m2, h2, gnames2 ] = multcompare( stats2, 'Dimension', 2 )

barcoords2 = [ 0.6220816029229,0.629884231529042; 0.715407045000381,0.723209673606523 ];

figure;

t5 = tiledlayout( 2, 1 ); 

t5.TileSpacing = 'compact';

ax = nexttile;

hold on;

for k = 1:12

    plot( barcoords( k, : ), [ k, k ], 'Color', [ 0, 0.4470, 0.7410 ], 'LineWidth', 2 );
    p1 = plot( m( k, 1 ), 13-k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [ 0, 0.4470, 0.7410 ], 'MarkerFaceColor', [ 0, 0.4470, 0.7410 ] )

end

hold off

ylim( [ 0, 13 ] )
yticks( 1:12 )
yticklabels( algorithms( 12:-1:1 ) )
ytickangle( 25 )

xlim( [ 0.6 0.75 ] )
xticks( 0.6:0.025:0.75)
xticklabels( { '0.60', '', '0.65', '', '0.70', '', '0.75' } )

title( 'Algorithms' )

grid on
box on

set( gca, 'FontSize', 14 )


ax = nexttile;

hold on

for k = 1:2

    plot( barcoords2( k, : ), [ k, k ], 'Color', [ 0, 0.4470, 0.7410 ], 'LineWidth', 2 );
    p1 = plot( m2( k, 1 ), k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [ 0, 0.4470, 0.7410 ], 'MarkerFaceColor', [ 0, 0.4470, 0.7410 ] )

end

hold off

ylim( [ 0, 3 ] )
yticks( 1:2 )
yticklabels( { 'Sample size \leq 500', 'Sample size \geq 600' } )
ytickangle( 25 )

xlim( [ 0.60 0.75 ] )
xticks( 0.60:0.025:0.75)
xticklabels( { '0.60', '', '0.65', '', '0.70', '', '0.75' } )

title( 'Sample size' )

grid on
box on

set( gca, 'FontSize', 14 )

ylabel( t5, '\bf{Factor levels}', 'FontSize', 16 )
xlabel( t5, '\bf{Mean test ROC-AUC}', 'FontSize', 16 )

% Figure 6.9

figure;
b = distributionPlot( results_mean_AUCs_test( :, 7:18 ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( results_mean_AUCs_test( :, 7:18 ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

xtickangle( 35 )

grid on;
box on;

ylabel( 'Mean test ROC-AUCs' )
ylim( [ 0.3, 1 ] )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 14 )

% Figure 6.10

tempcell = cell( 1, 2 );

tempcell{ :, 1 } = reshape( results_mean_AUCs_test( Y == 1, 7:18 ), 1, 217 * 12 );
tempcell{ :, 2 } = reshape( results_mean_AUCs_test( Y == 0, 7:18 ), 1, 48 * 12 );

temp = [ tempcell{ :, 1 }, tempcell{ :, 2 } ];
tempind = [ ones( 1, 2604 ), 2 * ones( 1, 576 ) ];

figure;

btemp = distributionPlot( tempcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( temp, tempind,  'PlotStyle', 'compact', 'Widths', 0.15, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

xticks( 1:2 )
xticklabels( { 'Sample size \geq 600', 'Sample size \leq 500' } )
ylabel( 'Mean test ROC-AUC' )
ylim( [ 0.3, 1 ] )
xlim( [ 0, 3 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure 6.11

figure;
[ c, m, h, gnames ] = multcompare( stats1, 'Dimension', [1 2] );

barcoords2 = [ 0.612267279634369,0.640750822859278; 0.589870788407506,0.618354331632415; 0.59696357807182,0.625447121296728; ...
               0.616219227268069,0.644702770492978; 0.616397743229727,0.644881286454635; 0.615774875281213,0.644258418506122; ...
               0.61558799791051,0.644071541135418; 0.615236489717889,0.643720032942798; 0.615782879978661,0.644266423203569; ...
               0.615732866488899,0.644216409713808; 0.615733593724091,0.644217136949; 0.615326427649441,0.64380997087435; ...
               0.70559272171185,0.734076264936759; 0.683196230484987,0.711679773709896; 0.690289020149301,0.71877256337421; ...
               0.70954466934555,0.738028212570459; 0.709723185307208,0.738206728532117; 0.709100317358694,0.737583860583603; ...
               0.708913439987991,0.7373969832129; 0.708561931795371,0.73704547502028; 0.709108322056142,0.737591865281051; ...
               0.709058308566381,0.737541851791289; 0.709059035801572,0.737542579026481; 0.708651869726922,0.737135412951831 ];

m( 2, 1 ) - m( 11, 1 )
m( 21, 1 ) - m( 23, 1 )

nanmean( results_mean_AUCs_test( Y == 0, 7:18 ) )
max( nanmean( results_mean_AUCs_test( Y == 0, 7:18 ) ) ) - min( nanmean( results_mean_AUCs_test( Y == 0, 7:18 ) ) )

nanmean( results_mean_AUCs_test( Y == 1, 7:18 ) )
max( nanmean( results_mean_AUCs_test( Y == 1, 7:18 ) ) ) - min( nanmean( results_mean_AUCs_test( Y == 1, 7:18 ) ) )

figure;

hold on

for k = 1:12

    plot( barcoords2( k, : ), [ k, k ], 'Color', [ 0, 0.4470, 0.7410 ], 'LineWidth', 2 );
    p1 = plot( m( 13 - k, 1 ), k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [ 0, 0.4470, 0.7410 ], 'MarkerFaceColor', [ 0, 0.4470, 0.7410 ] )

    plot( barcoords2( k + 12, : ), [ k, k ], 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2 );
    p2 = plot( m( 25 - k, 1 ), k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.6350 0.0780 0.1840], 'MarkerFaceColor', [0.6350 0.0780 0.1840] )

end

hold off

ylim( [ 0, 13 ] )
yticks( 1:12 )
yticklabels( algorithms( 12:-1:1 ) )
ytickangle( 25 )

xlim( [ 0.55 0.8 ] )
xticks( 0.55:0.025:0.8)
xticklabels( { '0.55', '', '0.6', '', '0.65', '', '0.7', '', '0.75', '', '0.8' } )

ylabel( '\bf{Algorithms}' )

grid on
box on

xlabel( 'Mean test ROC-AUC' )

set( gca, 'FontSize', 14 )

lg = legend( [ p1, p2 ], { 'Sample size \leq 500', 'Sample size \geq 600' }, 'Orientation', 'horizontal', 'FontSize', 14, 'Location', 'NorthOutside' )
legend( 'boxoff' )
title( lg, 'Sample size' )

% Figure C.17

figure;
b1 = distributionPlot( results_mean_AUCs_test( Y==1, 7:18 ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] );

hold on

b2 = distributionPlot( results_mean_AUCs_test( Y==0, 7:18 ), 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

boxplot( results_mean_AUCs_test( Y==1, 7:18 ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:11.95 );
boxplot( results_mean_AUCs_test( Y==0, 7:18 ), 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:12.05 );

hold off

xtickangle( 25 )

grid on;
box on;

ylabel( 'Mean test ROC-AUCs per model' )
ylim( [ 0, 1 ] )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 14 )

legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Sample size \geq 600', 'Sample size \leq 500' }, 'Location', 'southeast', 'FontSize', 14 );
legend( 'boxoff' )

% Table B.9

for k = 1:12

    k
    
    nanmedian( results_mean_AUCs_test( Y==1, k + 6 )  )
    nanmean( results_mean_AUCs_test( Y==1, k + 6 )  )
    nanmedian( results_mean_AUCs_test( Y==0, k + 6 )  )
    nanmean( results_mean_AUCs_test( Y==0, k + 6 )  )
    nanmean( results_mean_AUCs_test( Y==1, k + 6 )  ) - nanmean( results_mean_AUCs_test( Y==0, k + 6 )  )
    nanmedian( results_mean_AUCs_test( Y==1, k + 6 )  ) - nanmedian( results_mean_AUCs_test( Y==0, k + 6 )  )

end

% Table 6.3

drugs = {'ABL', 'DNArepl', 'EGFR', 'ERK MAPK', 'GenInt', 'IGFR', 'JNK p38', 'PI3K', 'RTK', 'TOR', 'WNT', 'ApopReg', 'CellCyc', 'ChrHistAcet', ...
         'ChrHistMeth', 'ChrOther', 'CytSkel', 'Met', 'Mit', 'Other', 'p53' };

temp = algs;

g1num = zeros( size( g1 ) );

for j = 1:size( temp, 2 )
    
    g1num( 1, strcmp( algs{ j }, g1 ) ) = j;
    
end

temp = drugs;

g2num = zeros( size( g2 ) );

for j = 1:size( temp, 2 )
    
    g2num( 1, strcmp( temp{ j }, g2 ) ) = j;
    
end

temp = unique( g3 );

g3num = zeros( size( g3 ) );

for j = 1:size( temp, 2 )
    
    g3num( 1, strcmp( temp{ j }, g3 ) ) = j;
    
end

[ p, F, df1, df2 ] = wanova ( y', g1num' )
[ p, F, df1, df2 ] = wanova ( y', g2num' )
[ p, F, df1, df2 ] = wanova ( y', g3num' )

% Figure 6.12

welchttest_two_classes = NaN( 21, 21 );

for j = 1:21
    
    for k = 1:21

            [ ~, p, ci, stats ] = ttest2( y( g2num == j ), y( g2num == k ), 'Vartype', 'unequal' );
            welchttest_two_classes( j, k ) = p;

    end

end

h = clustergram( log10( welchttest_two_classes.*210 ), 'RowLabels', drugtypes, 'ColumnLabels', drugtypes, 'Symmetric', ...
                 'false', 'Annotate', 'true', 'ColumnLabelsRotate', 25, 'RowLabelsRotate', 25, 'AnnotPrecision', 2 )
h.Colormap = turbo(100);

h.Dendrogram = 16;

for k = 1:21
    
    drugtypes{ k }
    nanmean( y( g2num == k ) )
    nanmedian( y( g2num == k ) )
    iqr( y( g2num == k ) )
    
end

% Figure 6.13

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y( g2num == k )';

end

neworder = [ 1, 9, 20, 6, 16, 10, 12, 8, 17, 19, 5, 2, 13, 7, 21, 15, 4, 14, 3, 11, 18 ];

gtemp = zeros( size( g2num ) );

for k = 1:21
    
    gtemp( g2num == neworder( k ) ) = k;
    
end

figure;

btemp = distributionPlot( tempcell( :, neworder ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( y, gtemp, 'PlotStyle', 'compact', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

hold off

xlim( [ 0, 22 ] )
xticks( 1:21 )
xticklabels( drugtypes( neworder ) )
xtickangle( 35 )
ylabel( 'Mean test ROC-AUC' )
ylim( [ 0.3, 1 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure C.18

figure;

t = tiledlayout( 4, 3 );

t.TileSpacing = 'compact';

for k = 1:12
   
    nexttile;
    histogram( mean_AUCs_test( :, k ) )
    xlabel( algorithms{ k } )
    
    ylim( [ 0, 80 ] )
    yticks( 0:20:80 )
    
    xlim( [ 0.4, 1  ] )
    xticks( 0.4:0.1:1 )
    
    grid on;
    box on;
    
end

title( t, 'Distribution of predictive performances across distinct algorithms' ) 

% Table B.10

for k = 1:12

    nanvar( mean_AUCs_test( :, k )  )

end

nanvar( mean_AUCs_test( :, 5 )  )/nanvar( mean_AUCs_test( :, 10 )  )

% Figure C.19

figure;

t2 = tiledlayout( 7, 3 );

t2.TileSpacing = 'compact';

for k = 1:21
   
    nexttile;
    
    temp = mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == k, 1:12 );
    histogram( reshape( temp, 1, size( temp, 1 ) * 12 ), 10 )
    xlabel( drugtypes{ k } )
    
    temp = yticks;
    
    if temp( end ) == 50
        
        yticks( 0:25:50 )
        
    end
    
    xlim( [ 0.4, 1 ] )
    xticks( 0.4:0.1:1 )
    
    grid on;
    box on;
    
end

title( t2, 'Distribution of predictive performance across distinct drug classes' ) 

% Table B.11

for k = 1:21

    k
    temp = mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == k, : );
    nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )

end

nanvar( reshape( mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 3, : ), 1, size( mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 3, : ), 1 ) * 12 ) )/ ...
        nanvar( reshape( mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 18, : ), 1, size( mean_AUCs_test( cell2mat( DrugOrder( :, 5 ) ) == 18, : ), 1 ) * 12 ) )

% Figure C.20

figure;

t3 = tiledlayout( 1, 2 );

t3.TileSpacing = 'compact';

ax = nexttile;

temp = mean_AUCs_test( Y == 1, 1:12 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \geq 600' )

set( ax, 'FontSize', 12 )

xlim( [ 0.4, 1 ] )
xticks( 0.4:0.1:1 )

grid on;
box on;

ax = nexttile;

temp = mean_AUCs_test( Y == 0, 1:12 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \leq 500' )

xlim( [ 0.4, 1 ] )
xticks( 0.4:0.1:1 )

grid on;
box on;

set( ax, 'FontSize', 12 )

title( t3, 'Distribution of predictive performance between drugs with large and small sample sizes' ) 

temp = mean_AUCs_test( Y == 1, : );
nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )

temp = mean_AUCs_test( Y == 0, : );
nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )

nanvar( reshape( mean_AUCs_test( Y == 0, : ), 1, size( mean_AUCs_test( Y == 0, : ), 1 ) * 12 ) )/...
nanvar( reshape( mean_AUCs_test( Y == 1, : ), 1, size( mean_AUCs_test( Y == 1, : ), 1 ) * 12 ) )

% Figure C.21

figure;

t5 = tiledlayout( 4, 3 ); 

t5.TileSpacing = 'compact';

for k = 1:12
   
    ax = nexttile;
    
    h1 = histogram( mean_AUCs_test( Y==1, k ) - nanmean( mean_AUCs_test( Y==1, k ) ) );
    h1.Normalization = 'probability';
    h1.BinWidth = 0.05;
    
    hold on;
    
    h2 = histogram( mean_AUCs_test( Y==0, k ) - nanmean( mean_AUCs_test( Y==0, k ) ), 'FaceColor', [0.6350 0.0780 0.1840] );
    h2.Normalization = 'probability';
    h2.BinWidth = 0.05;
    
    hold off
    
    xlabel( algorithms{ k } )
    
    if k == 12
        
        lgd = legend( ax, [ h1, h2 ], { 'Sample size \geq 600', 'Sample size \leq 500' }, 'Orientation', 'Horizontal', ...
            'NumColumns', 2, 'FontSize', 11 );
        lgd.Layout.Tile = 'North';
        
    end
    
    grid on;
    box on;
    xlim( [ -0.4, 0.4 ] )
    xticks( -0.4:0.1:0.4 )
    
end

title( t5, 'Distribution of predictive performances across distinct algorithms for large and small sample sizes' ) 

% Table B.12

for k = 1:12

    k
    nanvar( mean_AUCs_test( Y==1, k )  )
    nanvar( mean_AUCs_test( Y==0, k )  )

end

temp = [];

for k = 1:12
    
    temp = [ temp, nanvar( mean_AUCs_test( Y==1, k )  ), nanvar( mean_AUCs_test( Y==0, k )  ) ];
    
end

max( temp )
min( temp )
max( temp )/min( temp )

% Table 6.4

y2 = mean_AUCs_test;
y2 = reshape( y2, 1, 12*265 );

[ p, tbl, stats, terms ] = anovan( y2, { g1 g3 }, 'model', 'linear', 'varnames', { 'Alg','Num' } ); %#ok<*ASGLU>

% Figure 6.14

figure;
[ c, m, h, gnames ] = multcompare( stats, 'Dimension', 1 )

barcoords = [ 0.710650184134791,0.734904678683932; 0.656428370862053,0.680682865411194; ...
              0.69040262040705,0.714657114956191; 0.702141865137518,0.726396359686659; ...
              0.702830648263959,0.7270851428131; 0.698337332658843,0.722591827207984; ...
              0.695789422915886,0.720043917465027; 0.692520948935181,0.716775443484322; ...
              0.70332630604575,0.727580800594891; 0.701084251086119,0.72533874563526; ...
              0.700962954236194,0.725217448785335; 0.689077661672396,0.713332156221537 ];

figure;
[ c2, m2, h2, gnames2 ] = multcompare( stats, 'Dimension', 2 )

barcoords2 = [ 0.682393092118175,0.690102895946736; 0.724743692995028,0.732453496823589 ];

figure;

t5 = tiledlayout( 2, 1 ); 

t5.TileSpacing = 'compact';

ax = nexttile;

hold on;

for k = 1:12

    plot( barcoords( k, : ), [ k, k ], 'Color', [ 0, 0.4470, 0.7410 ], 'LineWidth', 2 );
    p1 = plot( m( k, 1 ), 13-k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [ 0, 0.4470, 0.7410 ], 'MarkerFaceColor', [ 0, 0.4470, 0.7410 ] )

end

hold off

ylim( [ 0, 13 ] )
yticks( 1:12 )
yticklabels( algorithms( 12:-1:1 ) )
ytickangle( 25 )

xlim( [ 0.65 0.75 ] )
xticks( 0.65:0.025:0.75)
xticklabels( { '0.65', '', '0.7', '', '0.75' } )

title( 'Algorithms' )

grid on
box on

set( gca, 'FontSize', 14 )


ax = nexttile;

hold on

for k = 1:2

    plot( barcoords2( k, : ), [ k, k ], 'Color', [ 0, 0.4470, 0.7410 ], 'LineWidth', 2 );
    p1 = plot( m2( k, 1 ), k, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [ 0, 0.4470, 0.7410 ], 'MarkerFaceColor', [ 0, 0.4470, 0.7410 ] )

end

hold off

ylim( [ 0, 3 ] )
yticks( 1:2 )
yticklabels( { 'Sample size \leq 500', 'Sample size \geq 600' } )
ytickangle( 25 )

xlim( [ 0.65 0.75 ] )
xticks( 0.65:0.025:0.75)
xticklabels( { '0.65', '', '0.7', '', '0.75' } )

title( 'Sample size' )

grid on
box on

set( gca, 'FontSize', 14 )

ylabel( t5, '\bf{Factor levels}', 'FontSize', 16 )
xlabel( t5, '\bf{Mean test ROC-AUC}', 'FontSize', 16 )

% Figure 6.15

figure;
b = distributionPlot( mean_AUCs_test, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 )

hold on

boxplot( mean_AUCs_test, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

hold off

xtickangle( 35 )

grid on;
box on;

ylabel( 'Mean test ROC-AUCs' )
ylim( [ 0, 1 ] )
xlim( [ 0, 13 ] )

xticklabels( algorithms )

set( gca, 'FontSize', 14 )

% Figure 6.16

tempcell = cell( 1, 2 );

tempcell{ :, 1 } = reshape( mean_AUCs_test( Y == 1, : ), 1, 217 * 12 );
tempcell{ :, 2 } = reshape( mean_AUCs_test( Y == 0, : ), 1, 48 * 12 );

temp = [ tempcell{ :, 1 }, tempcell{ :, 2 } ];
tempind = [ ones( 1, 2604 ), 2 * ones( 1, 576 ) ];

figure;

btemp = distributionPlot( tempcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( temp, tempind,  'PlotStyle', 'compact', 'Widths', 0.15, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

xticks( 1:2 )
xticklabels( { 'Sample size \geq 600', 'Sample size \leq 500' } )
ylabel( 'Mean test ROC-AUC' )
ylim( [ 0.35, 1 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )

xlim( [ 0, 3 ] )
ylim( [ 0.3, 1 ] )

nanmedian( tempcell{ :, 1 } ) - nanmedian( tempcell{ :, 2 } )

% Table B.13

for k = 1:12
    
    k
    
    m3 = nanmean( mean_AUCs_test( Y==0, k )  )
    m4 = nanmedian( mean_AUCs_test( Y==0, k )  )
    
    m1 = nanmean( mean_AUCs_test( Y==1, k )  )
    m2 = nanmedian( mean_AUCs_test( Y==1, k )  )
    
    m1 - m3
    m2 - m4
    
end

% Table 6.5

[ p, F, df1, df2 ] = wanova ( y2', g1num' )
[ p, F, df1, df2 ] = wanova ( y2', g2num' )
[ p, F, df1, df2 ] = wanova ( y2', g3num' )

% Figure 6.17

welchttest_one = NaN( 21, 21 );

for j = 1:21
    
    for k = 1:21

            [ ~, p, ci, stats ] = ttest2( y2( g2num == j ), y2( g2num == k ), 'Vartype', 'unequal' );
            welchttest_one( j, k ) = p;

    end

end


h = clustergram( log10( welchttest_one.*210 ), 'RowLabels', drugtypes, 'ColumnLabels', drugtypes, 'Symmetric', ...
                 'false', 'Annotate', 'true', 'ColumnLabelsRotate', 25, 'RowLabelsRotate', 25, 'DisplayRatio', 0.2 )
h.Colormap = turbo( 100 );

h.Dendrogram = 14;

% Figure 6.18

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y2( g2num == k )';

end

neworder = [ 6, 9, 20, 12, 2, 5, 8, 17, 19, 10, 16, 7, 21, 15, 1, 14, 13, 4, 11, 3, 18 ];

gtemp = zeros( size( g2num ) );

for k = 1:21
    
    gtemp( g2num == neworder( k ) ) = k;
    
end

figure;

btemp = distributionPlot( tempcell( :, neworder ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( y2, gtemp, 'PlotStyle', 'compact', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

hold off

xlim( [ 0, 22 ] )
xticks( 1:21 )
xticklabels( drugtypes( neworder ) )
xtickangle( 35 )
ylabel( 'Mean test ROC-AUC' )
ylim( [ 0.3, 1 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )


for k = 1:21
    
    drugtypes{ neworder( k ) }
    nanmean( tempcell{ 1, neworder( k ) } )
    nanmedian( tempcell{ 1, neworder( k ) } )
    iqr( tempcell{ 1, neworder( k ) } )
    nanvar( tempcell{ 1, neworder( k ) } )
    
end

% Rankings of drugs wrt. average model performance

% 1. Two-step models

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y( g2num == k )';

end

means2 = zeros( 21, 1 );
medians2 = zeros( 21, 1 );

for k = 1:21
    
    means2( k ) = nanmean( tempcell{ 1, k } );
    medians2( k ) = nanmedian( tempcell{ 1, k } );
    
end

[ ~, indmean2 ] = sort( means2, 'descend' );
[ ~, indmedian2 ] = sort( medians2, 'descend' );


% One-step models

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y2( g2num == k )';

end

means1 = zeros( 21, 1 );
medians1 = zeros( 21, 1 );

for k = 1:21
    
    means1( k ) = nanmean( tempcell{ 1, k } );
    medians1( k ) = nanmedian( tempcell{ 1, k } );
    
end

[ ~, indmean1 ] = sort( means1, 'descend' );
[ ~, indmedian1 ] = sort( medians1, 'descend' );

% Comparison

comp_mean = zeros( 21, 3 );

for k = 1:21
    
    comp_mean( k, 1 ) = find( indmean2 == k );
    comp_mean( k, 2 ) = find( indmean1 == k );
    comp_mean( k, 3 ) = comp_mean( k, 1 ) - comp_mean( k, 2 );
    
end

% Figure C.22

differences = results_mean_AUCs_test( :, 7:18 ) - mean_AUCs_test;

figure;

t = tiledlayout( 4, 3 );

t.TileSpacing = 'compact';

for k = 1:12
   
    nexttile;
    histogram( differences( :, k ) )
    xlabel( algorithms{ k } )
    
    ylim( [ 0, 100 ] )
    yticks( 0:50:100 )
    
    grid on;
    box on;
    
    xlim( [ -0.3, 0.2 ] )
    xticks( -0.3:0.1:0.2)
    
end

title( t, 'Distribution of differences in predictive performance across distinct algorithms' )

% Table B.14

for k = 1:12

    k
    nanvar( differences( :, k )  )

end

nanvar( differences( :, 7 )  )\nanvar( differences( :, 10 )  )

% Figure C.23

figure;

t2 = tiledlayout( 7, 3 );

t2.TileSpacing = 'compact';

ylimits = { [ 0, 30 ], [0, 50], [0, 40], [0, 100], [0, 30], [0,40], ...
    [0, 30], [0, 150], [0, 200], [0, 30], [0, 15], [0, 40], ...
    [0, 75], [0, 40], [0, 20], [0, 10], [0, 75], [0, 15], ...
    [0, 30], [0, 400], [0, 15] };

xlimits = { [ -0.3, 0.1 ], [-0.1, 0.1], [-0.3, 0.1], [-0.15, 0.1], [-0.15, 0.1], [-0.3, 0.1], ...
    [-0.1, 0.1], [-0.2, 0.2], [-0.3, 0.1], [-0.2, 0.1], [-0.1, 0.05], [-0.1, 0.1], ...
    [-0.3, 0.1], [-0.1, 0.15], [-0.06, 0.06], [-0.08, 0.06], [-0.1, 0.1], [-0.1, 0.1], ...
    [-0.1, 0.15], [-0.2, 0.1], [-0.25, 0.1] };

xtls = { -0.3:0.1:0.1, -0.1:0.05:0.1, -0.3:0.1:0.1, -0.15:0.05:0.1, -0.15:0.05:0.1, -0.3:0.1:0.1, ...
    -0.1:0.05:0.1, -0.2:0.1:0.2, -0.3:0.1:0.1, -0.2:0.1:0.1, -0.1:0.05:0.05, -0.1:0.05:0.1, ...
    -0.3:0.1:0.1, -0.1:0.05:0.15, -0.06:0.02:0.06, -0.08:0.02:0.06, -0.1:0.05:0.1, -0.1:0.05:0.1, ...
    -0.1:0.05:0.15, -0.2:0.1:0.1, -0.25:0.05:0.1 };

for k = 1:21

    nexttile;

    temp = differences( cell2mat( DrugOrder( :, 5 ) ) == k, 1:12 );
    histogram( reshape( temp, 1, size( temp, 1 ) * 12 ), 10 )
    xlabel( drugtypes{ k } )

    temp = yticks;

    if temp( end ) == 50

        yticks( 0:25:50 )

    end

    grid on;
    box on;

    xlim( xlimits{ k } )
    ylim( ylimits{ k } )
    xticks( xtls{ k } )

end

title( t2, 'Distribution of differences in predictive performance across distinct drug classes' )

% Table B.15

for k = 1:21

    temp = differences( cell2mat( DrugOrder( :, 5 ) ) == k, : );
    nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )

end

nanvar( reshape( differences( cell2mat( DrugOrder( :, 5 ) ) == 11, : ), 1, size( differences( cell2mat( DrugOrder( :, 5 ) ) == 11, : ), 1 ) * 12 ) )\...
nanvar( reshape( differences( cell2mat( DrugOrder( :, 5 ) ) == 6, : ), 1, size( differences( cell2mat( DrugOrder( :, 5 ) ) == 6, : ), 1 ) * 12 ) )

nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) )

% Figure C.24

figure;

t3 = tiledlayout( 1, 2 );

t3.TileSpacing = 'compact';

nexttile;

temp = differences( Y == 1, 1:12 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \geq 600' )
set(gca, 'FontSize',12)

grid on;
box on;

xlim( [ -0.3, 0.2 ] )
xticks( -0.3:0.1:0.2 )

nexttile;

temp = differences( Y == 0, 1:12 );
histogram( reshape( temp, 1, size( temp, 1 ) * 12 ) )
xlabel( 'Sample size \leq 500' )
set(gca, 'FontSize',12)

grid on;
box on;

xlim( [ -0.3, 0.2 ] )
xticks( -0.3:0.1:0.2 )

title( t3, 'Distribution of differences in predictive performance between drugs with high and low sample sizes' ) 

temp = differences( Y == 1, : );
temp1 = nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) );

temp = differences( Y == 0, : ); % low sample size
temp2 = nanvar( reshape( temp, 1, size( temp, 1 ) * 12 ) );

temp2/temp1

% Table 6.6

y3 = reshape( differences, 1, 12*265 );

[ p, F, df1, df2 ] = wanova ( y3', g1num' )
[ p, F, df1, df2 ] = wanova ( y3', g2num' )
[ p, F, df1, df2 ] = wanova ( y3', g3num' )

% Figure 6.19

welchttest_diff = NaN( 12, 12 );

for j = 1:12
    
    for k = 1:12

        [ ~, p, ci, stats ] = ttest2( y3( g1num == j ), y3( g1num == k ), 'Vartype', 'unequal' );
        welchttest_diff( j, k ) = p;

    end

end

h = clustergram( log10( welchttest_diff.*66 ), 'RowLabels', algorithms, 'ColumnLabels', algorithms, 'Symmetric', ...
                 'false', 'Annotate', 'true', 'ColumnLabelsRotate', 25, 'RowLabelsRotate', 25 )
h.Colormap = turbo( 100 );
h.Dendrogram = 5.8;

% Figure 6.20

neworder = [ 11, 1, 5, 6, 7, 2, 3, 9, 8, 4, 10, 12 ];

figure;

btemp = distributionPlot( differences( :, neworder ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( differences( :, neworder ), 'PlotStyle', 'compact', 'Widths', 0.15, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target' )

ylim( [ -0.3, 0.2 ] )
xlim( [ 0, 13 ] )
xticks( 1:12 )
xticklabels( algorithms( neworder ) )
xtickangle( 35 )
ylabel( 'Mean test ROC-AUC' )

grid on;
box on;

set( gca, 'FontSize', 14 )


for k = 1:12
    
    algorithms{ neworder( k ) }
    nanmean( differences( :, neworder( k ) ) )
    nanmedian( differences( :, neworder( k ) ) )
    
end

% Figure 6.21

welchttest_diff_classes = NaN( 21, 21 );

for j = 1:21
    
    for k = 1:21

            [ ~, p, ci, stats ] = ttest2( y3( g2num == j ), y3( g2num == k ), 'Vartype', 'unequal' );
            welchttest_diff_classes( j, k ) = p;

    end

end

h = clustergram( log10( welchttest_diff_classes.*210 ), 'RowLabels', drugtypes, 'ColumnLabels', drugtypes, 'Symmetric', ...
                 'false', 'Annotate', 'true', 'ColumnLabelsRotate', 25, 'RowLabelsRotate', 25 )
h.Colormap = turbo( 100 );
h.Dendrogram = 14.7;

% Figure 6.22

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = y3( g2num == k )';

end

neworder = [ 1, 13, 9, 4, 20, 8, 6, 21, 16, 10, 3, 12, 5, 2, 17, 11, 7, 18, 19, 15, 14 ];

gtemp = zeros( size( g2num ) );

for k = 1:21
    
    gtemp( g2num == neworder( k ) ) = k;
    
end

figure;

btemp = distributionPlot( tempcell( :, neworder ), 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( y3, gtemp, 'PlotStyle', 'compact', 'Widths', 0.03, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

hold off

xlim( [ 0, 22 ] )
xticks( 1:21 )
xticklabels( drugtypes( neworder ) )
xtickangle( 35 )
ylabel( 'Mean test ROC-AUC' )
ylim( [ -0.35, 0.15 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )


for k = 1:21
    
    drugtypes( neworder( k ) )
    nanmedian( tempcell{ 1, neworder( k ) } )
    nanmean( tempcell{ 1, neworder( k ) } )
    nanmin( tempcell{ 1, neworder( k ) } )
    iqr( tempcell{ 1, neworder( k ) } )
    
end

% Figure 6.23

tempcell = cell( 1, 2 );

tempcell{ :, 1 } = y3( g3num == 1 );
tempcell{ :, 2 } = y3( g3num == 2 );

figure;

btemp = distributionPlot( tempcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( y3, g3num, 'PlotStyle', 'compact', 'Widths', 0.15, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

xticks( 1:2 )
xticklabels( { 'Sample size \geq 600', 'Sample size \leq 500' } )
ylabel( 'Mean test ROC-AUC' )
ylim( [ -0.35, 0.2 ] )
xlim( [ 0, 3 ] )

grid on;
box on;

set( gca, 'FontSize', 14 )


for k = 1:2
   
    k
    nanvar( tempcell{ :, k } )
    iqr( tempcell{ :, k } )
    nanmean( tempcell{ :, k } )
    
end

% Figure 6.24

impscores_alg_two = NaN( 11, 265, 6 );

for k = 1:265 % run through the drugs
    
    load( [ 'res_', num2str( k ), '.mat' ],'coeffs' )
    
    for l = 2:12 % run through the algorithms
        
        index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;
        
        temp = NaN( 10, 6 );
        
        for m = 1:10 % run through the folds
            
            if l <= 9
                
                temptemp = coeffs{ m, l };
                temptemp = temptemp( 2:length( temptemp ) );
                
            else
                
                temptemp = coeffs{ m, l };
                
            end
            
            temp( m, coeffs{ m, 1 } == 1 ) = temptemp;
            
        end
        
        % convert into absolute values and normalize per row, then
        % compute the median
        
        temp = abs( temp );
        temp = temp./nanmax( temp, [], 2 );
        temp = nanmedian( temp );
        temp( ~index ) = NaN;
        
        impscores_alg_two( l - 1, k, : ) = temp;
        
    end
    
end

impscores_all_alg_two = cell( 11, 6 );

for k = 1:265 % run through the drugs
    
    load( [ 'res_', num2str( k ), '.mat' ],'coeffs' );
    
    for l = 2:12 % run through the algorithms
        
        index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;
        
        temp = cell( 1, 6 );
        
        for m = 1:10 % run through the folds
            
            if l <= 9 % remove constant
                
                temptemp = coeffs{ m, l };
                temptemp = abs( temptemp( 2:length( temptemp ) ) );
                
                if nanmax( temptemp )
                
                    temptemp = temptemp./nanmax( temptemp );

                end
                
            else
                
                temptemp = abs( coeffs{ m, l } );
                
                if nanmax( temptemp )
                    
                    temptemp = temptemp./nanmax( temptemp );
                
                end
                
            end
            
            % use all weights for any stable data type
            
            for n = 1:6
                
                temptempindex = coeffs{ m , 1 };
                
                if temptempindex( n ) && index( n )

                    temp{ 1, n } = [ temp{ 1, n }; temptemp( sum( temptempindex( 1:n ) ) ) ];

                end
                
            end
            
        end
        
        for n = 1:6
        
            impscores_all_alg_two{ l - 1, n } = [ impscores_all_alg_two{ l - 1, n }; temp{ 1, n } ];
        
        end
        
    end
    
end

impscores_alg_one = NaN( 11, 265, 6 );

for k = 1:265 % run through the drugs
    
    load( [ 'res_', num2str( k ), '.mat' ],'coeffs' )
    
    for l = 3:13 % run through the algorithms
        
        index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;
        
        temp = NaN( 10, 6 );
        
        for m = 1:10 % run through the folds
            
            if l <= 10 % remove the constant coefficient
                
                temptemp = coeffs{ m, l };
                temptemp = abs( temptemp( 2:length( temptemp ) ) );
                
            else
                
                temptemp = abs( coeffs{ m, l } );
                
            end
            
            % Only use the maximum weight per data type
            
            for n = 1:6
                
                temptempindex = coeffs{ m , 1 };
                
                if temptempindex( n )
                    
                    tempindex = coeffs{ m, 2 } == n;
                    temp( m, n ) = nanmax( temptemp( tempindex ) );
                    
                end
                
            end
            
        end
        
        % Normalise per row, then compute the median
        
        temp = temp./nanmax( temp, [], 2 );
        temp = nanmedian( temp );
        temp( ~index ) = NaN;
        
        impscores_alg_one( l - 2, k, : ) = temp;
        
    end
    
end

impscores_max_alg_one = cell( 11, 6 );

for k = 1:265 % run through the drugs
    
    load( [ 'res_', num2str( k ), '.mat' ],'coeffs' )
    
    for l = 3:13 % run through the algorithms
        
        index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5; % only data types with non-trivial models in at least five folds
        
        temp = cell( 1, 6 ); 
        
        for m = 1:10 % run through the folds
            
            if l <= 10 % remove constant
                
                temptemp = coeffs{ m, l };
                temptemp = abs( temptemp( 2:length( temptemp ) ) );
                
                if nanmax( temptemp )
                
                    temptemp = temptemp./nanmax( temptemp );

                end
                
            else
                
                temptemp = abs( coeffs{ m, l } );
                
                if nanmax( temptemp )
                    
                    temptemp = temptemp./nanmax( temptemp );
                
                end
                
            end
            
            % use the maximum weights for any stable data type
            
            for n = 1:6
                
                temptempindex = coeffs{ m , 1 };
                
                if temptempindex( n ) && index( n )
                    
                    tempindex = coeffs{ m, 2 } == n;
                    
                    temp{ 1, n } = [ temp{ 1, n }; nanmax( temptemp( tempindex ) ) ];
           
                end
                
            end
            
        end
        
        for n = 1:6
        
            impscores_max_alg_one{ l - 2, n } = [ impscores_max_alg_one{ l - 2, n }; temp{ 1, n } ];
        
        end
        
    end
    
end


figure;
t = tiledlayout( 6, 1 );
t.TileSpacing = 'compact';
t.Padding = 'compact';

for k = 1:6
    
    nexttile;
    
    temp = cell( 1, 11 );
    temptwo = cell( 1, 11 );
    temp2 = [];
    temp3 = [];
    g = [];
    g2 = [];
    
    if k >= 6
        
        for l = 1:11
            
            temptemp = impscores_max_alg_one{ l, k };
            temptemptwo = impscores_all_alg_two{ l, k };
            
            temp{ 1, l } = temptemp( ~isnan( temptemp ) );
            temptwo{ 1, l } = temptemptwo( ~isnan( temptemptwo ) );
            
            temp2 = [ temp2, temptemp ];
            g = [ g, repmat( algorithms(1, l + 1 ), 1, length( temptemp ) ) ];
            
            temp3 = [ temp3, temptemptwo ];
            g2 = [ g2, repmat( algorithms(1, l + 1 ), 1, length( temptemptwo ) ) ];
            
        end
 
        b = distributionPlot( temptwo, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1] );
        
        hold on
        
        b2 = distributionPlot( temp, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );
        
        b3 = boxplot( temp3, g2, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle','target', 'Positions', 0.95:1:10.95 );
        b4 = boxplot( temp2, g, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle','target', 'Positions', 1.05:1:11.05 );
        
        hold off
        
        xtickangle( 25 )
        
        xlim( [ 0, 12 ] )

    else
        
        for l = 1:11
            
            temptemp = impscores_max_alg_one{ l, k };
            temptemptwo = impscores_all_alg_two{ l, k };
            
            temp{ 1, l } = temptemp( ~isnan( temptemp ) );
            temptwo{ 1, l } = temptemptwo( ~isnan( temptemptwo ) );
            
            temp2 = [ temp2, temptemp ];
            g = [ g, l * ones( 1, length( temptemp ) ) ];
            
            temp3 = [ temp3, temptemptwo ];
            g2 = [ g2, l * ones( 1, length( temptemptwo ) ) ];
            
        end
        
        b = distributionPlot( temptwo, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1]  );
        
        hold on
        
        b2 = distributionPlot( temp, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

        b3 = boxplot( temp3, g2, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle','target', 'Positions', 0.95:1:10.95 );
        b4 = boxplot( temp2, g, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle','target', 'Positions', 1.05:1:11.05 );
        
        
        hold off
        
        xlim( [ 0, 12 ] )
        
        set( gca, 'XTickLabel', {' '} );
    
    end
    
    grid on;
    box on;
    
    ylim( [ -0.5, 1.5 ] )
    yticks( -0.5:0.25:1.5 )
    yticklabels( { '-0.5', '', '0', '', '0.5', '', '1', '', '1.5' } )
    title( datatypes{ k }, 'FontSize', 12 )
    
    set( gca, 'FontSize', 12 )
    
    ylabel( t, 'Normalised absolute weights or importance scores', 'FontSize', 14 )
    
end

lgd = legend( [ b{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Two-step models', 'One-step models' }, 'Location', 'southeast', 'Orientation', 'horizontal' );
legend( 'boxon' )

title( t, ' \newline ' )

% Figure C.25

p1 = zeros( 6, 11 ); % between one- and two-step models

for k = 1:6
    
    for l = 1:11
        
        x = NaN( max( length( impscores_all_alg_two{ l, k } ), length( impscores_max_alg_one{ l, k } ) ), 2 );
        x( 1:length( impscores_all_alg_two{ l, k } ), 1 ) = impscores_all_alg_two{ l, k };
        x( 1:length( impscores_max_alg_one{ l, k } ), 2 ) = impscores_max_alg_one{ l, k };
        
       p = kruskalwallis( x );
       p1( k, l ) = p;
        
    end
    
end

p1 = p1*66;
p1 = p1';

figure;
h = heatmap( datatypes, algorithms( 1, 2:12 ), log10( p1 ), 'Colormap',jet, 'MissingDataLabel', 'None found', 'MissingDataColor', [ 1 1 1 ], ...
    'ColorbarVisible','on', 'CellLabelFormat', '%0.3g', 'GridVisible', 'on' );

set( gca, 'FontSize', 12 )

% Figure 6.25

impscores_med2 = zeros( 11, 6 );

for k = 1:11 
    
    for j = 1:6
        
        temp = impscores_all_alg_two{ k, j };
        impscores_med2( k, j ) = nanmedian( temp );
        
    end
    
end

impscores_med1 = zeros( 11, 6 );

for k = 1:11 
    
    for j = 1:6
        
        temp = impscores_max_alg_one{ k, j };
        impscores_med1( k, j ) = nanmedian( temp );
        
    end
    
end

Algorithms = { 'Neural network', 'LASSO-regularised linear regression', 'Elastic net-regularised linear regression', 'Ridge-regularised linear regression', ...
               'LASSO-regularised logistic regression', 'Elastic net-regularised logistic regression', 'Ridge-regularised logistic regression', ...
           	   'Ridge-regularised SVM', 'LASSO-regularised SVM', 'Bagged decision tree ensemble', 'Boosted decision tree ensemble', 'Random forest' };

figure;

t = tiledlayout( 2, 2 );
t.TileSpacing = 'compact';

for k = 1:3
    
    nexttile;
    
    temp = [ impscores_med2( k, : ); impscores_med1( k, : ) ]; 
    spider_plot( temp, 'AxesLabels', datatypes, 'FillOption', 'on', 'FillTransparency', ...
        0.2, 'LabelFontSize', 10, 'AxesLimits', [ 0, 0, 0, 0, 0, 0; 1, 1, 1, 1, 1, 1 ], ...
        'Color', [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840 ] )
    title( Algorithms{ 1, k + 1 }, 'FontSize', 14 )
    
    if k == 3
        
        legend( { 'Two-step models', 'One-step models' }, 'FontSize', 14 )
        
    end
    
end

% Figure 6.26

figure;

t = tiledlayout( 2, 2 );
t.TileSpacing = 'compact';

for k = 4:6
    
    nexttile;
    
    temp = [ impscores_med2( k, : ); impscores_med1( k, : ) ]; 
    spider_plot( temp, 'AxesLabels', datatypes, 'FillOption', 'on', 'FillTransparency', ...
        0.2, 'LabelFontSize', 10, 'AxesLimits', [ 0, 0, 0, 0, 0, 0; 1, 1, 1, 1, 1, 1 ], ...
        'Color', [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840 ] )
    title( Algorithms{ 1, k + 1 }, 'FontSize', 14 )
    
    if k == 6
        
        legend( { 'Two-step models', 'One-step models' }, 'FontSize', 14 )
        
    end
    
end

% Figure 6.27

figure;

t = tiledlayout( 1, 2 );
t.TileSpacing = 'compact';

for k = 7:8
    
    nexttile;
    
    temp = [ impscores_med2( k, : ); impscores_med1( k, : ) ]; 
    spider_plot( temp, 'AxesLabels', datatypes, 'FillOption', 'on', 'FillTransparency', ...
        0.2, 'LabelFontSize', 10, 'AxesLimits', [ 0, 0, 0, 0, 0, 0; 1, 1, 1, 1, 1, 1 ], ...
        'Color', [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840 ] )
    title( Algorithms{ 1, k + 1 }, 'FontSize', 14 )
    
    if k == 8
        
        legend( { 'Two-step models', 'One-step models' }, 'FontSize', 14 )
        
    end
    
end

% Figure 6.28

figure;

t = tiledlayout( 2, 2 );
t.TileSpacing = 'compact';

for k = 9:11
    
    nexttile;
    
    temp = [ impscores_med2( k, : ); impscores_med1( k, : ) ]; 
    spider_plot( temp, 'AxesLabels', datatypes, 'FillOption', 'on', 'FillTransparency', ...
        0.2, 'LabelFontSize', 10, 'AxesLimits', [ 0, 0, 0, 0, 0, 0; 1, 1, 1, 1, 1, 1 ], ...
        'Color', [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840 ] )
    title( Algorithms{ 1, k + 1 }, 'FontSize', 14 )
    
    if k == 11
        
        legend( { 'Two-step models', 'One-step models' }, 'FontSize', 14 )
        
    end
    
end

% Figure 6.29

diffs = results_mean_AUCs_test( :, 7:18 ) - mean_AUCs_test;
[ ind1, ind2 ] = find( diffs > 0 ); % 979 out of 3180 cases

impscores_tbetter_two = zeros( 979, 6 );
impscores_tbetter_one = zeros( 979, 6 );

for k = 1:979
    
    % Exclude neural networks, because we do not compute weights for those
    if ind2( k ) == 1
        
        impscores_tbetter_two( k, : ) = NaN( 1, 6 );
        impscores_tbetter_one( k, : ) = NaN( 1, 6 );
    
    else
        
        impscores_tbetter_two( k, : ) = squeeze( impscores_alg_two( ind2( k ) - 1, ind1( k ), : ) );
        impscores_tbetter_one( k, : ) = squeeze( impscores_alg_one( ind2( k ) - 1, ind1( k ), : ) );
    
    end
    
end

% Remove the first 112 rows, which pertain to neural networks

impscores_tbetter_two = impscores_tbetter_two( 113:979, : );
impscores_tbetter_one = impscores_tbetter_one( 113:979, : );

% compute p-values

ps = zeros( 1, 6 );

for k = 1:6
    
    p = kruskalwallis( [ impscores_tbetter_two( :, k ), impscores_tbetter_one( :, k ) ] );
    ps( 1, k ) = p;
    
end

% Repeat for the complement set of models 

[ ind1, ind2 ] = find( diffs <= 0 ); % 2201 out of 3180 cases

impscores_tworse_two = zeros( 2201, 6 );
impscores_tworse_one = zeros( 2201, 6 );

for k = 1:2201
    
    % exclude neural networks, because we do not compute weights for those
    if ind2( k ) == 1
        
        impscores_tworse_two( k, : ) = NaN( 1, 6 );
        impscores_tworse_one( k, : ) = NaN( 1, 6 );
    
    else
        
        impscores_tworse_two( k, : ) = squeeze( impscores_alg_two( ind2( k ) - 1, ind1( k ), : ) );
        impscores_tworse_one( k, : ) = squeeze( impscores_alg_one( ind2( k ) - 1, ind1( k ), : ) );
    
    end
    
end

% p-values

ps = zeros( 1, 6 );

for k = 1:6
    
    p = kruskalwallis( [ impscores_tworse_two( :, k ), impscores_tworse_one( :, k ) ] );
    ps( 1, k ) = p;
    
end

pstext1 = { '2\cdot10^{-10}', '0.0001', '6\cdot10^{-10}', '0.0764', '2\cdot10^{-14}', '1\cdot10^{-7}' };
pstext2 = { '2\cdot10^{-54}', '1\cdot10^{-14}', '2\cdot10^{-18}', '1\cdot10^{-8}', '5\cdot10^{-49}', '4\cdot10^{-31}' };

figure;

t = tiledlayout( 2, 1 );
t.TileSpacing = 'compact';

nexttile;

b = distributionPlot( impscores_tbetter_two, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1]  );

hold on

b2 = distributionPlot( impscores_tbetter_one, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

b3 = boxplot( impscores_tbetter_two, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle','target', 'Positions', 0.95:1:5.95 );
b4 = boxplot( impscores_tbetter_one, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle','target', 'Positions', 1.05:1:6.05 );

for k = 1:6
    
    plot( [ k - 0.25, k + 0.25 ], [ -0.3, -0.3 ], 'k' )
    plot( [ k - 0.25, k - 0.25 ], [ -0.25, -0.3 ], 'k' )
    plot( [ k + 0.25, k + 0.25 ], [ -0.25, -0.3 ], 'k' )
    text( k, -0.5, [ 'p = ', pstext1( k ) ], 'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'tex' )
    
end

hold off

xticks( 1:6 ) 

grid on; box on;
ylim( [ -0.75, 1.5 ] )
yticks( 0:0.5:1 )
xlim( [ 0, 7 ] )
xticks( 1:6 )
xticklabels( { '', '', '', '', '', '' } )
set( gca, 'FontSize', 14 )

title( '\Delta_{perf} > 0' )

nexttile;

b = distributionPlot( impscores_tworse_two, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1]  );

hold on

b2 = distributionPlot( impscores_tworse_one, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

b3 = boxplot( impscores_tworse_two, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle','target', 'Positions', 0.95:1:5.95 );
b4 = boxplot( impscores_tworse_one, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle','target', 'Positions', 1.05:1:6.05 );

for k = 1:6
    
    plot( [ k - 0.25, k + 0.25 ], [ -0.3, -0.3 ], 'k' )
    plot( [ k - 0.25, k - 0.25 ], [ -0.25, -0.3 ], 'k' )
    plot( [ k + 0.25, k + 0.25 ], [ -0.25, -0.3 ], 'k' )
    text( k, -0.5, [ 'p = ', pstext2( k ) ], 'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'tex' )
    
end

hold off

xticks( 1:6 ) 

grid on; box on;

ylim( [ -0.75, 1.5 ] )
yticks( 0:0.5:1 )

xlim( [ 0, 7 ] )
xticks( 1:6 )
xticklabels( datatypes )
xtickangle( 25 )
set( gca, 'FontSize', 14 )

title( '\Delta_{perf} < 0' )

xlabel( t, '\bf{Data types}', 'FontSize', 14  )
ylabel( t, '\bf{Normalised weights and importance scores}', 'FontSize', 14 )
title( t, '' )

lgd = legend( [ b{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Two-step models', 'One-step models' }, 'Orientation', 'horizontal' );

%% Chapter 7 

% Retrieve data from the paper "SYSTEMATIC ASSESSMENT OF ANALYTICAL METHODS 
% FOR DRUG SENSITIVITY PREDICTION FROM CANCER CELL LINE DATA"
% by Jang et al., 2014.

% The performance results pertaining to classification models are stored in
% the file 'SangerperformanceTableAUCnew'. Import it as a table.

drugorder_jang = unique(SangerperformanceTableAUCnew.DrugName); 
drugorder = DrugOrder( :, 1 );

% Remove the signs '.', '-','(', ')' and '_' from both orderings

drugorder_jang = strrep( drugorder_jang, '.', '' );
drugorder_jang = strrep( drugorder_jang, ' ', '' );

drugorder = strrep( drugorder, '-', '' );
drugorder = strrep( drugorder, ' ', '' );
drugorder = strrep( drugorder, '_', '' );
drugorder = strrep( drugorder, '(', '' );
drugorder = strrep( drugorder, ')', '' );

% Match compound names between both files

indices = NaN( 138, 1 );

for k = 1:138
    
    temp = strcmpi( drugorder_jang( k ), drugorder );
    
    if sum( temp ) == 1
       
        indices( k ) = find( temp == 1 );
        
    end
    
end

sum( isnan( indices ) )

% 15 drugs could not be found - search for alternate names that do match
% and manually add them, taking into account duplicate drugs within the
% GDSC

find( isnan( indices ) )
drugorder_jang(  isnan( indices )  )

% "ABT263" = Navitoclax (177)
% "ABT888" = Veliparib (184)
% "AZD0530" = Saracatinib (14)
% "AZD2281" = Olaparib, Olaparib_2 (183,260)
% "AZD6244" = Selumetinib (217,261)
% "BID1870"
% "BIBW2992" = Afatinib (197,257)
% "Metformin"
%"NVPBEZ235" = BEZ235 (212)
%"NVPTAE684" = TAE684 (12)
% "OSI906" = Lisitinib (75)
% "PF02341066" = Crizotinib (13)
% "WO2009093972" = AZD6482_2 (54, 218)
% "X17AAG" = 17AAG (192)
% "X681640" = 681640 (205)

indices( 3, 1 ) = 177;
indices( 4, 1 ) = 184;
indices( 14, 1 ) = 14;
indices( 15, 1 ) = 183;
indices( 15, 2 ) = 260;
indices( 16, 1 ) = 217;
indices( 16, 2 ) = 261;
indices( 24, 1 ) = 197;
indices( 24, 2 ) = 257;
indices( 90, 1 ) = 212;
indices( 91, 1 ) = 12;
indices( 94, 1 ) = 75;
indices( 100, 1 ) = 13;
indices( 132, 1 ) = 54;
indices( 132, 2 ) = 218;
indices( 134, 1 ) = 192;
indices( 135, 1 ) = 205;

% Two compounds still cannot be identified.

% Reorder Jang's results for the performances of classification models 
% (best performance, mean performance, variance in performance)

Jang_AUC = NaN( 265, 3 ); 

% Only consider the models fit to the sensitivity measure 'activity area', 
% indexed 1 to 8280.

temp = SangerperformanceTableAUCnew.DrugName;
temp = temp( 1:8280 ); 

temporder = unique( temp );

for k = 1:138
    
    drug = temporder( k );
    ind = strcmpi( temp, drug ); % Find the indices of all models for the drug in question
    
    temptemp = SangerperformanceTableAUCnew.AUCclassification( ind );
    
    
    if ~isnan( indices( k, 1 ) )
        
        if indices( k, 2 ) == 0
        
            Jang_AUC( indices( k, 1 ), 1 ) = nanmax( temptemp );
            Jang_AUC( indices( k, 1 ), 2 ) = nanmean( temptemp );
            Jang_AUC( indices( k, 1 ), 3 ) = nanstd( temptemp );
            
        else
            
            Jang_AUC( indices( k, 1 ), 1 ) = nanmax( temptemp );
            Jang_AUC( indices( k, 1 ), 2 ) = nanmean( temptemp );
            Jang_AUC( indices( k, 1 ), 3 ) = nanstd( temptemp );
            
            Jang_AUC( indices( k, 2 ), 1 ) = nanmax( temptemp );
            Jang_AUC( indices( k, 2 ), 2 ) = nanmean( temptemp );
            Jang_AUC( indices( k, 2 ), 3 ) = nanstd( temptemp );
            
        end
    
    end
    
end

% Compute concatenated AUCs for a fair comparison between studies

AUCs_conc = zeros( 265, 18 );

for drug = 1:265
    
    str = strcat( 'res_', num2str( drug ), '.mat' )
    
    load( str );
    
    % Get the binarised response
    
    temp = Response_AUC( :, drug );
    
    resp = quantile( temp( ~isnan( temp ) ), 0.25 );
    nonresp = quantile( temp( ~isnan( temp ) ), 0.75 );
    
    responseb = NaN( 968, 1 );
    
    responseb( temp <= resp, 1 ) = 1;
    responseb( temp >= nonresp, 1 ) = 0;
    
    % Concatenate prediction vector from test results

    predictions = zeros( 968, 18 );
    
    for j = 1:10
        
        temp = prediction_test( j, : );
        tempts = find( test( cv_part, j ) == 1 );
        
        for k = 1:18
            
            temptemp = double( temp{ 1, k + 19 } );
            predictions( tempts, k ) = temptemp;
            
        end
        
    end
    
    % Calculate concatenated AUCs
    
    AUC_conc = zeros( 18, 1 );
    
    for j = 1:18
        
        [ ~, ~, ~, concauc ] = perfcurve( responseb, predictions( :, j ), 1 );
        
        AUC_conc( j, 1 ) = concauc;
        
    end
    
    AUCs_conc( drug, : ) = AUC_conc;
    
end

% Figure 7.1

figure;

hold on;

symbols = { '^', 'v', '>', '<', 's', 'd', 'o', 'p', 'h', '+', 'x', '*' };

colours =  [ 0.4940, 0.1840, 0.5560;
             0.1254, 0.3372, 0.4941;
             0, 0.4470, 0.7410;
             0.3010, 0.7450, 0.9330;
             0.4588, 0.5725, 0.8352;
             0.8745, 0.7725, 0.8352;
             0.8901, 0.4901, 0.4941;
             0.4660, 0.6740, 0.1880;
             0.0980, 0.2901, 0.0431;
             0.8500, 0.3250, 0.0980;
             0.6350, 0.0780, 0.1840;
             0.9290, 0.6940, 0.1250 ];

for k = 1:12
    
    plot( results_mean_AUCs_test( :, 6 + k ), AUCs_conc( :, 6 + k ), symbols{ k }, ...
        'MarkerEdgeColor', colours( k, : ), 'MarkerFaceColor', colours( k, : ) )
    
end

plot( [ 0.3, 1 ],[ 0.3, 1 ], 'k--' )

hold off;

xlim( [ 0.3, 1 ] )
ylim( xlim )

grid on; box on;

ylabel( 'Concatenated ROC-AUCs' )
xlabel( 'Averaged ROC-AUCs' )

legend( Algorithms, 'FontSize', 13 ) 
legend( 'boxoff' )

set( gca, 'FontSize', 14 )

%  Pearsson correlation coefficient, R

temp1 = reshape( results_mean_AUCs_test( :, 7:18 ), 1, 12*265 );
temp2 = reshape( AUCs_conc( :, 7:18 ), 1, 12*265 );

corrcoef( temp1, temp2 )

% Figure 7.2

bestconcval = nanmax( AUCs_conc, [], 2 );

figure;

plot( [ 0, 0 ],[ 0, 15 ], 'k--' )

hold on

histogram( bestconcval - Jang_AUC( :, 1 ), 30, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', [0 0.4470 0.7410] )

hold off

grid on;
box on;

xlim( [ -0.4, 0.2 ] )
xticks( -0.4:0.1:0.2 ) 

xlabel( '\Delta_{perf} between two-step models and models by Jang et al.' )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure 7.3a

better = DrugOrder( bestconcval >= Jang_AUC( :, 1 ) , : );

index = find( bestconcval >= Jang_AUC( :, 1 ) );

bettermodels = cell( length( index ), 1 );
bettermodelsmat = zeros( length( index ), 18 );

for k = 1:length( index )
    
    bettermodels{ k, 1 } = find( AUCs_conc( index( k ), : ) >= Jang_AUC( index( k ), 1 ) );
    bettermodelsmat( k, bettermodels{ k, 1 } ) = 1;
    
end 

index2 = find( bestconcval + 0.2 <= Jang_AUC( :, 1 ) );

worstmodels = cell( length( index2 ), 1 );

for k = 1:length( index2 )
    
    worstmodels{ k, 1 } = find( AUCs_conc( index2( k ), : ) + 0.2 <= Jang_AUC( index2( k ), 1 ) );
    
end

figure;

b = bar( sum( bettermodelsmat ), 'FaceAlpha', 0.7 )


xtips1 = b( 1 ).XEndPoints;
ytips1 = b( 1 ).YEndPoints;
labels1 = string( b( 1 ).YData );
text( xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12 )

xlim( [ 0, 19 ] )
xticks( 1:18 )
xticklabels( [ 'Somatic mutation-based model', 'CNV-based model', 'Hypermethylation-based model', 'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model', Algorithms ] )
xtickangle( 30 )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure 7.3b

plotauc = NaN( 23, 18 );

for k = 1:23
    
    temp = bettermodels{ k, 1 };
    
    for l = 1:length( temp )
        
        plotauc( k, temp( l ) ) = AUCs_conc( index( k ), temp( l ) );
        
    end
    
end

better2 = DrugOrder( bestconcval >= Jang_AUC( :, 1 ), 1 );
better2{ 23, 1 } = 'Selumetinib';
better2{ 17, 1 } = 'Nutlin-3a';
better2{ 1, 1 } = better{ 1, 2 };
better2{ 12, 1 } = better{ 12, 2 };
better2{ 12, 1 } = 'Tanespimycin';

[ val, ind ] = sort( Jang_AUC( index, 1 ), 'ascend' );

symbols = { '.', '|', '_', '^', 'v', '>', '<', 's', 'd', 'o', 'p', 'h', '+', 'x', '*' };

colours2 = [ 1, 1, 1;
             1, 1, 1;
             1, 1, 1;
             0, 0, 0;
             0.1254, 0.3372, 0.4941;
             0.1098, 0.3176, 0.0666;
             0.7411, 0.5333, 0.8352;
             0.4588, 0.5725, 0.8352;
             0.1411, 0.5686, 0.8352;
             0.1764, 0.4431, 0.8588;
             0.4274, 0.0431, 0.1333;
             0.8156, 0.2901, 0.4156;
             0.8901, 0.4901, 0.4941;
             0.6039, 0.6823, 0.3764;
             0.7372, 0.8666, 0.8352;
             0.9254, 0.2901, 0.0549;
             0.8980, 0.4431, 0.1803;
             0.9568, 0.6823, 0.1411 ];

figure;
b = [];
count = 2;

b( 1 ) = stem( val, 'k', 'LineWidth', 2, 'MarkerSize', 7, ...
    'BaseValue', 0.5 ) % order them

hold on

for k = 1:23
    
    temp = plotauc( ind( k ), : );
    
    for l = 1:length( temp )
        
        if ~isnan( temp( l ) )
            
            if ind( k ) == 13
        
                b( count ) = plot( k + randn( 1, 1 )/6, temp( l ), symbols{ l - 3 }, 'MarkerSize', ...
                    7, 'MarkerEdgeColor', colours2( l, : ), 'MarkerFaceColor', colours2( l, : ) );
                count = count + 1;
                
            else
                
                plot( k + randn( 1, 1 )/6, temp( l ), symbols{ l - 3 }, 'MarkerSize', ...
                    7, 'MarkerEdgeColor', colours2( l, : ), 'MarkerFaceColor', colours2( l, : ) );
            
            end
        
        end
        
    end
    
end

ylim( [ 0.5, 1 ] )
ylabel( 'Concatenated test ROC-AUC' )
yticks( 0.5:0.05:1 )
yticklabels( { '0.5', '', '0.6', '', '0.7', '', '0.8', '', '0.9', '', '1' } )

xticks( 1:23 )
xticklabels( better2( ind ) )
xtickangle( 35 )

xlim( [ 0, 24 ] )

legend( b, [ 'Model by Jang et al.', 'Tissue descriptor-based Model', 'Pathway activation-based Model', 'Gene expression-based Model', Algorithms ], 'Location', 'Eastoutside' )

grid on;
box on;

set( gca, 'FontSize', 14 )

% Figure 7.4

Drug_redset = DrugOrder( ~isnan( Jang_AUC( :, 1 ) ), : );
drugtypes_redset = unique( cell2mat( Drug_redset( :, 5 ) ) ); % 19/21 found

pvals_enr_jang = zeros( 19, 1 );

for k = 1:19
    
    x = sum( cell2mat( better( :, 5 ) ) == drugtypes_redset( k ) );
    K = sum( cell2mat( Drug_redset( :, 5 ) ) == drugtypes_redset( k ) );
    
    pvals_enr_jang( k, 1 ) = hygecdf( x, 139, K, 23, 'upper' );
    
end

drugtypes( drugtypes_redset( [ 3, 7, 11 ] ) )

% Bonferroni correction with factor 19: only 3 remains significant
% Benjamini-Hochberg: same

% plot as bar plot

[ val, ind ] = sort( pvals_enr_jang, 'ascend' );

b1 = log10( val );
b1( 1:4 ) = 0;

b2 = log10( val );
b2( 1 ) = 0;
b2( 5:19 ) = 0;

b3 = log10( val );
b3( 2:19 ) = 0;

figure;

b = barh( [ b1, b2, b3 ], 'stacked' );

hold on;

plot( [ -1.3010, -1.3010 ], [ 0, 20 ], 'k--' )
plot( [-2.5798, -2.5798 ], [ 0, 20 ], 'k--' )

yticks( 1:19 )

yticklabels( drugtypes( drugtypes_redset( ind ) ) )
xlabel( 'log_{10}(p)' )

grid on;
box on;

legend( b, { 'Not significant', 'Significant in single testing', 'Significant in multiple testing' } )

set( gca, 'FontSize', 14 )

index = bestconcval >= Jang_AUC( :, 1 );

x = sum( index .* high' ); % 20 out of 99
hygecdf( x, 139, sum( high.* ~isnan( Jang_AUC( :, 1 ) )' ), 23, 'upper' ) % 0.0137

x = sum( index .* low' ); % 3 out of 40
hygecdf( x, 139, sum( low.* ~isnan( Jang_AUC( :, 1 ) )' ), 23, 'upper' ) % 0.9476

% Table 7.1

better_impscores = zeros( 1, 6 );
count = 1;

for k = 1:23
    
    temp = bettermodels{ k, 1 };
    
    for l = 1:length( temp )
        
        if temp( l ) > 7

            better_impscores( count, : ) = impscores_alg_two( temp( l ) - 7, index( k ), : );
            count = count + 1;
            
        else
            
        end
    
    end

end

foundbetter = sum( ~isnan( better_impscores ) )./180;
% 180 models over 23 drugs

worst_impscores = zeros( 1, 6 );
count = 1;

for k = 1:19
    
    temp = worstmodels{ k, 1 };
    
    for l = 1:length( temp )
        
        if temp( l ) > 7

            worst_impscores( count, : ) = impscores_alg_two( temp( l ) - 7, index2( k ), : );
            count = count + 1;
        
        end
    
    end

end

foundworst = sum( ~isnan( worst_impscores ) )./208;
% 209 models over 19 drugs

sum( ~isnan( better_impscores( :, 1:4 ) ) )
sum( ~isnan( worst_impscores( :, 1:4 ) ) )

[ ~,p, stats] = fishertest( [ 141, 39; 88, 120 ] ) % 4.6992e-13

% Table 7.2

[ ~,p, stats] = fishertest( [ 129, 51; 33, 175 ] ) % 5.7695e-30

% Table 7.3

[ ~,p, stats] = fishertest( [ 125, 55; 33, 175 ] ) % 1.0773e-27

% Table 7.4

[ ~,p, stats] = fishertest( [ 169, 11; 121, 87 ] ) % 4.2240e-17

% Figure 7.5

figure;

b1 = distributionPlot( better_impscores, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0, 'histOri', 'left', 'widthDiv', [2 1]  );
b2 = distributionPlot( worst_impscores, 'color', [0.6350, 0.0780, 0.1840], 'histOpt', 1, 'showMM', 0, 'histOri', 'right', 'widthDiv', [2 2] );

hold on

b3 = boxplot( better_impscores, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'b.', 'MedianStyle', 'target', 'Positions', 0.95:1:5.95 )
b4 = boxplot( worst_impscores, 'PlotStyle','traditional', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0.6350, 0.0780, 0.1840], 'Symbol', 'r.', 'MedianStyle', 'target', 'Positions', 1.05:1:6.05 )

% Add p-values for Welch's t-test

for k = 1:6
    
    [ ~, p ] = ttest2( better_impscores( :, k ), worst_impscores( :, k ), 'Vartype', 'unequal' )
    
end

pvals = { 'p = 9\cdot10^{-7}', 'p = 0.1380', 'p = 1\cdot10^{-4}', 'p = 3\cdot10^{-9}', 'p = 0.1502', 'p = 4\cdot10^{-29}' };

for k = 1:6
    
    plot( [ k - 0.15, k + 0.15 ], [ -0.35, -0.35 ], 'k' )
    plot( [ k - 0.15, k - 0.15 ], [ -0.35, -0.3], 'k' )
    plot( [ k + 0.15, k + 0.15 ], [ -0.35, -0.3], 'k' )
    
    text( k, -0.45, pvals{ k }, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )
    
end

xticks( 1:6 )
xticklabels( { 'Somatic mutation-based model', 'CNV-based model', 'Hypermethylation-based model', 'Tissue descriptor-based model', 'Pathway activation-based model', 'Gene expression-based model' } )
xtickangle( 25 )

ylabel( 'Normalised mean weights and importance scores' )

xlim( [ 0, 7 ] )
ylim( [ -0.5, 1.5 ] )

grid on;
box on;

lgd = legend( [ b1{ 1, 1 }( 1, 1 ), b2{ 1, 1 }( 1, 1 ) ], { 'Overperforming two-step models', 'Underperforming two-step models' }, 'Location', 'southeast', 'FontSize', 14, 'Orientation', 'horizontal' );
legend( 'boxoff' )

set( gca, 'FontSize', 14 )

% Table B.16 and Figure C.26

jangdiff = bestconcval - Jang_AUC( :, 1 );

tempcell = cell( 1, 21 );

for k = 1:21

    tempcell{ 1, k } = jangdiff( cell2mat( DrugOrder( :, 5 ) ) == k )';
    nanvar( tempcell{ 1, k } )
    sum( ~isnan( tempcell{ 1, k } ) )

end

nanvar( tempcell{ 1, 9 } )/nanvar( tempcell{ 1, 11 } )

figure;

t2 = tiledlayout( 7, 3 );

t2.TileSpacing = 'compact';

xlimits = { [-0.2, 0.05], [-0.3, 0.1], [-0.1, 0.1], [-0.3, 0.2], [-0.2, 0.05], [-0.3, -0.1], ...
    [-0.3, 0.2], [-0.4, 0], [-0.4, 0.2], [-0.3, 0.1], [-0.06, 0.02], [-0.15, 0.05], ...
    [-0.4, 0.1], [-0.15, 0.05], [0, 4], [0, 4], [-0.4, 0.2], [-0.5, 0.5], ...
    [-0.25, -0.05], [-0.4, 0.1], [-0.1, 0.1] };

for k = 1:21
   
    nexttile;
    
    histogram( tempcell{ 1, k }, 5 )
    xlabel( drugtypes{ k } )
    
    temp = ylim;
    
    if temp( end ) == 5
       
        yticks( 0:2.5:5 )
        yticklabels( { '0', '', '5' } )
        
    end
    
    grid on;
    box on;
    
    xlim( xlimits{ k } )
    
end

title( t2, '\bf{Distribution of \Delta_{pred} in predictive performance across distinct drug classes}' ) 
set( gca, 'FontSize', 11 )

% Figure C.27

figure;

t3 = tiledlayout( 1, 2 );

t3.TileSpacing = 'compact';

nexttile;

histogram( jangdiff( high ), 7 )
xlabel( 'Sample size \geq 600' )

set( gca, 'FontSize', 13 )

grid on;
box on;

xlim( [ -0.3, 0.2 ] )

nexttile;

histogram( jangdiff( low ), 8 )
xlabel( 'Sample size \leq 500' )

grid on;
box on;

xlim( [ -0.4, 0.1 ] )

title( t3, '\bf{Distribution of \Delta_{perf} between drugs with high and low sample sizes}', 'FontSize', 13 ) 

set( gca, 'FontSize', 13 )

sum( ~isnan( jangdiff( low ) ) )
sum( ~isnan( jangdiff( high ) ) )
nanvar( jangdiff( low ) )
nanvar( jangdiff( high ) )
nanmedian( jangdiff( low ) )
nanmedian( jangdiff( high ) )

nanvar( jangdiff( low ) ) \ nanvar( jangdiff( high ) )

% Table 7.5

index_red = ~isnan( jangdiff );
jangdiff_red = jangdiff( index_red );
drugs_red = cell2mat( DrugOrder( :, 5 ) );
drugs_red = drugs_red( index_red );

find( drugs_red == 18 ) % only one item   

jangdiff_red = jangdiff_red( [ 1:78, 80:139 ] );
drugs_red = drugs_red( [ 1:78, 80:139 ] );

temp = unique( drugs_red );
drugs_red2 = drugs_red;

for k = 1:length( temp )
    
    drugs_red2( drugs_red == temp( k ) ) = k;
    
end

[ p, F, df1, df2 ] = wanova ( jangdiff_red, drugs_red2 );

% p = 0.0501
% F = 2.1738
% df1 = 17
% df2 = 19.7353

samplesize = ones( size( jangdiff ) );
samplesize( high ) = 2;

jangdiff_red = jangdiff( ~isnan( jangdiff ) );
samplesize_red = samplesize( ~isnan( jangdiff ) );

[ p2, F2, df12, df22 ] = wanova ( jangdiff_red, samplesize_red );

% p2 = 3.1289e-07
% F2 = 32.4122
% df12 = 1
% df22 = 65.9061

% Figure 7.6

tempcell = cell( 1, 2 );

tempcell{ :, 1 } = jangdiff( low );
tempcell{ :, 2 } = jangdiff( high );


samplesize = ones( size( jangdiff ) );
samplesize( high ) = 2;

figure;

h = distributionPlot( tempcell, 'color', [0 0.4470 0.7410], 'histOpt', 1, 'showMM', 0 );

hold on
    
boxplot( jangdiff, samplesize, 'PlotStyle', 'compact', 'Widths', 0.05, 'BoxStyle', 'filled', 'Colors', [0 0.4470 0.7410], 'Symbol', 'r.', 'MedianStyle', 'target', 'Jitter', 0 )

hold off

ylim( [ - 0.4, 0.3] )
xlim( [ 0, 3 ] )

xticks( 1:2 ) 
yticks( -0.4:0.1:0.3 )

xticklabels( { 'Sample size \leq 500','Sample size \geq 600' } )
ylabel( '\Delta_{perf}' )

box on; grid on;
set( gca, 'FontSize', 14 )
