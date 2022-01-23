function [ h, h2, meancoeffs_n ] = heatmapplot( drug, coeffs, DrugOrder )


meancoeffs = zeros( 6, 11 ); % no NN, no Bayes

for k = 2:12
    
    temp = coeffs( :, k );
    tempmat = NaN( 10, 6 );
    
    for l = 1:10
        
        temptemp = temp{ l, 1 };
        
        if sum( coeffs{ l, 1 } ) == length( temp{ l, 1 } ) % no constant
        
            tempmat( l, coeffs{ l, 1 } == 1 ) = temptemp;
            
        else % constant needs to be excluded
            
            tempmat( l, coeffs{ l, 1 } == 1 ) = temptemp( 2:length( temptemp ) );
        
        end
        
    end
    
    meancoeffs( :, k-1 ) = nanmean( abs( tempmat ) );
    
end

% calculate index

index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;

xvalues = { 'Somatic mutation', 'CNV', 'Hypermethylation', 'Tissue descriptors', 'Pathway activation', 'Gene expression' };
yvalues = { 'LASSO-regularised lin.regr.', 'E.net-regularised lin.regr.', 'Ridge-regularised lin.regr.','LASSO-regularised log.regr.', ...
'E.net-regularised log.regr.', 'Ridge-regularised log.regr.', 'Ridge-regularised SVM', 'LASSO-regularised SVM', 'Bagged DT ensemble', 'Boosted DT ensemble', ...
'Random forest' };

% normalise

meancoeffs_n = meancoeffs./nanmax( meancoeffs );

% fix number of decimal digits (4)

meancoeffs_str = string( meancoeffs_n );

for k = 1:numel( meancoeffs_str )
    
    meancoeffs_str( k ) = sprintf( '%.4f', meancoeffs_str( k ) );
    meancoeffs_n( k ) = str2num( meancoeffs_str( k ) );
    
end



% heatmap

figure;

h = heatmap( xvalues( index ), yvalues, meancoeffs_n( index, : )', 'Colormap', jet );

h.Title = [ 'Relative importance of data types in second-step models for the drug ', ' ', DrugOrder{ drug, 1 } ];
h.XLabel = 'Data types';
h.YLabel = 'Algorithms';

% clustergram

h2 = clustergram( meancoeffs_n( index, : )', 'Symmetric', 'false', 'Colormap', 'jet', 'ColumnLabels', xvalues( index ), ...
    'RowLabels', yvalues, 'ColumnLabelsRotate', 25, 'RowLabelsRotate', 25, 'Annotate', true, 'AnnotPrecision', 4, 'DisplayRatio', [ 0.2, 0.2 ] );

addTitle( h2, [ 'Relative importance of data types in second-step models for the drug', ' ', DrugOrder{ drug, 1 } ] );


end