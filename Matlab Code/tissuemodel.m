function [ modelfun4, model4, model4b, cutoff4, modelstats4, auc4 ] = tissuemodel( response_training, clustdata_training )

% Model the drug response based on clustering on tissue types

n = max( clustdata_training );
clustmean = zeros( n, 1 );


% Calculate the mean response per cluster

for i = 1:n
    
    clustmean(i,1) = nanmean( response_training( clustdata_training == i ) );
    
    if isnan( clustmean( i, 1 ) )
       
        clustmean( i, 1 ) = NaN ;
        
    end
    
end


% Write the model as a function to export with the cluster index as an input

model4 = clustmean( clustdata_training, 1 );

modelfun4 = @(x) clustmean( x, 1 );


% Binarize results

modelstats4 = zeros( 5, 101 );
index = isnan( model4 );

for i=1:101
    
    cutoff4 = ( i - 1 ) / 100;
    
    model4b = model4 >= cutoff4;
    model4b = +model4b;
    model4b( index ) = NaN;

    [ acc, prec, recall, f1, FDR ] = statev( model4b, response_training );
    
    modelstats4( 1, i ) = acc;
    modelstats4( 2, i ) = prec;
    modelstats4( 3, i ) = recall;
    modelstats4( 4, i ) = f1;
    modelstats4( 5, i ) = FDR;
    
end


[ ~, i ] = max( modelstats4( 1, : ) );

cutoff4 = ( i - 1 )/100;

model4b = model4 >= cutoff4;
model4b = +model4b;
model4b( index ) = NaN;

[ ~, ~, ~, auc4 ] = perfcurve( response_training, model4, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model4b, response_training );

modelstats4 = zeros( 5, 1 );

modelstats4( 1, 1 ) = acc;
modelstats4( 2, 1 ) = prec;
modelstats4( 3, 1 ) = recall;
modelstats4( 4, 1 ) = f1;
modelstats4( 5, 1 ) = FDR;


end