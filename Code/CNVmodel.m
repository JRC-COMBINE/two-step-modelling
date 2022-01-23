function [ modelfun2, model2, model2b, clustmean2, cutoff2, modelstats2, auc2 ] = CNVmodel( response_training, clustdata_training )

% Model the drug response based on clustering on CNV patterns

n = max( clustdata_training );
clustmean2 = zeros( n, 1 );


% Calculate the mean response per cluster

for i = 1:n
    
    clustmean2(i,1) = nanmean( response_training( clustdata_training == i ) ); %#ok<*NANMEAN> 
    
    if isnan( clustmean2( i, 1 ) )
       
        clustmean2( i, 1 ) = NaN ;
        
    end
    
end


% Write the model as a function to export with the cluster index as an input

model2 = clustmean2( clustdata_training, 1 );

modelfun2 = @(x) clustmean2( x, 1 );


% Find the best classification cutoff out of 100 potential candidates 

n_cutoff = 100;

modelstats2 = zeros( 5, n_cutoff + 1 );
index = isnan( model2 );

for i = 1:n_cutoff + 1
    
    cutoff2 = ( i - 1 ) / n_cutoff;
    
    model2b = model2 >= cutoff2;
    model2b = +model2b;
    model2b( index ) = NaN;

    [ acc, prec, recall, f1, FDR ] = statev( model2b, response_training );
    
    modelstats2( 1, i ) = acc;
    modelstats2( 2, i ) = prec;
    modelstats2( 3, i ) = recall;
    modelstats2( 4, i ) = f1;
    modelstats2( 5, i ) = FDR;
    
end


[ ~, i ] = max( modelstats2( 1, : ) );

cutoff2 = ( i - 1 )/n_cutoff;

model2b = model2 >= cutoff2;
model2b = +model2b;
model2b( index ) = NaN;

[ ~, ~, ~, auc2 ] = perfcurve( response_training, model2, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model2b, response_training );

modelstats2 = zeros( 5, 1 );

modelstats2( 1, 1 ) = acc;
modelstats2( 2, 1 ) = prec;
modelstats2( 3, 1 ) = recall;
modelstats2( 4, 1 ) = f1;
modelstats2( 5, 1 ) = FDR;


end