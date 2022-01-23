function [ modelfun3, model3, model3b, clustmean3, cutoff3, modelstats3, auc3 ] = methmodel( response_training, clustdata_training )

% Model the drug response based on clustering on hypermethylation patterns

n = max( clustdata_training );
clustmean3 = zeros( n, 1 );


% Calculate the mean response per cluster

for i = 1:n
    
    clustmean3( i, 1 ) = nanmean( response_training( clustdata_training == i ) ); %#ok<*NANMEAN> 
    
    if isnan( clustmean3( i, 1 ) )
       
        clustmean3( i, 1 ) = NaN ;
        
    end
    
end


% Write the model as a function to export with the cluster index as an input

model3 = clustmean3( clustdata_training, 1 );

modelfun3 = @(x) clustmean3( x, 1 );


% Find the best classification cutoff out of 100 potential candidates 

n_cutoff = 100;

modelstats3 = zeros( 5, n_cutoff + 1 );
index = isnan( model3 ); 

for i = 1:n_cutoff + 1
    
    cutoff3 = ( i - 1 ) / n_cutoff;
    
    model3b = model3 >= cutoff3;
    model3b = +model3b;
    model3b( index ) = NaN;

    [ acc, prec, recall, f1, FDR ] = statev( model3b, response_training );
    
    modelstats3( 1, i ) = acc;
    modelstats3( 2, i ) = prec;
    modelstats3( 3, i ) = recall;
    modelstats3( 4, i ) = f1;
    modelstats3( 5, i ) = FDR;
    
end


[ ~, i ] = max( modelstats3( 1, : ) );

cutoff3 = ( i - 1 )/n_cutoff;

model3b = model3 >= cutoff3;
model3b = +model3b;
model3b( index ) = NaN;

[ ~, ~, ~, auc3 ] = perfcurve( response_training, model3, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model3b, response_training );

modelstats3 = zeros( 5, 1 );

modelstats3( 1, 1 ) = acc;
modelstats3( 2, 1 ) = prec;
modelstats3( 3, 1 ) = recall;
modelstats3( 4, 1 ) = f1;
modelstats3( 5, 1 ) = FDR;


end