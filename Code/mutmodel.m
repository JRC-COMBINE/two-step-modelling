function [ modelfun1, model1, model1b, clustmean1, cutoff1, modelstats1, auc1 ] = mutmodel( response_training, clustdata_training )

% Model the drug response based on mutational clustering

n = max( clustdata_training );
clustmean1 = zeros( n, 1 );


% Calculate the mean response per cluster

for i = 1:n
    
    clustmean1( i, 1 ) = nanmean( response_training( clustdata_training == i ) );
    
    if isnan( clustmean1( i, 1 ) )
       
        clustmean1( i, 1 ) = NaN ;
        
    end
    
end


% Write the model as a function to export with the cluster index as an input

model1 = clustmean1( clustdata_training, 1 );

modelfun1 = @(x) clustmean1( x, 1 );


% Binarize results

modelstats1 = zeros( 1, 101 );
index = isnan( model1 );

for i=1:101
    
    cutoff1 = ( i - 1 ) / 100;

    model1b = model1 >= cutoff1;
    model1b = +model1b;
    model1b( index ) = NaN;

    [ acc, ~, ~, ~, ~ ] = statev( model1b, response_training );
    
    modelstats1( 1, i ) = acc;
    
end


[ ~, i ] = max( modelstats1( 1, : ) );

cutoff1 = ( i - 1 )/100;

model1b = model1 >= cutoff1;
model1b = +model1b;
model1b( index ) = NaN;

[ ~, ~, ~, auc1 ] = perfcurve( response_training, model1, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model1b, response_training );

modelstats1 = zeros( 5, 1 );

modelstats1( 1, 1 ) = acc;
modelstats1( 2, 1 ) = prec;
modelstats1( 3, 1 ) = recall;
modelstats1( 4, 1 ) = f1;
modelstats1( 5, 1 ) = FDR;


end