function [ modelfun6, model6, model6b, cutoff6, modelstats6, auc6 ] = PCmodel( response, PC )

% Model the drug response based on PCA scores

% Calculate the regression coefficients

modelfun6 = fitlm( PC, response, 'linear', 'RobustOpts', 'on', 'Intercept', false );
model6 = modelfun6.Fitted;

% Find the best classification cutoff out of 100 potential candidates 
n_cutoff = 100;

modelstats6 = zeros( 5, n_cutoff + 1 );
index = isnan( model6 );

for i = 1:n_cutoff + 1
    
    cutoff6 = ( i - 1 )/n_cutoff;
    
    model6b = model6 >= cutoff6;
    model6b = +model6b;
    model6b( index ) = NaN;
    
    [ acc, prec, recall, f1, FDR ] = statev( model6b, response );
    
    modelstats6( 1, i ) = acc;
    modelstats6( 2, i ) = prec;
    modelstats6( 3, i ) = recall;
    modelstats6( 4, i ) = f1;
    modelstats6( 5, i ) = FDR;
    
end

[ ~, i ] = max( modelstats6( 1, : ) );

cutoff6 = ( i - 1 )/n_cutoff;

model6b = model6 >= cutoff6;
model6b = +model6b;
model6b( index ) = NaN;

% Evaluate the performance

[ ~, ~, ~, auc6 ] = perfcurve( response, model6, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model6b, response );

modelstats6 = zeros( 5, 1 );

modelstats6( 1, 1 ) = acc;
modelstats6( 2, 1 ) = prec;
modelstats6( 3, 1 ) = recall;
modelstats6( 4, 1 ) = f1;
modelstats6( 5, 1 ) = FDR;


end