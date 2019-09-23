function [ modelfun5, model5, model5b, cutoff5, modelstats5, auc5 ] = pathwaymodel( response, pathways )

% Model the drug response based on pathway-activation regression

% Calculate the regression coefficients

modelfun5 = fitlm( pathways, response, 'linear', 'RobustOpts', 'on', 'Intercept', false );
model5 = modelfun5.Fitted;

% Find the optimal cutoff

modelstats5 = zeros( 5, 101 );
index = isnan( model5 );

for i = 1:101
    
    cutoff5 = ( i - 1 )/100;
    
    model5b = model5 >= cutoff5;
    model5b = +model5b;
    model5b( index ) = NaN;
    
    [ acc, prec, recall, f1, FDR ] = statev( model5b, response );
    
    modelstats5( 1, i ) = acc;
    modelstats5( 2, i ) = prec;
    modelstats5( 3, i ) = recall;
    modelstats5( 4, i ) = f1;
    modelstats5( 5, i ) = FDR;
    
end

[ ~, i ] = max( modelstats5( 1, : ) );

cutoff5 = ( i - 1 )/100;

model5b = model5 >= cutoff5;
model5b = +model5b;
model5b( index ) = NaN;

% Evaluate the performance

[ ~, ~, ~, auc5 ] = perfcurve( response, model5, 1 );

[ acc, prec, recall, f1, FDR ] = statev( model5b, response );

modelstats5 = zeros( 5, 1 );

modelstats5( 1, 1 ) = acc;
modelstats5( 2, 1 ) = prec;
modelstats5( 3, 1 ) = recall;
modelstats5( 4, 1 ) = f1;
modelstats5( 5, 1 ) = FDR;


end