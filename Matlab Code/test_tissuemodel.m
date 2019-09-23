function [ model4_test, model4b_test, teststats4, auc ] = test_tissuemodel( modelfun4, clustdata_test, cutoff4, response_test, n4 )

% Check if some test cell lines are from clusters that are not present in
% the training set of clusters

present = find( clustdata_test <= n4 ) ;


% Run the model on the test set

model4_test = NaN( length( clustdata_test ), 1 );

model4_test( present, 1 ) = modelfun4( clustdata_test( present ) );

index = isnan( model4_test );

model4b_test = model4_test >= cutoff4;
model4b_test = +model4b_test;
model4b_test( index ) = NaN;


% Evaluate the model performance

if sum(  model4b_test == 1 ) > 1 && sum( model4b_test == 0 ) > 1 
    
[ ~, ~, ~, auc ] = perfcurve( response_test, model4_test, 1 );

else
    
    auc = NaN;
    
end

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model4b_test, response_test );

teststats4 = zeros( 5, 1 );

teststats4( 1, 1 ) = acc_test;
teststats4( 2, 1 ) = prec_test;
teststats4( 3, 1 ) = recall_test;
teststats4( 4, 1 ) = f1_test;
teststats4( 5, 1 ) = FDR_test;


end