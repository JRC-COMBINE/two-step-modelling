function [ model3_test, model3b_test, teststats3, auc ] = test_methmodel( modelfun3, clustdata_test, cutoff3, response_test, n3 )

% Check if some test cell lines are from clusters that are not present in
% the training set of clusters

present = find( clustdata_test <= n3 ) ;


% Run the model on the test set

model3_test = NaN( length( clustdata_test ), 1 );

model3_test( present, 1 ) = modelfun3( clustdata_test( present ) );

index = isnan( model3_test );

model3b_test = model3_test >= cutoff3;
model3b_test = +model3b_test;
model3b_test( index ) = NaN;


% Evaluate the model performance

if sum(  model3b_test == 1 ) > 1 && sum( model3b_test == 0 ) > 1 
    
[ ~, ~, ~, auc ] = perfcurve( response_test, model3_test, 1 );

else
    
    auc = NaN;
    
end

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model3b_test, response_test );

teststats3 = zeros( 5, 1 );

teststats3( 1, 1 ) = acc_test;
teststats3( 2, 1 ) = prec_test;
teststats3( 3, 1 ) = recall_test;
teststats3( 4, 1 ) = f1_test;
teststats3( 5, 1 ) = FDR_test;


end