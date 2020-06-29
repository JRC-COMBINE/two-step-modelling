function [ model2_test, model2b_test, teststats2, auc ] = test_CNVmodel( modelfun2, clustdata_test, cutoff2, response_test, n2 )

% Check if some test cell lines are from clusters that are not present in
% the training set of clusters

present = find( clustdata_test <= n2 ) ;


% Run the model on the test set

model2_test = NaN( length( clustdata_test ), 1 );

model2_test( present, 1 ) = modelfun2( clustdata_test( present ) );

index = isnan( model2_test );

model2b_test = model2_test >= cutoff2;
model2b_test = +model2b_test;
model2b_test( index ) = NaN;


% Evaluate the model performance

if sum( model2b_test == 1 ) > 1 && sum( model2b_test == 0 ) > 1
    
[ ~, ~, ~, auc ] = perfcurve( response_test, model2_test, 1 );

else
    
    auc = NaN;
    
end

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model2b_test, response_test );

teststats2 = zeros( 5, 1 );

teststats2( 1, 1 ) = acc_test;
teststats2( 2, 1 ) = prec_test;
teststats2( 3, 1 ) = recall_test;
teststats2( 4, 1 ) = f1_test;
teststats2( 5, 1 ) = FDR_test;


end