function [ model1_test, model1b_test, teststats1, auc ] = test_mutmodel( modelfun1, clustdata_test, cutoff1, response_test, n1 )

% Check if some test cell lines are from clusters that are not present in
% the training set of clusters

present = find( clustdata_test <= n1 ) ;


% Run the model on the test set

model1_test = NaN( length( clustdata_test ), 1 );

model1_test( present, 1 ) = modelfun1( clustdata_test( present ) );

index = isnan( model1_test );

model1b_test = model1_test >= cutoff1;
model1b_test = +model1b_test;
model1b_test( index ) = NaN;


% Evaluate the model performance

if sum(  model1b_test == 1 ) > 1 && sum( model1b_test == 0 ) > 1 
    
[ ~, ~, ~, auc ] = perfcurve( response_test, model1_test, 1 );

else
    
    auc = NaN;
    
end

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model1b_test, response_test );

teststats1 = zeros( 5, 1 );

teststats1( 1, 1 ) = acc_test;
teststats1( 2, 1 ) = prec_test;
teststats1( 3, 1 ) = recall_test;
teststats1( 4, 1 ) = f1_test;
teststats1( 5, 1 ) = FDR_test;


end