function [ model6_test, model6b_test, teststats6, auc ] = test_PCmodel( modelfun6, PC_test, cutoff6, response_test )


% Run the model on the test set

[ model6_test, ~ ] = predict( modelfun6, PC_test );
model6b_test = model6_test >= cutoff6;


% Evaluate the model performance

[ ~, ~, ~, auc ] = perfcurve( response_test, model6_test, 1 );

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model6b_test, response_test );

teststats6 = zeros( 5, 1 );

teststats6( 1, 1 ) = acc_test;
teststats6( 2, 1 ) = prec_test;
teststats6( 3, 1 ) = recall_test;
teststats6( 4, 1 ) = f1_test;
teststats6( 5, 1 ) = FDR_test;


end