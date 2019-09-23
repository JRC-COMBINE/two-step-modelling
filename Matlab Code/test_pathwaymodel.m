function [ model5_test, model5b_test, teststats5, auc ] = test_pathwaymodel( modelfun5, pathways_test, cutoff5, response_test )


% Run the model on the test set

[ model5_test, ~ ] = predict( modelfun5, pathways_test );
model5b_test = model5_test >= cutoff5;


% Evaluate the model performance

[ ~, ~, ~, auc ] = perfcurve( response_test, model5_test, 1 );

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model5b_test, response_test );

teststats5 = zeros( 5, 1 );

teststats5( 1, 1 ) = acc_test;
teststats5( 2, 1 ) = prec_test;
teststats5( 3, 1 ) = recall_test;
teststats5( 4, 1 ) = f1_test;
teststats5( 5, 1 ) = FDR_test;


end