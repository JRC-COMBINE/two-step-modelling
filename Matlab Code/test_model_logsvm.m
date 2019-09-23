function [ model_test, teststats ] = test_model_logsvm( model, model_test, responseb_test )

% Run the model on the test set

model_test = predict( model, model_test );


% Calculate the performance metrics

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( model_test, responseb_test );

teststats = zeros( 5, 1 );

teststats( 1, 1 ) = acc_test;
teststats( 2, 1 ) = prec_test;
teststats( 3, 1 ) = recall_test;
teststats( 4, 1 ) = f1_test;
teststats( 5, 1 ) = FDR_test;


end