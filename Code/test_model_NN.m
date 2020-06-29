function [ model_test, modelb_test, teststats, auc ] = test_model_NN( model, model_test, cutoff, responseb_test, normfunction )

% Run the model on the test set

model_test = model( model_test' ); 
model_test = normfunction( model_test );


modelb_test = model_test >= cutoff;

index = isnan( model_test );

modelb_test = +modelb_test;
modelb_test( index ) = NaN;


% Calculate the performance metrics

[ ~, ~, ~, auc ] = perfcurve( responseb_test, model_test, 1 );

[ acc_test, prec_test, recall_test, f1_test, FDR_test ] = statev( modelb_test', responseb_test );

teststats = zeros( 5, 1 );

teststats( 1, 1 ) = acc_test;
teststats( 2, 1 ) = prec_test;
teststats( 3, 1 ) = recall_test;
teststats( 4, 1 ) = f1_test;
teststats( 5, 1 ) = FDR_test;


end