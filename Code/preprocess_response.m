function [ response_training, response_test, responseb_training, responseb_test ] = preprocess_response( response, trainingset, testset )

% Binarise the cellular drug sensitivity data, using the upper and lower
% quartile

threshold_responsive = 0.25;
threshold_resistant = 0.75;

response_training = response( trainingset );
response_test = response( testset );

% Compute the quantiles

resp = quantile( response_training, threshold_responsive );
nonresp = quantile( response_training, threshold_resistant );

% Binarise the data: 1 represents responders, 
%                    0 represents non-responders

responseb = NaN( length( response ), 1 );

responseb( response <= resp, 1 ) = 1; 
responseb( response >= nonresp, 1 ) = 0;

responseb_training = responseb( trainingset );
responseb_test = responseb( testset );


end