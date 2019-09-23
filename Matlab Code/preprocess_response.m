function [ response_training, response_test, responseb_training, responseb_test ] = preprocess_response( response, trainingset, testset )


response_training = response( trainingset );
response_test = response( testset );


resp = quantile( response_training, 0.25 );
nonresp = quantile( response_training, 0.75 );


responseb = NaN( length( response ), 1 );

responseb( response <= resp, 1 ) = 1; 
responseb( response >= nonresp, 1 ) = 0;

responseb_training = responseb( trainingset );
responseb_test = responseb( testset );


end