function set4 = tissueclusters( tissuedata_training, response )


%% Remove all tissues with less than ten cell lines 

index = find( sum( tissuedata_training, 1 ) >= 10 );


%% Check for tissues that change responsiveness  

relevanttissues = nan( 1, length( index ) );

for j = 1:length( index )
    
    tissue = response( tissuedata_training( :, index( j ) ) == 1 );
    rest = response( tissuedata_training( :, index( j ) ) == 0 );
    
    [ ~, p ] = ttest2( rest, tissue );
    relevanttissues( 1, j ) = p;
    
end


%% Correct for multiple testing

testcorrection = sum( ~isnan( relevanttissues ) );

relevanttissues = relevanttissues .* testcorrection;


%% Select statistically significant tissues

set4 = index(  relevanttissues <= 0.05  );


end