function [ clusters4, set4 ] = tissueclusters( tissuedata_training, tissuedata, response )


% Remove all tissues with less than ten cell lines 

index = find( sum( tissuedata_training, 1 ) >= 10 );


% Check for tissues that change responsiveness  

relevanttissues = nan( 1, length( index ) );

for j = 1:length( index )
    
    tissue = response( tissuedata_training( :, index( j ) ) == 1 );
    rest = response( tissuedata_training( :, index( j ) ) == 0 );
    
    [ ~, p ] = ttest2( rest, tissue );
    relevanttissues( 1, j ) = p;
    
end


% Correct for multiple testing

testcorrection = sum( ~isnan( relevanttissues ) );

relevanttissues = relevanttissues .* testcorrection;


% Select statistically significant tissues

set4 = index(  relevanttissues <= 0.05  );


% Save the resulting clusters of tissues

clusters4 = ones( 1, size( tissuedata, 1 ) );
count = 2;


for j = 1:length( set4 )
    
    clusters4( tissuedata( :, set4( j ) ) == 1 ) = count;
    count = count + 1;
    
end


end