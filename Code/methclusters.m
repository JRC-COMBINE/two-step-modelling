function [ clusters3, set3, set3_large, clusterinfo3 ] = methclusters( methdata_training, methdata, response )


% Remove all genes with not enough cell lines in any of the categories

index = zeros( size( methdata_training, 2 ), 1 );

for j = 1:size( methdata_training, 2 )
    
    hypermethylations = sum( methdata_training( :, j ) == 1 ) >= 15;
    normal = sum( methdata_training( :, j ) == 0 ) >= 15;
    
    index( j, 1 ) = normal * hypermethylations;
    
end

index = find( index == 1 );


% Check for single hypermethylations that change responsiveness 

relevantmethylations = nan( 1, length( index ) );

for j = 1:length( index )
    
    meth = response( methdata_training( :, index( j ) ) == 1 );
    normal = response( methdata_training( :, index( j ) ) == 0 );
    
    [ ~, p ] = ttest2( meth, normal );
    relevantmethylations( 1, j ) = p;
    
end


% Correct for multiple testing

testcorrection = sum( ~isnan( relevantmethylations ) );

relevantmethylations = relevantmethylations .* testcorrection;


% Select statistically significant methylations

[ clusters3, set3, set3_large ] = clusteringmethylations( methdata, relevantmethylations, index );


% Associate clusters with state of significant methylations

clusterinfo3 = NaN( max( clusters3 ), length( set3 ) );

for j = 1:max( clusters3 )
    
    ind = clusters3 == j;
    
    for k = 1:length( set3 )
        
        clusterinfo3( j, k ) = unique( methdata( ind, set3( k ) ) );
        
    end
    
end


end