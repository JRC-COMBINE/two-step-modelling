function [ clusters1, set1, set1_large, clusterinfo1 ] = mutclusters( mutdata_training, mutdata, response )


% Check for single mutations that change responsiveness 

[ ~, n2 ] = size( mutdata ); 

relevantmutations = nan( 1, n2 );

for i = 1:n2
    
    mut = response( mutdata_training( :, i ) == 1 );
    nonmut = response( mutdata_training( :, i ) == 0 );
    
    if min( length( mut ),length( nonmut ) ) >= 15
        
        [ ~, p ] = ttest2( mut, nonmut );
        relevantmutations( 1, i ) = p;
        
    end
    
end


% Correct for multiple testing

testcorrection = sum( ~isnan( relevantmutations ) );

relevantmutations = relevantmutations .* testcorrection;


% Select statistically significant mutations

[ clusters1, set1, set1_large ] = clusteringmutations( mutdata, relevantmutations );


% Associate clusters with state of significant mutations

clusterinfo1 = NaN( max( clusters1 ), length( set1 ) );

for j = 1:max( clusters1 )
    
    ind = clusters1 == j;
    
    for k = 1:length( set1 )
        
        clusterinfo1( j, k ) = unique( mutdata( ind, set1( k ) ) );
        
    end
    
end


end