function [ clusters2, set2, set2_large, clusterinfo2 ] = cnvclusters( cnvdata_training, cnvdata, response )


% Remove all CNVs with not enough cell lines for at least two of the categories (deletion, regular, amplification)
min_num = 15;

index = zeros( size( cnvdata_training, 2 ), 1 );

for j = 1:size( cnvdata_training, 2 )
    
    temp_del = sum( cnvdata_training( :, j ) == -1 ) >= min_num;
    temp_norm = sum( cnvdata_training( :, j ) == 0 ) >= min_num;
    temp_amp = sum( cnvdata_training( :, j ) == 1 ) >= min_num;
    
    index( j, 1 ) = temp_del * temp_norm || temp_norm * temp_amp || temp_del * temp_amp; %#ok<BDLOG>
    
end

index = find( index == 1 );


% Check for single CNVs that change responsiveness  

relevantcnvs = nan( 1, length( index ) );

for j = 1:length( index )
    
    deletion = response( cnvdata_training( :, index( j ) ) == -1 );
    normal = response( cnvdata_training( :, index( j ) ) == 0 );
    amplification = response( cnvdata_training( :, index( j ) ) == 1 );
    

    if isempty( deletion )
       
        [ ~, p ] = ttest2( normal, amplification );
        relevantcnvs( 1, j ) = p;
        
    end
    
    if isempty( amplification )
        
        [ ~, p ] = ttest2( normal, deletion );
        relevantcnvs( 1, j ) = p;
        
    end
    
end


% Correct for multiple testing

testcorrection = sum( ~isnan( relevantcnvs ) );

relevantcnvs = relevantcnvs .* testcorrection;


% Select statistically significant CNVs

[ clusters2, set2, set2_large ] = clusteringcnvs( cnvdata, relevantcnvs, index );


% Associate clusters with state of significant CNVs

clusterinfo2 = NaN( max( clusters2 ), length( set2 ) );

for j = 1:max( clusters2 )
    
    ind = clusters2 == j;
    
    for k = 1:length( set2 )
        
        clusterinfo2( j, k ) = unique( cnvdata( ind, set2( k ) ) );
        
    end
    
end


end