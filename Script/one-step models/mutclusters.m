function [ set1, set_large ] = mutclusters( mutdata_training, mutdata, response )


%% Check for single mutations that change responsiveness 

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


%% Correct for multiple testing

testcorrection = sum( ~isnan( relevantmutations ) );

relevantmutations = relevantmutations .* testcorrection;


%% Select statistically significant features

set1 = find( relevantmutations <= 0.05 );

set_large = set1;


%% Remove highly-correlated features that are present/absent in the exact same set of cell lines

set2 = set1;

count = 2;

while count <= length( set2 )
    
    temp = abs( mutdata( :, set2( count ) ) );
    
    for j = 1:count-1
        
        temp2 = abs( mutdata( :, set2( j ) ) );
        
        if isequal( temp, temp2 )
            
            set2 = set2( set2 ~= set2( count ) );
            count = count - 1;
            
            break;
            
        end
        
    end
    
    count = count + 1;
    
end

set1 = set2;


%% Reduce the set of relevant features to guarantee a reasonable run time

if length( set1 ) > 15
    
    % Define a new set by ordering the features with respect to relevance
    
    [ ~, ind ] = sort( relevantmutations, 'ascend' );
    set1 = ind( 1:15 );
    
    % Remove redundant features
    
    set2 = set1;
    
    count = 2;
    
    while count <= length( set2 )
        
        temp = abs( mutdata( :, set2( count ) ) );
        
        for j = 1:count-1
            
            temp2 = abs( mutdata( :, set2( j ) ) );
            
            if isequal( temp, temp2 )
                
                set2 = set2( set2 ~= set2( count ) );
                count = count - 1;
                
                break;
                
            end
            
        end
        
        count = count + 1;
        
    end
    
    set1 = set2;
    
end

end