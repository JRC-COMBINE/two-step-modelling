function [ set1, set1_large ] = methclusters( methdata_training, methdata, response )


%% Remove all genes with not enough cell lines in any of the categories

index = zeros( size( methdata_training, 2 ), 1 );

for j = 1:size( methdata_training, 2 )
    
    hypermethylations = sum( methdata_training( :, j ) == 1 ) >= 15;
    normal = sum( methdata_training( :, j ) == 0 ) >= 15;
    
    index( j, 1 ) = normal * hypermethylations;
    
end

index = find( index == 1 );


%% Check for single hypermethylations that change responsiveness 

relevantmethylations = nan( 1, length( index ) );

for j = 1:length( index )
    
    meth = response( methdata_training( :, index( j ) ) == 1 );
    normal = response( methdata_training( :, index( j ) ) == 0 );
    
    [ ~, p ] = ttest2( meth, normal );
    relevantmethylations( 1, j ) = p;
    
end


%% Correct for multiple testing

testcorrection = sum( ~isnan( relevantmethylations ) );

relevantmethylations = relevantmethylations .* testcorrection;


%% Select statistically significant features

set1 = find( relevantmethylations <= 0.05 );

set1_large = index( set1 );


%% Remove highly-correlated features that are present/absent in the exact same set of cell lines

set2 = set1;

count = 2;

while count <= length( set2 )
    
    temp = abs( methdata( :, index( set2( count ) ) ) );
    
    for j = 1:count-1
        
        temp2 = abs( methdata( :, index( set2( j ) ) ) );
        
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
    
    [ ~, ind ] = sort( relevantmethylations, 'ascend' );
    set1 = ind( 1:15 );
    
    % Remove redundant features
    
    set2 = set1;
    
    count = 2;
    
    while count <= length( set2 )
        
        temp = abs( methdata( :, index( set2( count ) ) ) );
        
        for j = 1:count-1
            
            temp2 = abs( methdata( :, index( set2( j ) ) ) );
            
            if isequal( temp, temp2 )
                
                set2 = set2( set2 ~= set2( count ) );
                count = count - 1;
                
                break;
                
            end
            
        end
        
        count = count + 1;
        
    end
    
    set1 = index( set2 );
    
else
    
    set1 = index( set1 );
    
end


end