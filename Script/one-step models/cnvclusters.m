function [ set1, set1_large ] = cnvclusters( cnvdata_training, cnvdata, response )


%% Remove all genes with not enough cell lines for at least two of the categories

index = zeros( size( cnvdata_training, 2 ), 1 );

for j = 1:size( cnvdata_training, 2 )
    
    temp_del = sum( cnvdata_training( :, j ) == -1 ) >= 15;
    temp_norm = sum( cnvdata_training( :, j ) == 0 ) >= 15;
    temp_amp = sum( cnvdata_training( :, j ) == 1 ) >= 15;
    
    index( j, 1 ) = temp_del * temp_norm || temp_norm * temp_amp || temp_del * temp_amp; %#ok<BDLOG>
    
end

index = find( index == 1 );


%% Check for single CNVs that change responsiveness  

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


%% Correct for multiple testing

testcorrection = sum( ~isnan( relevantcnvs ) );

relevantcnvs = relevantcnvs .* testcorrection;


%% Select statistically significant features

set1 = find( relevantcnvs <= 0.05 );

set1_large = index( set1 );


%% Remove highly-correlated features that are present/absent in the exact same set of cell lines

set2 = set1;

count = 2;

while count <= length( set2 )
    
    temp = abs( cnvdata( :, index( set2( count ) ) ) );
    
    for j = 1:count-1
        
        temp2 = abs( cnvdata( :, index( set2( j ) ) ) );
        
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
    
    [ ~, ind ] = sort( relevantcnvs, 'ascend' );
    set1 = ind( 1:15 );
    
    % Remove redundant features
    
    set2 = set1;
    
    count = 2;
    
    while count <= length( set2 )
        
        temp = abs( cnvdata( :, set2( count ) ) );
        
        for j = 1:count-1
            
            temp2 = abs( cnvdata( :, set2( j ) ) );
            
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