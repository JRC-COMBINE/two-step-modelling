function [ clusters, set, set_large ] = clusteringmethylations( data, relevantvariables, index1 )


%% Select statistically significant features

% Use a level of significance of 0.05
alpha = 0.05;
% Use no more than 15 significant methylation features for clustering
max_num_features = 25;
% Compute no more than 32 non-empty methylation-based clusters of cell lines
max_num_clusters = 32;

set1 = find( relevantvariables <= alpha );
n1 = size( data, 1 );

set_large = index1( set1 );


%% Remove highly-correlated features that are present/absent in the exact same set of cell lines

set2 = set1;

count = 2;

while count <= length( set2 )
    
    temp = abs( data( :, index1( set2( count ) ) ) );
    
    for j = 1:count-1
        
        temp2 = abs( data( :, index1( set2( j ) ) ) );
        
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

if length( set1 ) > max_num_features
    
    % Define a new set by ordering the features with respect to relevance
    
    [ ~, ind ] = sort( relevantvariables, 'ascend' );
    set1 = ind( 1:max_num_features );
    
    % Remove redundant features
    
    set2 = set1;
    
    count = 2;
    
    while count <= length( set2 )
        
        temp = abs( data( :, index1( set2( count ) ) ) );
        
        for j = 1:count-1
            
            temp2 = abs( data( :, index1( set2( j ) ) ) );
            
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


%% Compute the initial number of clusters 

clusters = zeros( 1, n1 );
run = 0;

for j = 1:length( set1 )
   
    var = index1( set1( j ) );
    clusters( 1, abs( data( :, var ) ) == 1 ) = clusters( 1, abs( data( :, var ) ) == 1 ) + 2^run; 
    
    run = run + 1;
 
end

clusters = clusters + 1;


%% Remove empty clusters

index = 1;
m = max( clusters );

while index <= m
    
    temp = sum( clusters == index );
    
    if temp == 0
        
        temp = find( clusters > index );
        clusters( 1, temp ) = clusters( 1, temp ) - 1;
        
        m = m - 1;
        
    else
        
        index = index + 1;
        
    end
    
end


%% Adjust the size of the set of relevant features if it yields too many non-empty clusters

mcl = max( clusters );

if mcl > max_num_clusters
    
    % Repeat the clustering procedure and remove features in a stepwise manner, starting with the least significant ones
    
    temp1 = min( sum( relevantvariables <= alpha ), max_num_features );
    [ ~, ind ] = sort( relevantvariables, 'ascend' );
    
    set1 = ind( 1:temp1 );
    
    % Remove redundant features

    set2 = set1;
    
    count = 2;
    
    while count <= length( set2 )
        
        temp = abs( data( :, index1( set2( count ) ) ) );
        
        for j = 1:count-1
            
            temp2 = abs( data( :, index1( set2( j ) ) ) );
            
            if isequal( temp, temp2 )
                
                set2 = set2( set2 ~= set2( count ) );
                count = count - 1;
                
                break;
                
            end
            
        end
        
        count = count + 1;
        
    end
    
    % Continue with a set of non-redundant significant features, sorted with respect to significance
    
    set1 = set2; 
    temp1 = length( set1 );
    
    while max( clusters ) > max_num_clusters
        
        temp1 = temp1 - 1;
        set1 = set1( 1:temp1 );
        
        % Repeat the clustering procedure with the reduced set of features
        
        clusters = zeros( 1, n1 );
        run = 0;
        
        for j = 1:length( set1 )
            
            var = index1( set1( j ) );
            clusters( 1, data( :, var ) == 1 ) = clusters( 1, data( :, var ) == 1 ) + 2^run;
            
            run = run + 1;
            
        end
        
        clusters = clusters + 1;
        
        % Remove empty clusters
        
        index = 1;
        m = max( clusters );
        
        while index <= m
            
            temp = sum( clusters == index );
            
            if temp == 0
                
                temp = find( clusters > index );
                clusters( 1, temp ) = clusters( 1, temp ) - 1;
                
                m = m - 1;
                
            else
                
                index = index + 1;
                
            end
            
        end
        
    end
    
      
    %% If modifying the set of relevant features didn't yield any changes, keep the original set
    
    if max( clusters ) == mcl
        
        set = find( relevantvariables <= alpha );
        
        set2 = set;
        
        count = 2;
        
        while count <= length( set2 )
            
            temp = abs( data( :, index1( set2( count ) ) ) );
            
            for j = 1:count-1
                
                temp2 = abs( data( :, index1( set2( j ) ) ) );
                
                if isequal( temp, temp2 )
                    
                    set2 = set2( set2 ~= set2( count ) );
                    count = count - 1;
                    
                    break;
                    
                end
                
            end
            
            count = count + 1;
            
        end
        
        set = index1( set2 );
        
    else
        
        set = index1( set1 );
        
    end
    
else
    
    set = index1( set1 );
    
end


end