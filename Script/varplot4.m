function [ completesets, importance ] = varplot4( drug, response, data )

%% load the data

str = strcat( 'res_', num2str( drug ), '.mat' );

load( str, 'var_sets', 'coeffs' );
load( 'GeneOrder.mat' ); %#ok<*LOAD>
load( 'TissueOrder.mat' );
load( 'DrugOrder.mat' );


%% get the feature sets

fsets = cell( 4, 2 );

fsets{ 1, 1 } = unique( cell2mat( var_sets( :, 1 )' ) )';

for l = 2:3
    
    fsets{ l, 1 } = unique( cell2mat( var_sets( :, l ) ) );
    
end

fsets{ 4, 1 } = unique( cell2mat( var_sets( :, 4 )' ) )';


% count how often features are found over all 10 folds

for l = 1:4
    
    temp = fsets{ l, 1 };
    counts = zeros( 1, length( temp ) );
    
    if l == 1 || l == 4
        
        for j = 1:length( temp )
            
            temptemp = temp( j );
            counts( 1, j ) = sum( cell2mat( var_sets( :, l )' ) == temptemp );
            
        end
        
    else
        
        for j = 1:length( temp )
            
            temptemp = temp( j );
            counts( 1, j ) = sum( cell2mat( var_sets( :, l ) ) == temptemp );
            
        end
        
    end
    
    fsets{ l, 2 } = counts;
    
end


%% note whether the features induce sensitivity (-1) or resistance (1)

for l = 1:4
    
    temp = fsets{ l, 1 };
    temptemp = NaN( 1, length( temp ) );
    
    for j = 1:length( temptemp )
        
        temptemp( 1, j ) = sign( nanmean( response( abs( data{ 1, l }( :, temp( j ) ) ) == 1 ) ) - nanmean( response( data{ 1, l }( :, temp( j ) ) == 0 ) ) );
        
    end
    
    fsets{ l, 4 } = temptemp;
    
end


%% sort the features with respect to occurrence

for l = 1:4
    
    [ val, ind ] =sort( fsets{ l, 2 }, 'descend' );
    
    fsets{ l, 2 } = val;
    fsets{ l, 1 } = fsets{ l, 1 }( ind );
    fsets{ l, 3 } = GeneOrder( fsets{ l, 1 } );
    fsets{ l, 4 } = fsets{ l, 4 }( ind );
    
end

fsets{ 4, 3 } = TissueOrder( fsets{ 4, 1 } );


%% calculate the mean importance of the feature sets

meancoeffs = zeros( 6, 11 ); % no NN, no Bayes

for i = 2:12
    
    temp = coeffs( :, i );
    tempmat = NaN( 10, 6 );
    
    for j = 1:10
        
        temptemp = temp{ j, 1 };
        
        if sum( coeffs{ j, 1 } ) == length( temp{ j, 1 } ) % no constant
            
            tempmat( j, coeffs{ j, 1 } == 1 ) = temptemp;
            
        else % constant needs to be excluded
            
            tempmat( j, coeffs{ j, 1 } == 1 ) = temptemp( 2:length( temptemp ) );
            
        end
        
    end
    
    meancoeffs( :, i-1 ) = nanmean( abs( tempmat ) );
    
end

meancoeffs_n = meancoeffs./nanmax( meancoeffs );


%% determine the index

index = sum( cell2mat( coeffs( :, 1 ) ) ) >= 5;
index = index( 1:4 );

datatype = { ' Somatic Mutation ', ' CNV ', ' Hypermethylation ', ' Tissue Descriptors ' };
datatype = datatype( index );

fsets = fsets( index, : );

meancoeffs_n = meancoeffs_n( 1:4, : );
meancoeffs_n = meancoeffs_n( index, : );

importance = nanmean( meancoeffs_n, 2 );
errors = nanstd( meancoeffs_n, 1, 2 );

[ ~, index ] = sort( importance, 'descend' );

datatype = datatype( index );
fsets = fsets( index, : );

importance = importance( index );
errors = errors( index );

completesets = cell( 1, 5 );
count = 0;

for k = 1:size( fsets, 1 )
   
    temp = fsets{ k, 2 };
    temp2 = fsets{ k, 4 };

    completesets( count + 1:count + size( temp, 2 ), 1 ) = fsets{ k, 3 };
    
    for l = 1:size( temp, 2 )
    
        completesets( count + l, 2 ) = datatype( 1, k ); 
        completesets{ count + l, 3 } = temp( l );
        completesets{ count + l, 4 } = temp2( l ) * importance( k );
        completesets{ count + l, 5 } = errors( k );
    
    end
    
    count = count + size( temp, 2 );
    
end

%% create the bar order and the horizontal placing

count = 1; 

temp = cell2mat( completesets( :,3 ) );
indices = [];
barheights = [];
bartypes = [];
barlabels = {};
barfolds = [];
barerrors = [];

for k = 1:10
    
    temp1 = find( temp == k );
    
    if ~isempty( temp1 )
        
        for l = 1:size( datatype, 2 )
            
            ind = strcmpi( datatype{ 1, l }, completesets( temp1, 2 ) );
            
            if sum( ind )
                
                indices = [ indices, temp1( ind )' ];
                
                barheights = [ barheights, cell2mat( completesets( temp1( ind ), 4 ) )' ];
                barerrors = [ barerrors, cell2mat( completesets( temp1( ind ), 5 ) )' ];
                bartypes = [ bartypes, l * ones( 1, sum( ind ) ) ];
                
                barlabels = [ barlabels, completesets{ temp1( ind ), 1 } ];
                
            end
            
        end

        count = count + length( temp1 ) + 1;
        
        barheights = [ barheights, 0 ]; 
        bartypes = [ bartypes, NaN ];
        barlabels = [ barlabels, ' ' ];
        barerrors = [ barerrors, NaN ]; 
        
        barfolds = [ barfolds, k * ones( 1, length( temp1 ) ), 0 ];
        
    else
        
        barheights( 1, count ) = 0; %#ok<*AGROW>
        bartypes( 1, count ) = NaN;
        barlabels = [ barlabels, ' ' ];
        barerrors = [ barerrors, NaN ]; 
        
        barfolds = [ barfolds, 0 ];
        
        count = count + 1;
    
    end
    
end


%% plot

colors = [ 0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560 ];


figure;

b = bar( barheights, 1 );

for k = 1:size( datatype, 2 )

    temptype = find( bartypes == k );
    
    if ~isempty( temptype )
        
        for l = 1:length( temptype )
            
            b.FaceColor = 'flat';
            b.CData( temptype( l ), : ) = colors( k, : );
            
        end
        
    end

end


hold on;


% add errorbars

er = errorbar( 1:length( barerrors ), barheights, barerrors/2, barerrors/2, 'HandleVisibility', 'off' );    
er.Color = [ 0 0 0 ];                            
er.LineStyle = 'none';  


% add invisible plots to correct legend

for k = 2:size( datatype, 2 )

    bar( zeros( 1, length( bartypes ) - 1 ), 'FaceColor', colors( k, : ) )

end


% add brackets and labels to brackets

for k = 1:10
    
    temp = find( barfolds == k );
    
    if ~isempty( temp )
    
        plot( [ min( temp ) - 0.5, max( temp ) + 0.5 ], [ -1, -1 ], 'k', 'HandleVisibility', 'off' );
        plot( [ min( temp ) - 0.5, min( temp ) - 0.5 ], [ -1, -0.95 ], 'k', 'HandleVisibility', 'off' );
        plot( [ max( temp ) + 0.5, max( temp ) + 0.5 ], [ -1, -0.95 ], 'k', 'HandleVisibility', 'off' );
        
        text( mean( temp ), -1.07, num2str( k ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12 )
    
    end
    
    
end


% vertical line for stability

xpos1 = find( barfolds >= 7 );
xpos2 = find( barfolds <= 6 & barfolds ~= 0 );

if ~ isempty( xpos1 ) && ~isempty( xpos2 )

    plot( [ mean( [ min( xpos1 ), max( xpos2 ) ] ), mean( [ min( xpos1 ), max( xpos2 ) ] ) ], [ -1.2, 1.2 ], '-.k', 'HandleVisibility', 'off' )

    text( mean( [ min( xpos1 ), max( xpos2 ) ] ) - length( barheights )/8.5, 1.1, '\leftarrow Unstable Features', 'FontSize', 12 )
    text( mean( [ min( xpos1 ), max( xpos2 ) ] ) + 0.5, 1.1, 'Stable Features \rightarrow', 'FontSize', 12 )
    
end

hold off;


legend( datatype, 'Location', 'northwest', 'FontSize', 12 )

set( gca, 'FontSize', 12 )

xlim( [ 0, length( bartypes ) ] );
xticks( 1:length( bartypes ) - 1 );
xticklabels( barlabels )
xtickangle( 35 )

ylim( [ -1.2, 1.2 ] )
yticks( -1.2:0.1:1.2 )
yticklabels( { 'Sensitivity', ' ', '-1', ' ', '-0.8', ' ', '-0.6', ' ', '-0.4', ' ', '-0.2', ' ', '0', ' ', '0.2', ' ', '0.4', ' ', '0.6', ' ', '0.8', ' ', '1', ' ', 'Resistance' } )

grid on;
box on;    


end