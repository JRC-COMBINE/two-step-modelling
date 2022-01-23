function [ b, b2, completesets ] = varplot6vert( completesets, drug, data )

%% load the data

str = strcat( 'res_', num2str( drug ), '.mat' );

load( str, 'var_sets' );
load( 'GeneOrder.mat' ); %#ok<*LOAD>
load( 'TissueOrder.mat' );
load( 'DrugOrder.mat' );


for k = 1:10
    
    var_sets{ k, 1 } = var_sets{ k, 1 }';
    var_sets{ k, 5 } = var_sets{ k, 5 }';
    
end

%% find the additional features, that are present in the extended set and absent in the normal set

ext = cell( 1, 3 );

% mutations

temp = unique( cell2mat( var_sets( :, 5 ) ) )';
ext{ 1, 1 } = setdiff( temp, unique( cell2mat( var_sets( :, 1 ) ) )' );

% CNVs

temp = unique( cell2mat( var_sets( :, 6 ) ) )';
ext{ 1, 2 } = setdiff( temp, unique( cell2mat( var_sets( :, 2 ) ) )' );

% methylations

temp = unique( cell2mat( var_sets( :, 7 ) ) )';
ext{ 1, 3 } = setdiff( temp, unique( cell2mat( var_sets( :, 3 ) ) )' );


%% for each new feature, find which other feature is distributed the same way, adapt the missing values from the paired feature and update completesets

for k = 1:3
    
    if ~isempty( ext{ 1, k } )
        
        temp = ext{ 1, k };
        
        for l = 1:length( temp )
            
            temp1 = unique( cell2mat( var_sets( :, k ) ) )';
            temp2 = data{ 1, k };
            
            temp3 = temp2( :, temp1 );
            temp2 = temp2( :, temp( l ) );
            
            idx = temp1( sum( temp3 - temp2 ) == 0 ); % index of the corresponding feature in the list
            
            if ~isempty( idx )
            
                name = GeneOrder{ idx, 1 }; %#ok<*USENS>
                newname = GeneOrder{ temp( l ), 1 };
                
                ind = find( strcmp( name, completesets( :, 1 ) ) == 1 );
                
                tempsets = cell( size( completesets, 1 ) + 1, 5 );
                
                tempsets( 1:ind, : ) = completesets( 1:ind, : );
                
                tempsets( ind + 1, : ) = completesets( ind, : );
                tempsets{ ind + 1, 1 } = newname;
                
                tempsets( ind + 2:size( tempsets, 1 ), : ) = completesets( ind + 1:size( completesets, 1  ), : );
                
                completesets = tempsets;
            
            end
            
        end
        
    end
    
end


%% find the importances and sort accordingly,also check which data types are present

datatype = unique( completesets( :, 2 ) );
importance = zeros( size( datatype ) );

for k = 1:size( datatype, 1 )
    
    importance( k, 1 ) = unique( abs( cell2mat( completesets( strcmp( datatype{ k, 1 }, completesets( :, 2 ) ), 4 ) ) ) );
    
end


[ ~, index ] = sort( importance, 'descend' );

datatype = datatype( index )';


%% create the bar order and the horizontal placing

count = 1; 

temp = cell2mat( completesets( :, 3 ) );
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

%colors = [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880 ]; % green
colors = [ 0, 0.4470, 0.7410; 0.6350, 0.0780, 0.1840; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560 ];

figure;

b = barh( barheights, 1, 'FaceAlpha', 0.75 );

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
  
er = errorbar( barheights, 1:length( barerrors ),  barerrors/2, 'horizontal', 'HandleVisibility', 'off' );    
er.Color = [ 0 0 0 ];                            
er.LineStyle = 'none';  


% add invisible plots to correct legend 

b2 = [];

for k = 1:size( datatype, 2 )

    b2( 1, k ) = barh( zeros( 1, length( bartypes ) - 1 ), 'FaceColor', colors( k, : ), 'FaceAlpha', 0.8 );

end


% add brackets and labels to brackets

for k = 1:10
    
    temp = find( barfolds == k );
    
    if ~isempty( temp )
    
        plot( [ -1, -1 ], [ min( temp ) - 0.5, max( temp ) + 0.5 ], 'k', 'HandleVisibility', 'off' );
        plot( [ -1, -0.95 ], [ min( temp ) - 0.5, min( temp ) - 0.5 ], 'k', 'HandleVisibility', 'off' );
        plot( [ -1, -0.95 ], [ max( temp ) + 0.5, max( temp ) + 0.5 ], 'k', 'HandleVisibility', 'off' );
        
        text( -1.07,  mean( temp ) - 0.75, num2str( k ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14 )
    
    end
    
    
end


% horizontal line for stability

ypos1 = find( barfolds >= 7 );
ypos2 = find( barfolds <= 6 & barfolds ~= 0 );

if ~ isempty( ypos1 ) && ~isempty( ypos2 )
    
    ttemp = mean( [ min( ypos1 ), max( ypos2 ) ] );
    plot( [ -1.2, 1.2 ], [ ttemp, ttemp ], '-.k', 'HandleVisibility', 'off' )

    text( 1.25, mean( [ min( ypos1 ), max( ypos2 ) ] ) - 0.5, ...
        {'\downarrow'; 'Unstable'; 'Features'}, 'FontSize', 14, 'HorizontalAlignment', ...
        'left', 'VerticalAlignment', 'top' )
    text( 1.25, mean( [ min( ypos1 ), max( ypos2 ) ] ) + 0.8, {'Stable'; 'Features'; '\uparrow' },...
        'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom' )
    
end

text( 0.05, length( bartypes ) + 1.5, ' \rightarrow Resistance', 'FontSize', 14 )
text(  -0.05, length( bartypes ) + 1.5, 'Sensitivity \leftarrow', 'FontSize', 14, 'HorizontalAlignment', 'right' )

hold off;


legend( b2, datatype, 'Location', 'northeast', 'FontSize', 14 )

set( gca, 'FontSize', 12 )

ylim( [ 0, length( bartypes ) ] );
yticks( 1:length( bartypes ) - 1 );
yticklabels( barlabels )
ytickangle( 25 )

xlim( [ -1.2, 1.2 ] )
xticks( -1.2:0.1:1.2 )
xticklabels( { ' ', ' ', '-1', ' ', '-0.8', ' ', '-0.6', ' ', '-0.4', ' ', '-0.2', ' ', '0', ' ', '0.2', ' ', '0.4', ' ', '0.6', ' ', '0.8', ' ', '1', ' ', ' ' } )
xlabel( 'Mean relative importance', 'FontSize', 14 )

grid on;
box on;  


end