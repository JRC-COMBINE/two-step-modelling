function [ h, t ] = featureplot( freqs )

% preprocess the required information

% mutations

temp = unique( freqs{ 4, 1 } );
mutfeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 1 };
    temp2 = find( freqs{ 4, 1 } == temp( k ) );
    
    mutfeatures{ 1, k } = temp( k );
    mutfeatures{ 2, k } = length( temp2 );
    mutfeatures{ 3, k } =  temp1( temp2 );
    
end

% CNVs

temp = unique( freqs{ 4, 2 } );
cnvfeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 2 };
    temp2 = find( freqs{ 4, 2 } == temp( k ) );
    
    cnvfeatures{ 1, k } = temp( k );
    cnvfeatures{ 2, k } = length( temp2 );
    cnvfeatures{ 3, k } =  temp1( temp2 );
    
end

% methylations

temp = unique( freqs{ 4, 3 } );
methfeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 3 };
    temp2 = find( freqs{ 4, 3 } == temp( k ) );
    
    methfeatures{ 1, k } = temp( k );
    methfeatures{ 2, k } = length( temp2 );
    methfeatures{ 3, k } =  temp1( temp2 );
    
end

% tissues

temp = unique( freqs{ 4, 4 } );
tissuefeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 4 };
    temp2 = find( freqs{ 4, 4 } == temp( k ) );
    
    tissuefeatures{ 1, k } = temp( k );
    tissuefeatures{ 2, k } = length( temp2 );
    tissuefeatures{ 3, k } =  temp1( temp2 );
    
end

% pathways

temp = unique( freqs{ 4, 5 } );
pwfeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 5 };
    temp2 = find( freqs{ 4, 5 } == temp( k ) );
    
    pwfeatures{ 1, k } = temp( k );
    pwfeatures{ 2, k } = length( temp2 );
    pwfeatures{ 3, k } =  temp1( temp2 );
    
end

% PCs

temp = unique( freqs{ 4, 6 } );
pcfeatures = cell( 3, length( temp ) );

for k = 1:length( temp )
    
    temp1 = freqs{ 5, 6 };
    temp2 = find( freqs{ 4, 6 } == temp( k ) );
    
    pcfeatures{ 1, k } = temp( k );
    pcfeatures{ 2, k } = length( temp2 );
    pcfeatures{ 3, k } =  temp1( temp2 );
    
end


count = 1;
t = [];

h = figure;
hold on;

for k = 1:size( mutfeatures, 2 )

    plot( mutfeatures{ 1, k }, 1, 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 6 + mutfeatures{ 2, k } )
    
    if length( mutfeatures{ 3, k } ) <= 6
        
        text( mutfeatures{ 1, k } + 1, 1, mutfeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( mutfeatures{ 1, k } + 2, 1, [ '#', '_', num2str( count ) ] )
        count = count + 1;
        
    end
    
end


for k = 1:size( cnvfeatures, 2 )

    plot( cnvfeatures{ 1, k }, 2, 'o', 'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6 + cnvfeatures{ 2, k } )
    
    if length( cnvfeatures{ 3, k } ) <= 6
        
        text( cnvfeatures{ 1, k } + 1, 2, cnvfeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( cnvfeatures{ 1, k } + 2, 2, [ '#', '_', num2str( count ) ] )
        
        dim = [ 0.2, 0.5, 0.3, 0.3 ];
        str = [ '#_', num2str( count ), cnvfeatures{ 3, k } ] ;
        annotation( 'textbox', dim, 'String', str, 'FitBoxToText', 'on' );
         
        count = count + 1;
        
    end

end


for k = 1:size( methfeatures, 2 )

    plot( methfeatures{ 1, k }, 3, 'o', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], ...
        'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerSize', 6 + methfeatures{ 2, k } )
    
    if length( methfeatures{ 3, k } ) <= 6
        
        text( methfeatures{ 1, k } + 1, 3, methfeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( methfeatures{ 1, k } + length( methfeatures{ 3, k } )/10, 3, [ '#', '_', num2str( count ) ] )
        
        dim = [ 0.83, 0.45, 0.3, 0.3 ];
        str = [ [ '#_', num2str( count ) ]; methfeatures{ 3, k } ] ;
        tt = annotation( 'textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'black', 'FaceAlpha', .1 );
        t = [ t, tt ];
        
        count = count + 1;
        
    end

end


for k = 1:size( tissuefeatures, 2 )

    plot( tissuefeatures{ 1, k }, 4, 'o', 'MarkerEdgeColor', [0.4940 0.1840 0.5560], ...
        'MarkerFaceColor', [0.4940 0.1840 0.5560], 'MarkerSize', 6 + tissuefeatures{ 2, k } )
    
    if length( tissuefeatures{ 3, k } ) <= 6
        
        text( tissuefeatures{ 1, k } + 1, 4, tissuefeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( tissuefeatures{ 1, k } + 2, 4, [ '#', '_', num2str( count ) ] )
        count = count + 1;
        
    end

end


for k = 1:size( pwfeatures, 2 )

    plot( pwfeatures{ 1, k }, 5, 'o', 'MarkerEdgeColor', [0.4660 0.6740 0.1880], ...
        'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerSize', 6 + pwfeatures{ 2, k } )
    
    if length( pwfeatures{ 3, k } ) <= 6
        
        text( pwfeatures{ 1, k } + 1, 5, pwfeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( pwfeatures{ 1, k } + 2, 5, [ '#', '_', num2str( count ) ] )
        count = count + 1;
        
    end

end


for k = 1:size( pcfeatures, 2 )

    plot( pcfeatures{ 1, k }, 6, 'o', 'MarkerEdgeColor', [0.8500 0.3250 0.0980], ...
        'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 6 + pcfeatures{ 2, k } )
    
    if length( pcfeatures{ 3, k } ) <= 6
        
        text( pcfeatures{ 1, k } + 1, 6, pcfeatures{ 3, k }, 'FontSize', 12 )
        
    else
        
        text( pcfeatures{ 1, k } + 2, 6, [ '#', '_', num2str( count ) ] )
        count = count + 1;
        
    end

end

% adjust widths of textboxes

temp = [];

for k = 1:length( t )
   
    temp = [ temp, t( k ).Position( 3 ) ];
    
end

m = max( temp );

for k = 1:length( t )
   
    t( k ).Position( 3 ) = m;
    
end


grid on; box on;
ylim( [ 0, 7 ] )
xlim( [ 0, 110 ] )
xticks( 10:10:100 )

set( gca, 'FontSize', 12 )
xlabel( 'Occurrence of stable, significant features [%]' )

yticks( 1:6 )
yticklabels( { 'Somatic mutation', 'CNV', 'Hypermethylation', 'Tissue descriptors', 'Pathway activation', 'Gene expression' } )
ytickangle( 25 )


end