function h = stepplot2( freqs, ind )

temp = cell2mat( freqs( 2, : ) );
temp2 = unique( temp );
perc2 = zeros( size( temp2 ) );

%custyticks = { 0:20:80, 0:20:100, 0:20:120, 0:20:140, 0:20:120, 0:10:20, 0:20:80, 0:20:120, 0:20:160, ...
%0:20:40, 0:20:40, 0:20:80, 0:20:100, 0:20:40, 0:20:120, 0:20:40, 0:20:100, 0:20:40, 0:20:60, 0:20:40 };

ylims = [ 80, 100, 120, 140, 120, 20, 80, 120, 160, 40, 40, 80, 100, 40, 120, 40, 100, 40, 60, 40 ];

for k = 1:length( temp2 )
   
    perc2( k ) = length( find( temp >= temp2( k ) ) );
    
end

temp2 = [ 0, temp2 ];
perc2 = [ perc2, 0 ];

h = 0;

hold on;

for k = 1:length( temp2 ) - 1 
    
    patch( [ temp2( k ), temp2( k ), temp2( k + 1 ), temp2( k + 1 ) ], [ 0, perc2( k ), perc2( k ), 0 ], ...
         [ 0 0.4470 0.7410 ], 'EdgeColor', 'none', 'FaceAlpha', '0.5' )
    
    plot( [ temp2( k ), temp2( k + 1 ) ], [ perc2( k ), perc2( k ) ], 'Color', ...
        [0 0.4470 0.7410], 'LineWidth', 1.5  )
    
    plot( [ temp2( k + 1 ), temp2( k + 1 ) ], [ perc2( k ), perc2( k + 1 ) ], ...
        'Color', [0 0.4470 0.7410], 'LineWidth', 1.5 )

end

set( gca, 'FontSize', 12 )

grid on;
box on;
xlim( [ 0, 100 ] )
xticks( 0:20:100 )

ylim( [ 0, ylims( ind ) ] )

end