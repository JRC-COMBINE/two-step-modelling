function [ model_training, modelb, auc, modelstats, cutoff, normfunction ] = model_postprocess( model_training, response_training )


%% 1. Normalize the model

if max( model_training )

    normfunction = @(x) ( x - nanmin( model_training ) )./nanmax( model_training ); %#ok<*NANMAX,*NANMIN> 
    model_training = ( model_training - nanmin( model_training ) )./ nanmax( model_training );
    
else
    
    normfunction = @(x) x - nanmin( model_training );
    model_training = model_training - nanmin( model_training );
    
end


%% 2. Binarize the results and assess model performance

modelstats = zeros( 1, 101 );
index = isnan( model_training );

for i = 1:101
    
    cutoff = ( i - 1 ) / 100;
    
    modelb = model_training >= cutoff;
    modelb = +modelb;
    modelb( index ) = NaN;
    
    [ acc, ~, ~, ~, ~ ] = statev( modelb, response_training );
    
    modelstats( 1, i ) = acc;

end

[ ~, i ] = max( modelstats( 1, : ) );

cutoff = ( i - 1 )/100;

modelb = model_training >= cutoff;
modelb = +modelb;
modelb( index ) = NaN;

[ ~, ~, ~, auc ] = perfcurve( response_training, model_training, 1 );

[ acc, prec, recall, f1, FDR ] = statev( modelb, response_training );

modelstats = zeros( 5, 1 );

modelstats( 1, 1 ) = acc;
modelstats( 2, 1 ) = prec;
modelstats( 3, 1 ) = recall;
modelstats( 4, 1 ) = f1;
modelstats( 5, 1 ) = FDR;


end