function [ submodel, normfunction ] = normsubmodels( submodel )

% Normalize the outputs of the models on continuous data

if ~isnan( nanmin( submodel ) ) && ~isnan( nanmax( submodel ) )
    
    if nanmax( submodel ) - nanmin( submodel ) ~= 0 
       
        normfunction = @(x) ( x - nanmin( submodel ) )./( nanmax( submodel ) - nanmin( submodel ) );
        submodel = ( submodel - nanmin( submodel ) )./( nanmax( submodel ) - nanmin( submodel ) );
        
    else
        
        normfunction = @(x) x - nanmin( submodel );
        submodel = submodel - nanmin( submodel );
        
    end
    
end


end