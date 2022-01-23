function [ accuracy, precision, recall, f1, FDR ] = statev( model, measurement )

% Compute metrics to evaluate the classification performance
% 1 represents responsive cell lines
% 0 represents resistant cell lines

tp = nansum( model == 1 & measurement == 1 ); %#ok<*NANSUM> 
tn = nansum( model == 0 & measurement == 0 );
fp = nansum( model == 1 & measurement == 0 );
fn = nansum( model == 0 & measurement == 1 );

accuracy = ( tp + tn )/( tp + tn + fp + fn );

precision = tp/( tp + fp ); 
recall = tp/( tp + fn );

f1 = 2 * ( precision * recall )/( precision + recall );

FDR = fp/( fp + tp );

end