function [ amodels, trainingstats, teststats ] = ablationstudies( model, model_testing, responseb_training, responseb_test, index, indexl, type )


combi2 = nchoosek( 1:6, 2 );
combi3 = nchoosek( 1:6, 3 );

combinations = cell( 41, 1 );

for j = 1:6
   
    combinations{ j, 1 } = j;
    
end

for j = 1:size( combi2, 1 )
   
    combinations{ j+6, 1 } = combi2( j, : );
    
end

for j = 1:size( combi3, 1 )
   
    combinations{ j+6+size( combi2, 1 ), 1 } = combi3( j, : );
    
end


%% Neural Network

if type == 1
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ]; %#ok<*AGROW> 
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            anet = fitnet( 5, 'trainbr' );
            
            % Setup Division of Data for Training, Validation, Testing
            
            anet.divideParam.trainRatio = 80/100;
            anet.divideParam.valRatio = 0/100;
            anet.divideParam.testRatio = 20/100;
            
            % Train the Network
            
            [ anet, ~ ] = train( anet, amodel', responseb_training' );
            
            % Apply the Network
            
            model_training = anet( amodel' );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( model_training', responseb_training );
            
            amodels{ 1, k } = anet;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_NN( anet, amodel_test, cutoff, responseb_test, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 
    

%% Linear regression with LASSO regularization

if type == 2
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lasso( amodel, responseb_training );
            
            [ ~, ind ] = min( fitinfo.MSE );
            coef = model_training( :, ind );
            constant = fitinfo.Intercept( ind );
            y_l = amodel*coef + constant;
            
            coeffs_l = [ constant; coef ];
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_l, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [~, ~, teststat, auc ] = test_model_linreg( coeffs_l, amodel_test, cutoff, responseb_test, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Linear regression with elastic net regularization

if type == 3 
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lasso( amodel, responseb_training, 'Alpha', 0.1 );
            
            [ ~, ind ] = min( fitinfo.MSE );
            coef = model_training( :, ind );
            constant = fitinfo.Intercept( ind );
            y_en = amodel*coef + constant;
            
            coeffs_en = [ constant; coef ];
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_en, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [~, ~, teststat, auc ] = test_model_linreg( coeffs_en, amodel_test, cutoff, responseb_test, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Linear regression with ridge regularization

if type == 4 
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lasso( amodel, responseb_training, 'Alpha', 10^(-3) );
            
            [ ~, ind ] = min( fitinfo.MSE );
            coef = model_training( :, ind );
            constant = fitinfo.Intercept( ind );
            y_r = amodel*coef + constant;
            
            coeffs_r = [ constant; coef ];
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_r, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [~, ~, teststat, auc ] = test_model_linreg( coeffs_r, amodel_test, cutoff, responseb_test, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
      
end 


%% Logistic regression with LASSO regularization

if type == 5 
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lassoglm( amodel, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5 );
            
            constant = fitinfo.Intercept( fitinfo.Index1SE );
            coeffs = model_training( :, fitinfo.Index1SE );
            coefficients = [ constant; coeffs ];
            
            y_log1 = glmval( coefficients, amodel, 'logit' );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_log1, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_log( coefficients, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Logistic regression with elastic net regularization

if type == 6 
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lassoglm( amodel, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5, 'Alpha', 0.5 );
            
            constant = fitinfo.Intercept( fitinfo.Index1SE );
            coeffs = model_training( :, fitinfo.Index1SE );
            coefficients = [ constant; coeffs ];
            
            y_log2 = glmval( coefficients, amodel, 'logit' );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_log2, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_log( coefficients, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Logistic regression with ridge regularization

if type == 7 
    
    amodels = cell( 2, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, fitinfo ] = lassoglm( amodel, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5, 'Alpha', 10^(-3) );
            
            constant = fitinfo.Intercept( fitinfo.Index1SE );
            coeffs = model_training( :, fitinfo.Index1SE );
            coefficients = [ constant; coeffs ];
            
            y_log3 = glmval( coefficients, amodel, 'logit' );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_log3, responseb_training );
            
            amodels{ 1, k } = model_training;
            amodels{ 2, k } = fitinfo;
            
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_log( coefficients, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% SVM with ridge regularization

if type == 8 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, ~ ] = fitrlinear( amodel, responseb_training, 'Learner', 'svm', 'Regularization', 'ridge' );
            
            y_svm1 = predict( model_training, amodel );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_svm1, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_svm( model_training, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% SVM with LASSO regularization

if type == 9 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ model_training, ~ ] = fitrlinear( amodel, responseb_training, 'Learner', 'svm', 'Regularization', 'lasso' );
            
            y_svm2 = predict( model_training, amodel );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_svm2, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_svm( model_training, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Naïve Bayes classifier

if type == 10 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 1, 41 );
    teststats = NaN( 1, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            model_training = fitcnb( amodel, responseb_training );
            
            y_nbayes = predict( model_training, amodel );
            
            [ ~, ~, ~, trainingstat, ~, ~ ] = model_postprocess( y_nbayes, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, teststat ] = test_model_nbayes( model_training, amodel_test, responseb_test );
            
            teststats( 1, k ) = teststat( 1, 1 );
            
        end
        
    end
       
end 


%% Bagged decision tree ensemble

if type == 11 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            treeStump = templateTree( 'MaxNumSplits', 1 );
    
            model_training = fitrensemble( amodel, responseb_training, 'Method', 'Bag', 'Learners', treeStump );
            
            y_treebag = predict( model_training, amodel );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_treebag, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_dectrees( model_training, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Boosted decision tree ensemble

if type == 12 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            treeStump = templateTree( 'MaxNumSplits', 1 );
    
            model_training = fitrensemble( amodel, responseb_training, 'Method', 'LSBoost', 'Learners', treeStump );
            
            y_treeboost = predict( model_training, amodel );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_treeboost, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_dectrees( model_training, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


%% Random forest

if type == 13 
    
    amodels = cell( 1, 41 );
    trainingstats = NaN( 2, 41 );
    teststats = NaN( 2, 41 );
    
    for k = 1:41
        
        if prod( index( combinations{ k, 1 } ) ) && sum( index ) > length( index( combinations{ k, 1 } ) )
            
            remove = [];
            
            for l = 1:length( index( combinations{ k, 1 } ) )
                
                remove = [ remove, find( indexl == combinations{ k, 1 }( l ) ) ];
                
            end
            
            amodel = model( :, setdiff( 1:size( model, 2 ), remove ) );
            
            model_training = TreeBagger( 50, amodel, responseb_training, 'Method', 'Regression', 'MinLeafSize', 20, 'OOBPredictorImportance','on', 'PredictorSelection', 'curvature' );
            
            y_rforest = predict( model_training, amodel );
            
            [ ~, ~, auc, trainingstat, cutoff, normfunction ] = model_postprocess( y_rforest, responseb_training );
            
            amodels{ 1, k } = model_training;
            trainingstats( 1, k ) = trainingstat( 1, 1 );
            trainingstats( 2, k ) = auc;
            
            amodel_test = model_testing( :, setdiff( 1:size( model, 2 ), remove ) );
            
            [ ~, ~, teststat, auc ] = test_model_rforest( model_training, amodel_test, responseb_test, cutoff, normfunction );
            
            teststats( 1, k ) = teststat( 1, 1 );
            teststats( 2, k ) = auc;
            
        end
        
    end
    
end 


end