function [ models, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part ] = twostepmodel( drug )


load( 'SomaticMutation.mat', 'SomaticMutation' );
load( 'CNV.mat', 'CNV' );
load( 'Hypermethylation.mat', 'Hypermethylation' );
load( 'TissueType.mat', 'TissueType' );
load( 'GeneExpression.mat', 'GeneExpression' );
load( 'PathwayActivation.mat', 'PathwayActivation' );
load( 'Response_AUC.mat', 'Response_AUC' );


models = cell( 10, 12 );

prediction_training = cell( 10, 28 );
prediction_test = cell( 10, 28 );

response = cell( 10, 4 );

trainingstats = cell( 10, 13 );
teststats = cell( 10, 13 );

AUCs = cell( 10, 2 );

var_sets = cell( 10, 7 );
coeffs = cell( 10, 8 );

ablstudies_results = cell( 10, 3 );


%% 1. Compute a training and test set partitition for cross validation

[ n1, ~ ] = size( GeneExpression ) ;

cv_part = cvpartition( n1, 'kFold', 10 );


for j = 1:10
    
    trainingset = training( cv_part, j );
    testset = test( cv_part, j );
    
    [ pc_coeff, score, ~ ] = pca( GeneExpression( trainingset, : ) );
    PC_training = score( :, 1:7 );
    
    % Project test gene expression data onto PC space
    
    PC_test = mrdivide( GeneExpression( testset, : ), pc_coeff' );
    PC_test = PC_test( :, 1:7 );
    
    pathway_training = PathwayActivation( trainingset, : );
    pathway_test = PathwayActivation( testset, : );
    
    [ response_training, response_test, responseb_training, responseb_test ] = preprocess_response( Response_AUC( :, drug ), trainingset, testset );
    
    
    %% 2. Compute discrete feature profiles related to shifts in responsiveness
    
    [ clustdata1, mutset, mutset_l ] = mutclusters( SomaticMutation( trainingset, : ), SomaticMutation, response_training );
    
    clustdata1_training = clustdata1( trainingset );
    clustdata1_test = clustdata1( testset );
    n1 = max( clustdata1_training );
    
    
    [ clustdata2, cnvset, cnvset_l ] = cnvclusters( CNV( trainingset, : ), CNV, response_training );
    
    clustdata2_training = clustdata2( trainingset );
    clustdata2_test = clustdata2( testset );
    n2 = max( clustdata2_training );
    
    
    [ clustdata3, methset, methset_l ] = methclusters( Hypermethylation( trainingset, : ), Hypermethylation, response_training );
    
    clustdata3_training = clustdata3( trainingset );
    clustdata3_test = clustdata3( testset );
    n3 = max( clustdata3_training );
    
    
    [ clustdata4, tissueset ] = tissueclusters( TissueType( trainingset, : ), TissueType, response_training );
    
    clustdata4_training = clustdata4( trainingset );
    clustdata4_test = clustdata4( testset );
    n4 = max( clustdata4_training );
    
    
    %% 3. Compute the first-step models separately on the training data set
    
    [ modelfun1, model1_training, model1b_training, cutoff1, trainingstats1, auc1 ] = mutmodel( responseb_training, clustdata1_training) ;
    [ modelfun2, model2_training, model2b_training, cutoff2, trainingstats2, auc2 ] = CNVmodel( responseb_training, clustdata2_training) ;
    [ modelfun3, model3_training, model3b_training, cutoff3, trainingstats3, auc3 ] = methmodel( responseb_training, clustdata3_training) ;
    [ modelfun4, model4_training, model4b_training, cutoff4, trainingstats4, auc4 ] = tissuemodel( responseb_training, clustdata4_training );
    [ modelfun5, model5_training, model5b_training, cutoff5, trainingstats5, auc5 ] = pathwaymodel( responseb_training, pathway_training );
    [ modelfun6, model6_training, model6b_training, cutoff6, trainingstats6, auc6 ] = PCmodel( responseb_training, PC_training );
    
    
    % Run the first-step models on the test set
    
    [ model1_test, model1b_test, teststats1, auc_test_1 ] = test_mutmodel( modelfun1, clustdata1_test, cutoff1, responseb_test, n1 );
    [ model2_test, model2b_test, teststats2, auc_test_2 ] = test_CNVmodel( modelfun2, clustdata2_test, cutoff2, responseb_test, n2 );
    [ model3_test, model3b_test, teststats3, auc_test_3 ] = test_methmodel( modelfun3, clustdata3_test, cutoff3, responseb_test, n3 );
    [ model4_test, model4b_test, teststats4, auc_test_4 ] = test_tissuemodel( modelfun4, clustdata4_test, cutoff4, responseb_test, n4 );
    [ model5_test, model5b_test, teststats5, auc_test_5 ] = test_pathwaymodel( modelfun5, pathway_test, cutoff5, responseb_test );
    [ model6_test, model6b_test, teststats6, auc_test_6 ] = test_PCmodel( modelfun6, PC_test, cutoff6, responseb_test );
   
    
    %% 4. Train second-step models on the output vectors of the first-step models
    
    model = [ model1_training, model2_training, model3_training, model4_training, model5_training, model6_training ];

    
    % Normalize the regression submodels to [min(response), max(response)]
    
    [ nmodel, nmodelf5 ] = normsubmodels( model( :, 5 ) );
    model( :, 5 ) = nmodel;
    
    [ nmodel, nmodelf6 ] = normsubmodels( model( :, 6 ) );
    model( :, 6 ) = nmodel;
    
    
    % Remove constant submodels
    
    index = ones( 1, 6 );
    
    for k = 1:6
        
        if length( unique( model( :, k ) ) ) == 1
            
            index( 1, k ) = 0;
            
        end
        
    end
    
    model = model( :, index == 1 );
    
    
    %% Neural Network model
    
    trainFcn = 'trainbr';  
    
    hiddenLayerSize = 5;
    model_nn = fitnet( hiddenLayerSize, trainFcn );
    
    model_nn.divideParam.trainRatio = 80/100;
    model_nn.divideParam.valRatio = 0/100;
    model_nn.divideParam.testRatio = 20/100;
    
    
    % Train and apply the network
    
    [ model_nn, ~ ] = train( model_nn, model', responseb_training' );
    
    model_training_nn = model_nn( model' );
    
    [ model_training_nn, modelb_training_nn, auc_nn, trainingstats_nn, cutoff_nn, normfunction_nn ] = model_postprocess( model_training_nn', responseb_training );
    
    
    %% Linear regression models with regularization
    
    [ model_l, fitinfo_l ] = lasso( model, responseb_training );
    [ model_en, fitinfo_en ] = lasso( model, responseb_training, 'Alpha', 0.0001 ); % elastic net regularization
    [ model_r, fitinfo_r ] = lasso( model, responseb_training, 'Alpha', 0.0000000001 ); % approximation of ridge regularization

    
    [ ~, ind ] = min( fitinfo_l.MSE );
    coef = model_l( :, ind );
    coef0 = fitinfo_l.Intercept( ind );
    y_l = model*coef + coef0;
    
    coeffs_l = [ coef0; coef ];
    
    [ ~, ind ] = min( fitinfo_en.MSE );
    coef = model_en( :, ind );
    coef0 = fitinfo_en.Intercept( ind );
    y_en = model*coef + coef0;
    
    coeffs_en = [ coef0; coef ];
    
    [ ~, ind ] = min( fitinfo_r.MSE );
    coef = model_r( :, ind );
    coef0 = fitinfo_r.Intercept( ind );
    y_r = model*coef + coef0;
    
    coeffs_r = [ coef0; coef ];
   
    
    [ model_training_l, modelb_training_l, auc_l, trainingstats_l, cutoff_l, normfunction_l ] = model_postprocess( y_l, responseb_training );
    [ model_training_en, modelb_training_en, auc_en, trainingstats_en, cutoff_en, normfunction_en ] = model_postprocess( y_en, responseb_training );
    [ model_training_r, modelb_training_r, auc_r, trainingstats_r, cutoff_r, normfunction_r ] = model_postprocess( y_r, responseb_training );
    
    
    %% Logistic regression models with regularization
    
    [ model_log1, ~ ] = fitclinear( model, responseb_training, 'Learner', 'logistic', 'Regularization', 'ridge' ); 
    [ model_log2, ~ ] = fitclinear( model, responseb_training, 'Learner', 'logistic', 'Regularization', 'lasso' ); 
    
    y_log1 = predict( model_log1, model );
    y_log2 = predict( model_log2, model );
    
    coeffs_log1 = model_log1.Beta; 
    coeffs_log2 = model_log2.Beta;
   
    [ model_training_log1, ~, ~, trainingstats_log1, ~, ~ ] = model_postprocess( y_log1, responseb_training );
    [ model_training_log2, ~, ~, trainingstats_log2, ~, ~ ] = model_postprocess( y_log2, responseb_training );
    
    
    %% Support vector Machines with regularization
    
    [ model_svm1, ~ ] = fitclinear( model, responseb_training, 'Learner', 'svm', 'Regularization', 'ridge' ); 
    [ model_svm2, ~ ] = fitclinear( model, responseb_training, 'Learner', 'svm', 'Regularization', 'lasso' ); 
 
    y_svm1 = predict( model_svm1, model );
    y_svm2 = predict( model_svm2, model );
    
    coeffs_svm1 = model_svm1.Beta;
    coeffs_svm2 = model_svm2.Beta;
    
    [ model_training_svm1, ~, ~, trainingstats_svm1, ~, ~ ] = model_postprocess( y_svm1, responseb_training );
    [ model_training_svm2, ~, ~, trainingstats_svm2, ~, ~ ] = model_postprocess( y_svm2, responseb_training );
    
    
    %% Naïve Bayes Classifier
    
    index2 = ones( 1, size( model, 2 ) );
    
    
    % Remove first-step models lacking variation in their prediction for sensitive or resistant cell lines
    
    for k = 1:size( model, 2 )
        
        if length( unique( model( responseb_training == 0, k ) ) ) == 1 || length( unique( model( responseb_training == 1, k ) ) ) == 1
            
            index2( 1, k ) = 0;
            
        end
        
    end
    
    model_red = model( :, index2 == 1 );
    
    model_nbayes = fitcnb( model_red, responseb_training ); 
    y_nbayes = predict( model_nbayes, model_red );
    
    [ model_training_nbayes, ~, ~, trainingstats_nbayes, ~, ~ ] = model_postprocess( y_nbayes, responseb_training );
    
    
    %% Decision Trees/Random Forests
    
    treeStump = templateTree( 'MaxNumSplits', 1 );
    
    model_treebag = fitcensemble( model, responseb_training, 'Method', 'Bag', 'Learners', treeStump );
    model_treeboost = fitcensemble( model, responseb_training, 'Method', 'AdaBoostM1','Learners', treeStump );
    model_rforest = TreeBagger( 50, model, responseb_training, 'Method', 'Classification', 'MinLeafSize', 20, 'OOBPredictorImportance','on', 'PredictorSelection', 'curvature' );
    
    y_treebag = predict( model_treebag, model );
    y_treeboost = predict( model_treeboost, model );
    y_rforest = str2num( cell2mat( predict( model_rforest, model ) ) ); %#ok<ST2NM>
    
    [ model_training_treebag, ~, trainingstats_treebag, ~, ~ ] = model_postprocess( y_treebag, responseb_training );
    [ model_training_treeboost, ~, trainingstats_treeboost, ~, ~ ] = model_postprocess( y_treeboost, responseb_training );
    [ model_training_rforest, ~, trainingstats_rforest, ~, ~ ] = model_postprocess( y_rforest, responseb_training );
    
    
    auc_training = [ auc1, auc2, auc3, auc4, auc5, auc6, auc_nn, auc_l, auc_en, auc_r ];
  
    
    %% 5. Apply the second-step models to the test set
    
    model_testing = [ model1_test, model2_test, model3_test, model4_test, model5_test, model6_test ];
    model_testing = model_testing( :, index == 1 );
    
    
    % Apply the normalization routines to continuous first-step models
    
    if index( 5 ) == 1
        
        model_testing( :, sum( index( 1:5 ) ) ) = nmodelf5( model_testing( :, sum( index( 1:5 ) ) ) );
        
    end
    
    
    if index( 6 ) == 1
        
        model_testing( :, sum( index( 1:6 ) ) ) = nmodelf6( model_testing( :, sum( index( 1:6 ) ) ) );
        
    end
    
    
    %% Neural Network
    
    [ model_test_nn, modelb_test_nn, teststats_nn, auc_nn ] = test_model_NN( model_nn, model_testing, cutoff_nn, responseb_test, normfunction_nn );
    
    
    %% Linear regression models with regularization
    
    [ model_test_l, modelb_test_l, teststats_l, auc_test_l ] = test_model_linreg( coeffs_l, model_testing, cutoff_l, responseb_test, normfunction_l );
    [ model_test_en, modelb_test_en, teststats_en, auc_test_en ] = test_model_linreg( coeffs_en, model_testing, cutoff_en, responseb_test, normfunction_en );
    [ model_test_r, modelb_test_r, teststats_r, auc_test_r ] = test_model_linreg( coeffs_r, model_testing, cutoff_r, responseb_test, normfunction_r );

    
    %% Logistic regression models with regularization
    
    [ model_test_log1, teststats_log1 ] = test_model_logsvm( model_log1, model_testing, responseb_test );
    [ model_test_log2, teststats_log2 ] = test_model_logsvm( model_log2, model_testing, responseb_test );
    
    
    %% Support Vector Machines with regularization
    
    [ model_test_svm1, teststats_svm1 ] = test_model_logsvm( model_svm1, model_testing, responseb_test );
    [ model_test_svm2, teststats_svm2 ] = test_model_logsvm( model_svm2, model_testing, responseb_test );
    
    
    %% Naive Bayes Classifier
    
    model_testing_red = model_testing( :, index2 == 1 );
    
    [ model_test_nbayes, teststats_nbayes ] = test_model_nbayes_dectrees( model_nbayes, model_testing_red, responseb_test );
    
    
    %% Decision trees/random forests
    
    [ model_test_treebag, teststats_treebag ] = test_model_nbayes_dectrees( model_treebag, model_testing, responseb_test );
    [ model_test_treeboost, teststats_treeboost ] = test_model_nbayes_dectrees( model_treeboost, model_testing, responseb_test );
    [ model_test_rforest, teststats_rforest ] = test_model_rforest( model_rforest, model_testing, responseb_test );
  
    
    importance_bag = predictorImportance( model_treebag );
    importance_boost = predictorImportance( model_treeboost );
    importance_rforest = model_rforest.OOBPermutedPredictorDeltaError; 
    
    
  
    auc_test = [ auc_test_1, auc_test_2, auc_test_3, auc_test_4, auc_test_5, auc_test_6, auc_nn, auc_test_l, auc_test_en, auc_test_r ];
    
    
    %% 6. Perform ablation studies
  
    models_abl = cell( 12, 1 );
    trainingstats_abl = cell( 12, 1 );
    teststats_abl = cell( 12, 1 );
    
    for l = 1:12
       
        % Naïve Bayes models pose additional requirements on first-step models

        if l == 9 
            
            red_index = index;
            temp = find( index == 1 );
            
            for k = 1:length( index2 )
               
                if index2( 1, k ) == 0
                
                    red_index( 1, temp( k ) ) = 0;
                    
                end
                
            end
            
            red_model = model( :, index2 == 1 );
            red_model_testing = model_testing( :, index2 == 1 );
            
        else
            
            red_index = index;
            red_model = model;
            red_model_testing = model_testing;
            
        end

        
        [ model_abl, trainingstat_abl, teststat_abl ] = ablationstudies( red_model, red_model_testing, responseb_training, responseb_test, red_index, l );
        
        models_abl{ l, 1 } = model_abl;
        trainingstats_abl{ l, 1 } = trainingstat_abl;
        teststats_abl{ l, 1 } = teststat_abl;
        
    end
    
    
    %% 7. Save the results 
    
    models{ j, 1 } = model_nn; 
    models{ j, 2 } = model_l;
    models{ j, 3 } = model_en;
    models{ j, 4 } = model_r;
    models{ j, 5 } = model_log1;
    models{ j, 6 } = model_log2;
    models{ j, 7 } = model_svm1;
    models{ j, 8 } = model_svm2;
    models{ j, 9 } = model_nbayes;
    models{ j, 10 } = model_treebag;
    models{ j, 11 } = model_treeboost;
    models{ j, 12 } = model_rforest;
    
    
    prediction_training{ j, 1 } = model1b_training;
    prediction_training{ j, 2 } = model2b_training;
    prediction_training{ j, 3 } = model3b_training;
    prediction_training{ j, 4 } = model4b_training;
    prediction_training{ j, 5 } = model5b_training;
    prediction_training{ j, 6 } = model6b_training;
    prediction_training{ j, 7 } = modelb_training_nn;
    prediction_training{ j, 8 } = modelb_training_l;
    prediction_training{ j, 9 } = modelb_training_en;
    prediction_training{ j, 10 } = modelb_training_r;
    prediction_training{ j, 11 } = model_training_log1;
    prediction_training{ j, 12 } = model_training_log2;
    prediction_training{ j, 13 } = model_training_svm1;
    prediction_training{ j, 14 } = model_training_svm2;
    prediction_training{ j, 15 } = model_training_nbayes;
    prediction_training{ j, 16 } = model_training_treebag;
    prediction_training{ j, 17 } = model_training_treeboost;
    prediction_training{ j, 18 } = model_training_rforest;
    prediction_training{ j, 19 } = model1_training;
    prediction_training{ j, 20 } = model2_training;
    prediction_training{ j, 21 } = model3_training;
    prediction_training{ j, 22 } = model4_training;
    prediction_training{ j, 23 } = model5_training;
    prediction_training{ j, 24 } = model6_training;
    prediction_training{ j, 25 } = model_training_nn;
    prediction_training{ j, 26 } = model_training_l;
    prediction_training{ j, 27 } = model_training_en;
    prediction_training{ j, 28 } = model_training_r;
    
    
    prediction_test{ j, 1 } = model1b_test;
    prediction_test{ j, 2 } = model2b_test;
    prediction_test{ j, 3 } = model3b_test;
    prediction_test{ j, 4 } = model4b_test;
    prediction_test{ j, 5 } = model5b_test;
    prediction_test{ j, 6 } = model6b_test;
    prediction_test{ j, 7 } = modelb_test_nn;
    prediction_test{ j, 8 } = modelb_test_l;
    prediction_test{ j, 9 } = modelb_test_en;
    prediction_test{ j, 10 } = modelb_test_r;
    prediction_test{ j, 11 } = model_test_log1;
    prediction_test{ j, 12 } = model_test_log2;
    prediction_test{ j, 13 } = model_test_svm1;
    prediction_test{ j, 14 } = model_test_svm2;
    prediction_test{ j, 15 } = model_test_nbayes;
    prediction_test{ j, 16 } = model_test_treebag;
    prediction_test{ j, 17 } = model_test_treeboost;
    prediction_test{ j, 18 } = model_test_rforest;
    prediction_test{ j, 19 } = model1_test;
    prediction_test{ j, 20 } = model2_test;
    prediction_test{ j, 21 } = model3_test;
    prediction_test{ j, 22 } = model4_test;
    prediction_test{ j, 23 } = model5_test;
    prediction_test{ j, 24 } = model6_test;
    prediction_test{ j, 25 } = model_test_nn;
    prediction_test{ j, 26 } = model_test_l;
    prediction_test{ j, 27 } = model_test_en;
    prediction_test{ j, 28 } = model_test_r;
    
    
    response{ j, 1 } = response_training;
    response{ j, 2 } = responseb_training;
    response{ j, 3 } = response_test;
    response{ j, 4 } = responseb_test;
    
    
    trainingstats{ j, 1 } = trainingstats1;
    trainingstats{ j, 2 } = trainingstats2;
    trainingstats{ j, 3 } = trainingstats3;
    trainingstats{ j, 4 } = trainingstats4;
    trainingstats{ j, 5 } = trainingstats5;
    trainingstats{ j, 6 } = trainingstats6;
    trainingstats{ j, 7 } = trainingstats_nn;
    trainingstats{ j, 8 } = trainingstats_l;
    trainingstats{ j, 9 } = trainingstats_en;
    trainingstats{ j, 10 } = trainingstats_r;
    trainingstats{ j, 11 } = trainingstats_log1;
    trainingstats{ j, 12 } = trainingstats_log2;
    trainingstats{ j, 13 } = trainingstats_svm1;
    trainingstats{ j, 14 } = trainingstats_svm2;
    trainingstats{ j, 15 } = trainingstats_nbayes;
    trainingstats{ j, 16 } = trainingstats_treebag;
    trainingstats{ j, 17 } = trainingstats_treeboost;
    trainingstats{ j, 18 } = trainingstats_rforest;
    

    teststats{ j, 1 } = teststats1;
    teststats{ j, 2 } = teststats2;
    teststats{ j, 3 } = teststats3;
    teststats{ j, 4 } = teststats4;
    teststats{ j, 5 } = teststats5;
    teststats{ j, 6 } = teststats6;
    teststats{ j, 7 } = teststats_nn;
    teststats{ j, 8 } = teststats_l;
    teststats{ j, 9 } = teststats_en;
    teststats{ j, 10 } = teststats_r;
    teststats{ j, 11 } = teststats_log1;
    teststats{ j, 12 } = teststats_log2;
    teststats{ j, 13 } = teststats_svm1;
    teststats{ j, 14 } = teststats_svm2;
    teststats{ j, 15 } = teststats_nbayes;
    teststats{ j, 16 } = teststats_treebag;
    teststats{ j, 17 } = teststats_treeboost;
    teststats{ j, 18 } = teststats_rforest;
    
    
    AUCs{ j, 1 } = auc_training;
    AUCs{ j, 2 } = auc_test;
    
    
    var_sets{ j, 1 } = mutset;
    var_sets{ j, 2 } = cnvset;
    var_sets{ j, 3 } = methset;
    var_sets{ j, 4 } = tissueset;
    var_sets{ j, 5 } = mutset_l;
    var_sets{ j, 6 } = cnvset_l;
    var_sets{ j, 7 } = methset_l;

    
    coeffs{ j, 1 } = index;
    coeffs{ j, 2 } = coeffs_l;
    coeffs{ j, 3 } = coeffs_en;
    coeffs{ j, 4 } = coeffs_r;
    coeffs{ j, 5 } = coeffs_log1;
    coeffs{ j, 6 } = coeffs_log2;
    coeffs{ j, 7 } = coeffs_svm1;
    coeffs{ j, 8 } = coeffs_svm2;
    coeffs{ j, 9 } = importance_bag;
    coeffs{ j, 10 } = importance_boost;
    coeffs{ j, 11 } = importance_rforest;
    

    ablstudies_results{ j, 1 } = models_abl;
    ablstudies_results{ j, 2 } = trainingstats_abl;
    ablstudies_results{ j, 3 } = teststats_abl;
    
    
end


end