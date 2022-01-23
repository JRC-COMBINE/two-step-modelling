function [ models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part ] = onestepmodel( drug )

models2 = cell( 10, 13 );

prediction_training = cell( 10, 25 );
prediction_test = cell( 10, 25 );

response = cell( 10, 4 );

trainingstats = cell( 10, 13 );
teststats = cell( 10, 13 );

AUCs = cell( 10, 2 );

var_sets = cell( 10, 7 );
coeffs = cell( 10, 13 );

ablstudies_results = cell( 10, 3 );

%% 0. Load the input feature data und the response data of the GDSC data base

load( 'SomaticMutation.mat', 'SomaticMutation' );
load( 'CNV.mat', 'CNV' );
load( 'Hypermethylation.mat', 'Hypermethylation' );
load( 'TissueType.mat', 'TissueType' );
load( 'GeneExpression.mat', 'GeneExpression' );
load( 'PathwayActivation.mat', 'PathwayActivation' );
load( 'Response_AUC.mat', 'Response_AUC' );

%% 1. Compute a training and test set partition for cross validation

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
    
    % normalise predictors to [0,1] in order to facilitate the comparison
    % to the two-step models wrt. coeffs
    
    norm1 = nanmin(PC_training); %#ok<*NANMIN> 
    PC_training = PC_training - norm1;
    PC_test = PC_test - norm1;
    
    norm2 = nanmax( PC_training ); %#ok<*NANMAX> 
    PC_training = PC_training ./ norm2;
    PC_test = PC_test ./ norm2;
    
    norm1 = nanmin(pathway_training);
    pathway_training = pathway_training - norm1;
    pathway_test = pathway_test - norm1;
    
    norm2 = nanmax( pathway_training );
    pathway_training = pathway_training ./ norm2;
    pathway_test = pathway_test ./ norm2;
    
    
    %% 2. Compute discrete feature profiles related to shifts in responsiveness
    
    [ mutset, mutset_l ] = mutclusters( SomaticMutation( trainingset, : ), SomaticMutation, response_training );
    mut_features = SomaticMutation( :, mutset );
     
    [ cnvset, cnvset_l ] = cnvclusters( CNV( trainingset, : ), CNV, response_training );
    cnv_features = CNV( :, cnvset );
    
    [ methset, methset_l ] = methclusters( Hypermethylation( trainingset, : ), Hypermethylation, response_training );
    meth_features = Hypermethylation( :, methset );
   
    tissueset = tissueclusters( TissueType( trainingset, : ), response_training );
    tissue_features = TissueType( :, tissueset );
    
    
    %% 3. Set up the model input feature set
    
    index = ones( 1, 6 );
    indexl = zeros( 1, length( mutset ) + length( cnvset ) + length( methset ) + length( tissueset ) + 18 );
    count = 1;
    
    if isempty( mutset )
        
        index( 1, 1 ) = 0;
        
    else
        
        indexl( 1, count:count + length( mutset ) - 1 ) = 1;
        count = count + length( mutset );
        
    end
    
    if isempty( cnvset )
        
        index( 1, 2 ) = 0;
        
    else
        
        indexl( 1, count:count + length( cnvset ) - 1 ) = 2;
        count = count + length( cnvset );
        
    end
    
    if isempty( methset )
        
        index( 1, 3 ) = 0;
        
    else
        
        indexl( 1, count:count + length( methset ) - 1 ) = 3;
        count = count + length( methset );
        
    end
    
    if isempty( tissueset )
        
        index( 1, 4 ) = 0;
        
    else
        
        indexl( 1, count:count + length( tissueset ) - 1 ) = 4;
        count = count + length( tissueset );
        
    end
    
    indexl( 1, count:count + 10 ) = 5;
    indexl( 1, count + 11:length( indexl ) ) = 6;
    
    model = [ mut_features( trainingset, : ), cnv_features( trainingset, : ), meth_features( trainingset, : ), tissue_features( trainingset, : ), pathway_training, PC_training ];

    %% 4. Train the models
    
    %% 4.1 Neural Network
    
    trainFcn = 'trainbr';
    
    hiddenLayerSize = 5;
    model_nn = fitnet( hiddenLayerSize, trainFcn );
    
    model_nn.divideParam.trainRatio = 80/100;
    model_nn.divideParam.valRatio = 0/100;
    model_nn.divideParam.testRatio = 20/100;
    
    [ model_nn, ~ ] = train( model_nn, model', responseb_training' );
    
    model_training_nn = model_nn( model' );
    
    [ model_training_nn, modelb_training_nn, auc_nn, trainingstats_nn, cutoff_nn, normfunction_nn ] = model_postprocess( model_training_nn', responseb_training );
    
    
    %% 4.2 Regularized linear regression models
    
    [ model_l, fitinfo_l ] = lasso( model, responseb_training );
    [ model_en, fitinfo_en ] = lasso( model, responseb_training, 'Alpha', 0.1 ); % elastic net regularization
    [ model_r, fitinfo_r ] = lasso( model, responseb_training, 'Alpha', 10^(-3) ); % approximation of ridge regularization
    
    
    [ ~, ind ] = min( fitinfo_l.MSE );
    
    coef = model_l( :, ind );
    constant = fitinfo_l.Intercept( ind );
    y_l = model*coef + constant;
    
    coeffs_l = [ constant; coef ];
    
    
    [ ~, ind ] = min( fitinfo_en.MSE );
    
    coef = model_en( :, ind );
    constant = fitinfo_en.Intercept( ind );
    y_en = model*coef + constant;
    
    coeffs_en = [ constant; coef ];
    
    
    [ ~, ind ] = min( fitinfo_r.MSE );
    
    coef = model_r( :, ind );
    constant = fitinfo_r.Intercept( ind );
    y_r = model*coef + constant;
    
    coeffs_r = [ constant; coef ];
    
    
    [ model_training_l, modelb_training_l, auc_l, trainingstats_l, cutoff_l, normfunction_l ] = model_postprocess( y_l, responseb_training );
    [ model_training_en, modelb_training_en, auc_en, trainingstats_en, cutoff_en, normfunction_en ] = model_postprocess( y_en, responseb_training );
    [ model_training_r, modelb_training_r, auc_r, trainingstats_r, cutoff_r, normfunction_r ] = model_postprocess( y_r, responseb_training );
    
    
    %% 4.3 Regularized logistic regression models
    
    [ model_log1, fitinfo_log1 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5 ); % lasso
    [ model_log2, fitinfo_log2 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5, 'Alpha', 0.5 ); % EN
    [ model_log3, fitinfo_log3 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', 100, 'CV', 5, 'Alpha', 10^(-3) ); % ridge
    
    
    index_log1 = fitinfo_log1.Index1SE;
    index_log2 = fitinfo_log2.Index1SE;
    index_log3 = fitinfo_log3.Index1SE;
    
    
    constant_log1 = fitinfo_log1.Intercept( index_log1 );
    coeffs_log1 = model_log1( :, index_log1 );
    coefficients_log1 = [ constant_log1; coeffs_log1 ];
    
    
    constant_log2 = fitinfo_log2.Intercept( index_log2 );
    coeffs_log2 = model_log2( :, index_log2 );
    coefficients_log2 = [ constant_log2; coeffs_log2 ];
    
    
    constant_log3 = fitinfo_log3.Intercept( index_log3 );
    coeffs_log3 = model_log3( :, index_log3 );
    coefficients_log3 = [ constant_log3; coeffs_log3 ];
    
    
    y_log1 = glmval( coefficients_log1, model, 'logit' );
    y_log2 = glmval( coefficients_log2, model, 'logit' );
    y_log3 = glmval( coefficients_log3, model, 'logit' );
    
    
    [ model_training_log1, modelb_training_log1, auc_log1, trainingstats_log1, cutoff_log1, normfunction_log1 ] = model_postprocess( y_log1, responseb_training );
    [ model_training_log2, modelb_training_log2, auc_log2, trainingstats_log2, cutoff_log2, normfunction_log2 ] = model_postprocess( y_log2, responseb_training );
    [ model_training_log3, modelb_training_log3, auc_log3, trainingstats_log3, cutoff_log3, normfunction_log3 ] = model_postprocess( y_log3, responseb_training );
    
    
    %% 4.4 Regularized Support Vector Machines
    
    [ model_svm1, ~ ] = fitrlinear( model, responseb_training, 'Learner', 'svm', 'Regularization', 'ridge' );
    [ model_svm2, ~ ] = fitrlinear( model, responseb_training, 'Learner', 'svm', 'Regularization', 'lasso' );
    
    
    y_svm1 = predict( model_svm1, model );
    y_svm2 = predict( model_svm2, model );
    
    coeffs_svm1 = [ model_svm1.Bias, model_svm1.Beta' ];
    coeffs_svm2 = [ model_svm2.Bias, model_svm2.Beta' ];
    
    
    [ model_training_svm1, modelb_training_svm1, auc_svm1, trainingstats_svm1, cutoff_svm1, normfunction_svm1 ] = model_postprocess( y_svm1, responseb_training );
    [ model_training_svm2, modelb_training_svm2, auc_svm2, trainingstats_svm2, cutoff_svm2, normfunction_svm2 ] = model_postprocess( y_svm2, responseb_training );
    
    
    %% 4.5 Naïve Bayes Classifier
    
    index2 = ones( 1, size( model, 2 ) );
    
    % Remove features lacking variation in their prediction for sensitive or resistant cell lines
    
    for k = 1:size( model, 2 )
        
        if length( unique( model( responseb_training == 0, k ) ) ) == 1 || length( unique( model( responseb_training == 1, k ) ) ) == 1
            
            index2( 1, k ) = 0;
            
        end
        
    end
    
    model_red = model( :, index2 == 1 );
    
    
    model_nbayes = fitcnb( model_red, responseb_training );
    y_nbayes = predict( model_nbayes, model_red );
    
    [ model_training_nbayes, ~, ~, trainingstats_nbayes, ~, ~ ] = model_postprocess( y_nbayes, responseb_training );
    
    
    %% 4.6 Decision Tree Ensembles and Random Forests
    
    treeStump = templateTree( 'MaxNumSplits', 1 );
    
    model_treebag = fitrensemble( model, responseb_training, 'Method', 'Bag', 'Learners', treeStump );
    model_treeboost = fitrensemble( model, responseb_training, 'Method', 'LSBoost','Learners', treeStump );
    model_rforest = TreeBagger( 50, model, responseb_training, 'Method', 'Regression', 'MinLeafSize', 20, 'OOBPredictorImportance','on', 'PredictorSelection', 'curvature' );
    
    
    y_treebag = predict( model_treebag, model );
    y_treeboost = predict( model_treeboost, model );
    y_rforest = predict( model_rforest, model );
    
    
    [ model_training_treebag, modelb_training_treebag, auc_treebag, trainingstats_treebag, cutoff_treebag, normfunction_treebag ] = model_postprocess( y_treebag, responseb_training );
    [ model_training_treeboost, modelb_training_treeboost, auc_treeboost, trainingstats_treeboost, cutoff_treeboost, normfunction_treeboost ] = model_postprocess( y_treeboost, responseb_training );
    [ model_training_rforest, modelb_training_rforest, auc_rforest, trainingstats_rforest, cutoff_rforest, normfunction_rforest ] = model_postprocess( y_rforest, responseb_training );
    
    
    auc_training = [ auc_nn, auc_l, auc_en, auc_r, auc_log1, auc_log2, auc_log3, auc_svm1, auc_svm2, auc_treebag, auc_treeboost, auc_rforest ];
    
    
    %% 5. Apply the second-step models to the test set
    
    model_testing = [ mut_features( testset, : ), cnv_features( testset, : ), meth_features( testset, : ), tissue_features( testset, : ), pathway_test, PC_test ];
    
    %% 5.1 Neural Network
    
    [ model_test_nn, modelb_test_nn, teststats_nn, auc_test_nn ] = test_model_NN( model_nn, model_testing, cutoff_nn, responseb_test, normfunction_nn );
    
    
    %% 5.2 Regularized linear regression models
    
    [ model_test_l, modelb_test_l, teststats_l, auc_test_l ] = test_model_linreg( coeffs_l, model_testing, cutoff_l, responseb_test, normfunction_l );
    [ model_test_en, modelb_test_en, teststats_en, auc_test_en ] = test_model_linreg( coeffs_en, model_testing, cutoff_en, responseb_test, normfunction_en );
    [ model_test_r, modelb_test_r, teststats_r, auc_test_r ] = test_model_linreg( coeffs_r, model_testing, cutoff_r, responseb_test, normfunction_r );
    
       
    %% 5.3 Regularized logistic regression models
    
    [ model_test_log1, modelb_test_log1, teststats_log1, auc_test_log1 ] = test_model_log( coefficients_log1, model_testing, responseb_test, cutoff_log1, normfunction_log1 );
    [ model_test_log2, modelb_test_log2, teststats_log2, auc_test_log2 ] = test_model_log( coefficients_log2, model_testing, responseb_test, cutoff_log2, normfunction_log2 );
    [ model_test_log3, modelb_test_log3, teststats_log3, auc_test_log3 ] = test_model_log( coefficients_log3, model_testing, responseb_test, cutoff_log3, normfunction_log3 );

    
    %% 5.4 Regularized Support Vector Machines
    
    [ model_test_svm1, modelb_test_svm1, teststats_svm1, auc_test_svm1 ] = test_model_svm( model_svm1, model_testing, responseb_test, cutoff_svm1, normfunction_svm1 );
    [ model_test_svm2, modelb_test_svm2, teststats_svm2, auc_test_svm2 ] = test_model_svm( model_svm2, model_testing, responseb_test, cutoff_svm2, normfunction_svm2 );
    
    
    %% 5.5 Naïve Bayes Classifier
    
    model_testing_red = model_testing( :, index2 == 1 );
    
    [ model_test_nbayes, teststats_nbayes ] = test_model_nbayes( model_nbayes, model_testing_red, responseb_test );
    
    
    %% 5.6 Decision Tree Ensemblesand Random Forests
    
    [ model_test_treebag, modelb_test_treebag, teststats_treebag, auc_test_treebag ] = test_model_dectrees( model_treebag, model_testing, responseb_test, cutoff_treebag, normfunction_treebag );
    [ model_test_treeboost, modelb_test_treeboost, teststats_treeboost, auc_test_treeboost ] = test_model_dectrees( model_treeboost, model_testing, responseb_test, cutoff_treeboost, normfunction_treeboost );
    [ model_test_rforest, modelb_test_rforest, teststats_rforest, auc_test_rforest ] = test_model_rforest( model_rforest, model_testing, responseb_test, cutoff_rforest, normfunction_rforest );
    
    
    importance_bag = predictorImportance( model_treebag );
    importance_boost = predictorImportance( model_treeboost );
    importance_rforest = model_rforest.OOBPermutedPredictorDeltaError;
    
    
    auc_test = [ auc_test_nn, auc_test_l, auc_test_en, auc_test_r, auc_test_log1, auc_test_log2, auc_test_log3, auc_test_svm1, auc_test_svm2, auc_test_treebag, auc_test_treeboost, auc_test_rforest ];
  
    %% 6. Perform ablation studies
    
    models_abl = cell( 13, 1 );
    trainingstats_abl = cell( 13, 1 );
    teststats_abl = cell( 13, 1 );
    
    for l = 1:13
        
        % Naïve Bayes models pose additional requirements on features
        
        if l == 10
            
            red_indexl = indexl( :, index2 == 1 );
            
            red_index = zeros( 1, 6 );
            red_index( unique( red_indexl ) ) = 1; 
            
            red_model = model( :, index2 == 1 );
            red_model_testing = model_testing( :, index2 == 1 );
            
        else
            
            red_index = index;
            red_indexl = indexl;
            red_model = model;
            red_model_testing = model_testing;
            
        end
        
        
        [ model_abl, trainingstat_abl, teststat_abl ] = ablationstudies( red_model, red_model_testing, responseb_training, responseb_test, red_index, red_indexl, l );
        
        models_abl{ l, 1 } = model_abl;
        trainingstats_abl{ l, 1 } = trainingstat_abl;
        teststats_abl{ l, 1 } = teststat_abl;
        
    end

    
    %% Sorting the values to be returned
    
    models2{ j, 1 } = model_nn;
    models2{ j, 2 } = model_l;
    models2{ j, 3 } = model_en;
    models2{ j, 4 } = model_r;
    models2{ j, 5 } = model_log1;
    models2{ j, 6 } = model_log2;
    models2{ j, 7 } = model_log3;
    models2{ j, 8 } = model_svm1;
    models2{ j, 9 } = model_svm2;
    models2{ j, 10 } = model_nbayes;
    models2{ j, 11 } = model_treebag;
    models2{ j, 12 } = model_treeboost;
    models2{ j, 13 } = model_rforest;
    
    
    prediction_training{ j, 1 } = modelb_training_nn;
    prediction_training{ j, 2 } = modelb_training_l;
    prediction_training{ j, 3 } = modelb_training_en;
    prediction_training{ j, 4 } = modelb_training_r;
    prediction_training{ j, 5 } = modelb_training_log1;
    prediction_training{ j, 6 } = modelb_training_log2;
    prediction_training{ j, 7 } = modelb_training_log3;
    prediction_training{ j, 8 } = modelb_training_svm1;
    prediction_training{ j, 9 } = modelb_training_svm2;
    prediction_training{ j, 10 } = model_training_nbayes;
    prediction_training{ j, 11 } = modelb_training_treebag;
    prediction_training{ j, 12 } = modelb_training_treeboost;
    prediction_training{ j, 13 } = modelb_training_rforest;
    prediction_training{ j, 14 } = model_training_nn;
    prediction_training{ j, 15 } = model_training_l;
    prediction_training{ j, 16 } = model_training_en;
    prediction_training{ j, 17 } = model_training_r;
    prediction_training{ j, 18 } = model_training_log1;
    prediction_training{ j, 19 } = model_training_log2;
    prediction_training{ j, 20 } = model_training_log3;
    prediction_training{ j, 21 } = model_training_svm1;
    prediction_training{ j, 22 } = model_training_svm2;
    prediction_training{ j, 23 } = model_training_treebag;
    prediction_training{ j, 24 } = model_training_treeboost;
    prediction_training{ j, 25 } = model_training_rforest;
    
    
    prediction_test{ j, 1 } = modelb_test_nn;
    prediction_test{ j, 2 } = modelb_test_l;
    prediction_test{ j, 3 } = modelb_test_en;
    prediction_test{ j, 4 } = modelb_test_r;
    prediction_test{ j, 5 } = modelb_test_log1;
    prediction_test{ j, 6 } = modelb_test_log2;
    prediction_test{ j, 7 } = modelb_test_log3;
    prediction_test{ j, 8 } = modelb_test_svm1;
    prediction_test{ j, 9 } = modelb_test_svm2;
    prediction_test{ j, 10 } = model_test_nbayes;
    prediction_test{ j, 11 } = modelb_test_treebag;
    prediction_test{ j, 12 } = modelb_test_treeboost;
    prediction_test{ j, 13 } = modelb_test_rforest;
    prediction_test{ j, 14 } = model_test_nn;
    prediction_test{ j, 15 } = model_test_l;
    prediction_test{ j, 16 } = model_test_en;
    prediction_test{ j, 17 } = model_test_r;
    prediction_test{ j, 18 } = model_test_log1;
    prediction_test{ j, 19 } = model_test_log2;
    prediction_test{ j, 20 } = model_test_log3;
    prediction_test{ j, 21 } = model_test_svm1;
    prediction_test{ j, 22 } = model_test_svm2;
    prediction_test{ j, 23 } = model_test_treebag;
    prediction_test{ j, 24 } = model_test_treeboost;
    prediction_test{ j, 25 } = model_test_rforest;
    
    
    response{ j, 1 } = response_training;
    response{ j, 2 } = responseb_training;
    response{ j, 3 } = response_test;
    response{ j, 4 } = responseb_test;

    
    trainingstats{ j, 1 } = trainingstats_nn;
    trainingstats{ j, 2 } = trainingstats_l;
    trainingstats{ j, 3 } = trainingstats_en;
    trainingstats{ j, 4 } = trainingstats_r;
    trainingstats{ j, 5 } = trainingstats_log1;
    trainingstats{ j, 6 } = trainingstats_log2;
    trainingstats{ j, 7 } = trainingstats_log3;
    trainingstats{ j, 8 } = trainingstats_svm1;
    trainingstats{ j, 9 } = trainingstats_svm2;
    trainingstats{ j, 10 } = trainingstats_nbayes;
    trainingstats{ j, 11 } = trainingstats_treebag;
    trainingstats{ j, 12 } = trainingstats_treeboost;
    trainingstats{ j, 13 } = trainingstats_rforest;
    
    
    teststats{ j, 1 } = teststats_nn;
    teststats{ j, 2 } = teststats_l;
    teststats{ j, 3 } = teststats_en;
    teststats{ j, 4 } = teststats_r;
    teststats{ j, 5 } = teststats_log1;
    teststats{ j, 6 } = teststats_log2;
    teststats{ j, 7 } = teststats_log3;
    teststats{ j, 8 } = teststats_svm1;
    teststats{ j, 9 } = teststats_svm2;
    teststats{ j, 10 } = teststats_nbayes;
    teststats{ j, 11 } = teststats_treebag;
    teststats{ j, 12 } = teststats_treeboost;
    teststats{ j, 13 } = teststats_rforest;
    
    
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
    coeffs{ j, 2 } = indexl;
    coeffs{ j, 3 } = coeffs_l;
    coeffs{ j, 4 } = coeffs_en;
    coeffs{ j, 5 } = coeffs_r;
    coeffs{ j, 6 } = coefficients_log1;
    coeffs{ j, 7 } = coefficients_log2;
    coeffs{ j, 8 } = coefficients_log3;
    coeffs{ j, 9 } = coeffs_svm1;
    coeffs{ j, 10 } = coeffs_svm2;
    coeffs{ j, 11 } = importance_bag;
    coeffs{ j, 12 } = importance_boost;
    coeffs{ j, 13 } = importance_rforest;
    
    
    ablstudies_results{ j, 1 } = models_abl;
    ablstudies_results{ j, 2 } = trainingstats_abl;
    ablstudies_results{ j, 3 } = teststats_abl;
    
end

end