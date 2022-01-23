function [ models1, models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part ] = twostepmodel( drug )

% Utilise a cross-validation framework with 10 folds
nfolds = 10;

models1 = cell( nfolds, 9 );
models2 = cell( nfolds, 13 );

prediction_training = cell( nfolds, 37 );
prediction_test = cell( nfolds, 37 );

response = cell( nfolds, 4 );

trainingstats = cell( nfolds, 19 );
teststats = cell( nfolds, 19 );

AUCs = cell( nfolds, 2 );

var_sets = cell( nfolds, 7 );
coeffs = cell( nfolds, 12 );

ablstudies_results = cell( nfolds, 3 );


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

cv_part = cvpartition( n1, 'kFold', nfolds );


for j = 1:nfolds
    
    trainingset = training( cv_part, j );
    testset = test( cv_part, j );
    
    % Utilize the first 7 principal components
    nPCs = 7;

    [ pc_coeff, score, ~ ] = pca( GeneExpression( trainingset, : ) );
    PC_training = score( :, 1:nPCs );
    
    % Project the test gene expression data onto the PC space
    
    PC_test = mrdivide( GeneExpression( testset, : ), pc_coeff' );
    PC_test = PC_test( :, 1:nPCs );
    
    pathway_training = PathwayActivation( trainingset, : );
    pathway_test = PathwayActivation( testset, : );
    
    [ response_training, response_test, responseb_training, responseb_test ] = preprocess_response( Response_AUC( :, drug ), trainingset, testset );
    
    
    %% 2. Compute discrete feature profiles related to shifts in responsiveness
    
    [ clustdata1, mutset, mutset_l, clusterinfo1 ] = mutclusters( SomaticMutation( trainingset, : ), SomaticMutation, response_training );
    
    clustdata1_training = clustdata1( trainingset );
    clustdata1_test = clustdata1( testset );
    n1 = max( clustdata1_training );
    
    
    [ clustdata2, cnvset, cnvset_l, clusterinfo2 ] = cnvclusters( CNV( trainingset, : ), CNV, response_training );
    
    clustdata2_training = clustdata2( trainingset );
    clustdata2_test = clustdata2( testset );
    n2 = max( clustdata2_training );
    
    
    [ clustdata3, methset, methset_l, clusterinfo3 ] = methclusters( Hypermethylation( trainingset, : ), Hypermethylation, response_training );
    
    clustdata3_training = clustdata3( trainingset );
    clustdata3_test = clustdata3( testset );
    n3 = max( clustdata3_training );
    
    
    [ clustdata4, tissueset ] = tissueclusters( TissueType( trainingset, : ), TissueType, response_training );
    
    clustdata4_training = clustdata4( trainingset );
    clustdata4_test = clustdata4( testset );
    n4 = max( clustdata4_training );
    
    
    %% 3. Train and test the six single-data first-step models
    
    [ modelfun1, model1_training, model1b_training, clustmean1, cutoff1, trainingstats1, auc1 ] = mutmodel( responseb_training, clustdata1_training) ;
    [ modelfun2, model2_training, model2b_training, clustmean2, cutoff2, trainingstats2, auc2 ] = CNVmodel( responseb_training, clustdata2_training) ;
    [ modelfun3, model3_training, model3b_training, clustmean3, cutoff3, trainingstats3, auc3 ] = methmodel( responseb_training, clustdata3_training) ;
    [ modelfun4, model4_training, model4b_training, clustmean4, cutoff4, trainingstats4, auc4 ] = tissuemodel( responseb_training, clustdata4_training );
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
    
    
    % Normalize the continuous-valued first-step models to [min(response), max(response)]
    
    [ nmodel, nmodelf5 ] = normsubmodels( model( :, 5 ) );
    model( :, 5 ) = nmodel;
    
    [ nmodel, nmodelf6 ] = normsubmodels( model( :, 6 ) );
    model( :, 6 ) = nmodel;
    
    
    % Store the indices of non-constant first-step models and remove the constant first-step models
    
    index = ones( 1, 6 );
    
    for k = 1:6
        
        if length( unique( model( :, k ) ) ) == 1
            
            index( 1, k ) = 0;
            
        end
        
    end
    
    model = model( :, index == 1 );
    
    
    %% 4.1 Neural Network
    
    % Use Bayesian regularisation backpropagation
    trainFcn = 'trainbr';
    
    % Use five neurons in the hidden layer
    hiddenLayerSize = 5;

    model_nn = fitnet( hiddenLayerSize, trainFcn );
    
    model_nn.divideParam.trainRatio = 80/100;
    model_nn.divideParam.valRatio = 0/100;
    model_nn.divideParam.testRatio = 20/100;
    
    [ model_nn, ~ ] = train( model_nn, model', responseb_training' );
    
    model_training_nn = model_nn( model' );
    
    [ model_training_nn, modelb_training_nn, auc_nn, trainingstats_nn, cutoff_nn, normfunction_nn ] = model_postprocess( model_training_nn', responseb_training );
    
    
    %% 4.2 Regularized linear regression models
    
    % LASSO-regularised linear regression
    [ model_l, fitinfo_l ] = lasso( model, responseb_training );

    % Elastic net-regularised linear regression
    ratio_LASSO_ridge = 0.1;
    [ model_en, fitinfo_en ] = lasso( model, responseb_training, 'Alpha', ratio_LASSO_ridge ); 

    % Approximation of ridge-regularised linear regression
    ratio_LASSO_ridge = 10^(-3);
    [ model_r, fitinfo_r ] = lasso( model, responseb_training, 'Alpha', ratio_LASSO_ridge ); 
    
    
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
    
    % LASSO-regularised logistic regression via a generalised linear model

    % Use 5 cross-validation folds
    nCV = 5;
    % Use 100 regularisation coefficients
    num_regcoeffs = 100;

    [ model_log1, fitinfo_log1 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', num_regcoeffs, 'CV', nCV ); 

    % Elastic net-regularised logistic regression via a generalised linear model
    
    % The LASSO- and the ridge-regularisation term are assigned the same weight
    ratio_LASSO_ridge = 0.5;

    [ model_log2, fitinfo_log2 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', num_regcoeffs, 'CV', nCV, 'Alpha', ratio_LASSO_ridge ); 

    % Ridge-regularised logistic regression via a generalised linear model

    ratio_LASSO_ridge = 10^(-3);

    [ model_log3, fitinfo_log3 ] = lassoglm( model, responseb_training, 'binomial', 'NumLambda', num_regcoeffs, 'CV', nCV, 'Alpha', ratio_LASSO_ridge ); 
    
    
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
    
    
    %% 4.6 Decision Tree Ensembles and Random Forests
    
    % To avoid overfitting, use default decision tree learner templates with only one decision split
    treeStump = templateTree( 'MaxNumSplits', 1 );
    
    model_treebag = fitrensemble( model, responseb_training, 'Method', 'Bag', 'Learners', treeStump ); 
    model_treeboost = fitrensemble( model, responseb_training, 'Method', 'LSBoost','Learners', treeStump ); 
    
    % To avoid overfitting, tree leafs are required to contain at least 20 observations
    num_obs_per_leaf = 20;

    % Compute 50 trees
    num_trees = 50;

    model_rforest = TreeBagger( num_trees, model, responseb_training, 'Method', 'Regression', 'MinLeafSize', num_obs_per_leaf, 'OOBPredictorImportance','on', 'PredictorSelection', 'curvature' );
    
    
    y_treebag = predict( model_treebag, model );
    y_treeboost = predict( model_treeboost, model );
    y_rforest = predict( model_rforest, model );
    
    
    [ model_training_treebag, modelb_training_treebag, auc_treebag, trainingstats_treebag, cutoff_treebag, normfunction_treebag ] = model_postprocess( y_treebag, responseb_training );
    [ model_training_treeboost, modelb_training_treeboost, auc_treeboost, trainingstats_treeboost, cutoff_treeboost, normfunction_treeboost ] = model_postprocess( y_treeboost, responseb_training );
    [ model_training_rforest, modelb_training_rforest, auc_rforest, trainingstats_rforest, cutoff_rforest, normfunction_rforest ] = model_postprocess( y_rforest, responseb_training );
    
    
    auc_training = [ auc1, auc2, auc3, auc4, auc5, auc6, auc_nn, auc_l, auc_en, auc_r, auc_log1, auc_log2, auc_log3, auc_svm1, auc_svm2, auc_treebag, auc_treeboost, auc_rforest ];
    
    
    %% 5. Apply the second-step models to the test set
    
    model_testing = [ model1_test, model2_test, model3_test, model4_test, model5_test, model6_test ];
    model_testing = model_testing( :, index == 1 );
    
    
    % Apply the normalization routines to the continuous first-step models
    
    if index( 5 ) == 1
        
        model_testing( :, sum( index( 1:5 ) ) ) = nmodelf5( model_testing( :, sum( index( 1:5 ) ) ) );
        
    end
    
    
    if index( 6 ) == 1
        
        model_testing( :, sum( index( 1:6 ) ) ) = nmodelf6( model_testing( :, sum( index( 1:6 ) ) ) );
        
    end
    
    
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
    
    
    %% 5.6 Decision Tree Ensembles and Random Forests
    
    [ model_test_treebag, modelb_test_treebag, teststats_treebag, auc_test_treebag ] = test_model_dectrees( model_treebag, model_testing, responseb_test, cutoff_treebag, normfunction_treebag );
    [ model_test_treeboost, modelb_test_treeboost, teststats_treeboost, auc_test_treeboost ] = test_model_dectrees( model_treeboost, model_testing, responseb_test, cutoff_treeboost, normfunction_treeboost );
    [ model_test_rforest, modelb_test_rforest, teststats_rforest, auc_test_rforest ] = test_model_rforest( model_rforest, model_testing, responseb_test, cutoff_rforest, normfunction_rforest );
    
    
    importance_bag = predictorImportance( model_treebag );
    importance_boost = predictorImportance( model_treeboost );
    importance_rforest = model_rforest.OOBPermutedPredictorDeltaError;
    
    
    auc_test = [ auc_test_1, auc_test_2, auc_test_3, auc_test_4, auc_test_5, auc_test_6, auc_test_nn, auc_test_l, auc_test_en, auc_test_r, auc_test_log1, auc_test_log2, auc_test_log3, auc_test_svm1, auc_test_svm2, auc_test_treebag, auc_test_treeboost, auc_test_rforest ];
    
    
    %% 6. Perform ablation studies
    
    models_abl = cell( 13, 1 );
    trainingstats_abl = cell( 13, 1 );
    teststats_abl = cell( 13, 1 );
    
    for l = 1:13
        
        % Naïve Bayes models pose additional requirements on first-step models
        
        if l == 10
            
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
    
    
    %% 7. Sorting the values to be returned
    
    models1{ j, 1 } = clustmean1;
    models1{ j, 2 } = clustmean2;
    models1{ j, 3 } = clustmean3;
    models1{ j, 4 } = clustmean4;
    models1{ j, 5 } = modelfun5;
    models1{ j, 6 } = modelfun6;
    models1{ j, 7 } = clusterinfo1;
    models1{ j, 8 } = clusterinfo2;
    models1{ j, 9 } = clusterinfo3;
    
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
    prediction_training{ j, 11 } = modelb_training_log1;
    prediction_training{ j, 12 } = modelb_training_log2;
    prediction_training{ j, 13 } = modelb_training_log3;
    prediction_training{ j, 14 } = modelb_training_svm1;
    prediction_training{ j, 15 } = modelb_training_svm2;
    prediction_training{ j, 16 } = model_training_nbayes;
    prediction_training{ j, 17 } = modelb_training_treebag;
    prediction_training{ j, 18 } = modelb_training_treeboost;
    prediction_training{ j, 19 } = modelb_training_rforest;
    prediction_training{ j, 20 } = model1_training;
    prediction_training{ j, 21 } = model2_training;
    prediction_training{ j, 22 } = model3_training;
    prediction_training{ j, 23 } = model4_training;
    prediction_training{ j, 24 } = model5_training;
    prediction_training{ j, 25 } = model6_training;
    prediction_training{ j, 26 } = model_training_nn;
    prediction_training{ j, 27 } = model_training_l;
    prediction_training{ j, 28 } = model_training_en;
    prediction_training{ j, 29 } = model_training_r;
    prediction_training{ j, 30 } = model_training_log1;
    prediction_training{ j, 31 } = model_training_log2;
    prediction_training{ j, 32 } = model_training_log3;
    prediction_training{ j, 33 } = model_training_svm1;
    prediction_training{ j, 34 } = model_training_svm2;
    prediction_training{ j, 35 } = model_training_treebag;
    prediction_training{ j, 36 } = model_training_treeboost;
    prediction_training{ j, 37 } = model_training_rforest;
    
    
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
    prediction_test{ j, 11 } = modelb_test_log1;
    prediction_test{ j, 12 } = modelb_test_log2;
    prediction_test{ j, 13 } = modelb_test_log3;
    prediction_test{ j, 14 } = modelb_test_svm1;
    prediction_test{ j, 15 } = modelb_test_svm2;
    prediction_test{ j, 16 } = model_test_nbayes;
    prediction_test{ j, 17 } = modelb_test_treebag;
    prediction_test{ j, 18 } = modelb_test_treeboost;
    prediction_test{ j, 19 } = modelb_test_rforest;
    prediction_test{ j, 20 } = model1_test;
    prediction_test{ j, 21 } = model2_test;
    prediction_test{ j, 22 } = model3_test;
    prediction_test{ j, 23 } = model4_test;
    prediction_test{ j, 24 } = model5_test;
    prediction_test{ j, 25 } = model6_test;
    prediction_test{ j, 26 } = model_test_nn;
    prediction_test{ j, 27 } = model_test_l;
    prediction_test{ j, 28 } = model_test_en;
    prediction_test{ j, 29 } = model_test_r;
    prediction_test{ j, 30 } = model_test_log1;
    prediction_test{ j, 31 } = model_test_log2;
    prediction_test{ j, 32 } = model_test_log3;
    prediction_test{ j, 33 } = model_test_svm1;
    prediction_test{ j, 34 } = model_test_svm2;
    prediction_test{ j, 35 } = model_test_treebag;
    prediction_test{ j, 36 } = model_test_treeboost;
    prediction_test{ j, 37 } = model_test_rforest;
    
    
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
    trainingstats{ j, 13 } = trainingstats_log3;
    trainingstats{ j, 14 } = trainingstats_svm1;
    trainingstats{ j, 15 } = trainingstats_svm2;
    trainingstats{ j, 16 } = trainingstats_nbayes;
    trainingstats{ j, 17 } = trainingstats_treebag;
    trainingstats{ j, 18 } = trainingstats_treeboost;
    trainingstats{ j, 19 } = trainingstats_rforest;
    
    
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
    teststats{ j, 13 } = teststats_log3;
    teststats{ j, 14 } = teststats_svm1;
    teststats{ j, 15 } = teststats_svm2;
    teststats{ j, 16 } = teststats_nbayes;
    teststats{ j, 17 } = teststats_treebag;
    teststats{ j, 18 } = teststats_treeboost;
    teststats{ j, 19 } = teststats_rforest;
    
    
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
    coeffs{ j, 5 } = coefficients_log1;
    coeffs{ j, 6 } = coefficients_log2;
    coeffs{ j, 7 } = coefficients_log3;
    coeffs{ j, 8 } = coeffs_svm1;
    coeffs{ j, 9 } = coeffs_svm2;
    coeffs{ j, 10 } = importance_bag;
    coeffs{ j, 11 } = importance_boost;
    coeffs{ j, 12 } = importance_rforest;
    
    
    ablstudies_results{ j, 1 } = models_abl;
    ablstudies_results{ j, 2 } = trainingstats_abl;
    ablstudies_results{ j, 3 } = teststats_abl;
    
    
end


end