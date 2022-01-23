function parsave( fname, models1, models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part )

    save( fname, 'models1', 'models2', 'prediction_training', 'prediction_test', 'response', 'trainingstats', 'teststats', 'AUCs', 'var_sets', 'coeffs', 'ablstudies_results', 'cv_part' );
    
end