function parsave_one( fname, models2, prediction_training, prediction_test, response, trainingstats, teststats, AUCs, var_sets, coeffs, ablstudies_results, cv_part )

    save( fname, 'models2', 'prediction_training', 'prediction_test', 'response', 'trainingstats', 'teststats', 'AUCs', 'var_sets', 'coeffs', 'ablstudies_results', 'cv_part' );
    
end