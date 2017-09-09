for K=35:5:50
    [X, class, coefficients, components] = LeeNMF('../digits_train.csv', K, 300); 
    writeCoefficientArray(sprintf('digits_Lee1999_coefficientsAll%d.arff', K), [class coefficients], ',')
    writeComponentArray(sprintf('digits_Lee1999_componentsAll%d.arff', K), components, ',')
end