[data, betas, factor_rtns, idiortns] = dataMaker(500,10,4,0.001);
[eigenValsPCA, eigenVecsPCA] = baseRolling(data, 490, 'PCA', 0.0001, 20);
[eigenValsPAF, eigenVecsPAF] = baseRolling(data, 490, 'PAF', 0.0001, 20);
eigenVecsPCA{1}
eigenVecsPAF{1}