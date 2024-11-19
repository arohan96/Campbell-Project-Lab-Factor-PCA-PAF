[data, betas, factor_rtns, idiortns] = dataMaker(500,10,4,0.001);
[eigenValsPCA, eigenVecsPCA] = rolling(data, 490, 'PCA', 0.0001, 20, 4);
[eigenValsPAF, eigenVecsPAF] = rolling(data, 490, 'PAF', 0.0001, 20, 4);
eigenValsPCA{1}
eigenValsPAF{1}