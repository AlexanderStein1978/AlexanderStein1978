load fdData;

[U, S] = pca(data);
K = 100;
Z = projectData(data, U, K);


