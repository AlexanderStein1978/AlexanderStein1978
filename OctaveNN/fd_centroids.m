dirname = uigetdir();
rawFileList = readdir(dirname);
%celldisp(fileList);
count=0;
for i=1:size(rawFileList, 1)
  if 0<size(strfind(rawFileList{i}, "png"), 1)
    count=count+1;
  endif
end
fileList = cell(count);
count=1;
for i=1:size(rawFileList, 1)
  if 0<size(strfind(rawFileList{i}, "png"), 1)
    fileList{count}=rawFileList{i};
    count=count+1;
  endif
end
classVec = strncmp ("No", fileList, 2);
%fprintf('%d', classVec);
number = size(classVec, 1);
K = 16;
max_iters = 10;
data = zeros(number, 10000);
Gcontroids = zeros(number, K, 3);
for i=1:number
  fprintf("%d of %d, %s\n", i, number, fileList{i});
  A = double(imread([dirname '\\' fileList{i}]));
  A = A / 255;
  img_size = size(A);
  X = reshape(A, img_size(1) * img_size(2), 3);
  initial_centroids = kMeansInitCentroids(X, K);
  [centroids, idx] = runkMeans(X, initial_centroids, max_iters);
  data(i, :) = findClosestCentroids(X, centroids);
  Gcontroids(i,:,:) = centroids;
end

[U, S] = pca(data);
K = 100;
Z = projectData(data, U, K);

fprintf("Calculations finished.");

save -binary classVecF classVec
fprintf("Saved classVec.");
save -binary centroidsF Gcontroids
fprintf("Saved Gcontroids.");
save -binary dataF data
fprintf("Saved data.");
save -binary UF U
fprintf("Saved U.");
save -binary SF S
fprintf("Saved S.");
save -binary ZF Z
fprintf("Saved Z.");


