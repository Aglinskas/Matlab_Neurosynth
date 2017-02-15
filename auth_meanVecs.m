clear all
load('/Users/aidasaglinskas/Desktop/Matlab_Neurosynth/Neurosynth_all.mat')
% get all the unique authors
disp('Looking up unique authors')
t.s  = cellfun(@(x) strsplit(x,','),Neurosynth_all.Database.authors,'UniformOutput',0);
unique_authors = unique([t.s{:}]);
% Get author vecs
clear auth
tic
disp('Computing author mean vectors')
for i = 1:length(unique_authors)%:10%length(unique_authors)
    if ismember(10,1:10:length(unique_authors));
    p = i / length(unique_authors) * 100;
    disp(sprintf('%4f percent complete',p));
    end
a = unique_authors{i};
vecmat(i,:) = mean(Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),unique(Neurosynth_all.Database.id(find(cellfun(@isempty,cellfun(@(x) strfind(x,a),Neurosynth_all.Database.authors,'UniformOutput',0)) == 0))))),2:end),1);
end
toc
save('/Users/aidasaglinskas/Desktop/Matlab_Neurosynth/vecmat.mat')
disp('all done')