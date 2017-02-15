au = 'Fairhall SL';
%au = 'Esposito';
%au  = 'Aguirre'
%strfind
%strmatch
%ind = find(cellfun(@(x) strcmp(x, 'Ishai A'),unique_authors))
load('/Users/aidasaglinskas/Desktop/Matlab_Neurosynth/Neurosynth_all.mat')
%% get all the unique authors
t.s  = cellfun(@(x) strsplit(x,','),Neurosynth_all.Database.authors,'UniformOutput',0);
unique_authors = unique([t.s{:}]);
ind = find(cellfun(@(x) strcmp(x,au),unique_authors));
%% Get author vecs
clear auth
tic
disp('Computing author mean vectors')
w_a = [[ind - 10 : ind+10] [2537 26025 24094]];
for i = w_a %1:2%:10%length(unique_authors)
    if ismember(i,w_a(1):round(length(w_a)/10):w_a(end));
    %p = i / length(unique_authors) * 100;
    p = find(w_a == i) / length(w_a) * 100;
    formatSpec = '%4f percent complete';
    %fprintf(formatSpec,p);
    disp(sprintf(formatSpec,p))
    end
a = unique_authors{i};
vecmat(i,:) = mean(Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),unique(Neurosynth_all.Database.id(find(cellfun(@isempty,cellfun(@(x) strfind(x,a),Neurosynth_all.Database.authors,'UniformOutput',0)) == 0))))),2:end),1);
end
toc
% One line method (dope, innit? yeah, too bad it doesnt work :( )
% tic
% mvecs = cellfun(@(y) mean(Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),unique(Neurosynth_all.Database.id(find(cellfun(@isempty,cellfun(@(x) strfind(x,y),Neurosynth_all.Database.authors,'UniformOutput',0)) == 0))))),2:end),1),unique_authors,'UniformOutput',0);
% disp('done')
% toc
%%
plot = 1
if plot
a = [];
a = vecmat(w_a,:);
cmat = corr(a');
lbls = unique_authors(w_a');
f = figure(1)
newVec = get_triu(cmat);
Z = linkage(1-newVec,'ward');
d = subplot(1,2,1);
[h x perm] = dendrogram(Z,'labels',lbls,'orientation','left')
    [h(1:end).LineWidth] = deal(3)
    d.FontSize =  12;
    d.FontWeight = 'bold'
m = subplot(1,2,2)
add_numbers_to_mat(cmat(perm(end:-1:1),perm(end:-1:1)),lbls(perm(end:-1:1)))
m.FontSize =  12
m.FontWeight = 'bold'
m.XTickLabelRotation = 45
end
%% Get author vec
t.a = cellfun(@(x) strfind(x,au),Neurosynth_all.Database.authors,'UniformOutput',0);
t.auth_idns = find(cellfun(@isempty,t.a) == 0);
%Neurosynth_all.Database.authors(t.auth_idns)
auth.idns = t.auth_idns;
auth.studies_ID = unique(Neurosynth_all.Database.id(t.auth_idns));
t.v = Neurosynth_all.features(:,1)
auth.mvec = mean(Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),auth.studies_ID)),2:end))
%auth.vec = Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),auth.studies_ID)),2:end);
%%




