load('/Users/aidasaglinskas/Desktop/Matlab_Neurosynth/Neurosynth_all.mat');
n = Neurosynth_all;
%% Individual Authors;
tic
ua.lines = unique(n.Database.authors,'stable');
ua.studyID = unique(n.Database.id,'stable');
ua.linesSplit = cellfun(@(x) strsplit(x,', '),ua.lines,'UniformOutput',0);
ua.individual_authors = [ua.linesSplit{:}]';
ua.unique_authors = unique(ua.individual_authors,'stable')
toc
%%
clc
auth_mat = [];
tic
for auth_ind = 1:length(ua.unique_authors);
    repvec = [1:round(length(ua.unique_authors)/100):length(ua.unique_authors)];
    if ismember(auth_ind,repvec)
       disp([num2str(auth_ind / length(ua.unique_authors) * 100) '% done']);
       toc
    end
st.auth_inds = arrayfun(@(x) ismember(ua.unique_authors{auth_ind},ua.linesSplit{x}),[1:length(ua.linesSplit)]);
auth_mat.mat(auth_ind,:) = mean(n.features(st.auth_inds,:),1);
auth_mat.auth{auth_ind} = ua.unique_authors{auth_ind};
end
%% Save
disp('saving')
tic
save('/Users/aidasaglinskas/Desktop/auth_mats_wrkspc')
save('/Users/aidasaglinskas/Desktop/auth_mat','auth_mat')
toc
%% Compute similarity
tic
disp('loading mat')
load('/Users/aidasaglinskas/Desktop/auth_mat.mat')
%
toc
disp('Computing Correlations')

% 10029 % FAIRHALL
% 3191 Peelem%
% 5456 Cross ES
% 13121 Caramzza 
% 2908 Piazza
% 3191 Ishai

rng(1)
%auth_mat.mat = zscore(auth_mat.mat,[],1);

CIMeC = [10029 3191 5456 13121 2908 672];
other = [18972       10316        2213        6492        3744        6846 10721       12842       11143       21325]
author_inds = [CIMeC randi(length(auth_mat.auth),1,1000)];  
a = find(auth_mat.mat(10029,:))
feature_inds = 2:size(auth_mat.mat,2); %a(2:end)%
auth_mat.corr = corr(auth_mat.mat(author_inds,feature_inds)','type','Spearman');

%find(ismember(auth_mat.auth,'Meyer M'))
tic
disp('Computing Triu')
auth_mat.triu = get_triu(auth_mat.corr);
toc
dend_fig = figure(1)
dend_fig.Visible = 'on'
disp('Computing Linkage')
auth_mat.Z = linkage(1-auth_mat.triu,'ward');

%linkage(1-auth_mat.corr,'complete')

[Y I] = sort(auth_mat.mat(1382,2:end),'descend');
n.labels(I(1:20)+1)
toc
disp('Computing Dendrogram')
[h x auth_mat.perm] = dendrogram(auth_mat.Z,0,'orientation','left','labels',auth_mat.auth(author_inds))
toc
auth_mat.ord = auth_mat.perm(end:-1:1);
disp('all done')

figure(2)
imagesc(auth_mat.corr(auth_mat.ord,auth_mat.ord))


%%

a1 = auth_mat.mat(10029,:);
[Y I] = sort(a1,'descend');
Ys = n.labels(I);

a2 = auth_mat.mat(3191,:);
[Y I] = sort(a2,'descend');
Yc = n.labels(I);

% 10029 % FAIRHALL
% 3191 Peelem%
% 5456 Cross ES
% 13121 Caramzza 
% 2908 Piazza
% 3191 Ishai


% Find minimul communality
m = [];
for i = 1:length(Ys)
m(i,1) = i;

if numel(find(strcmp(Yc,Ys{i}))) > 1 
    disp(Yc(find(strcmp(Yc,Ys{i}))))
m(i,2)= max(find(strcmp(Yc,Ys{i})));
else
    m(i,2)= find(strcmp(Yc,Ys{i}));
end

end
ind =sum(m,2);
ind(1) = nan
Yc(find(ind == min(ind)))