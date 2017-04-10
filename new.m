load('/Users/aidasaglinskas/Desktop/Matlab_Neurosynth/Neurosynth_all.mat')
n = Neurosynth_all;
%%
all_include_w_inds = [];
include = {'face' 'faces' 'facial' 'person' 'people'};
include_drop = {'surface'};

exclude = {'social'};
exclude_drop = {};
% get all word inds for 'include' list
for w_ind = 1:length(include)
f  = cellfun(@(x) strfind(x,include{w_ind}),n.labels,'UniformOutput',0);
ff = cellfun(@isempty,f) == 0;
ff_inds = find(ff);
all_include_w_inds = [all_include_w_inds;ff_inds];
end

% include drop
all_include_drop_w_inds = [];
for w_ind = 1:length(include_drop);
f  = cellfun(@(x) strfind(x,include_drop{w_ind}),n.labels,'UniformOutput',0);
ff = cellfun(@isempty,f) == 0;
ff_inds = find(ff);
all_include_drop_w_inds = [all_include_drop_w_inds;ff_inds];
end
all_include_w_inds(all_include_w_inds == ff_inds) = [];

    inc.inds = all_include_w_inds;
    inc.words = {n.labels{all_include_w_inds}}';

    
    
% Exclude drop
all_exclude_w_inds = [];
for w_ind = 1:length(exclude)
f  = cellfun(@(x) strfind(x,exclude{w_ind}),n.labels,'UniformOutput',0);
ff = cellfun(@isempty,f) == 0;
ff_inds = find(ff);
all_exclude_w_inds = [all_exclude_w_inds;ff_inds];
end
n.labels(all_exclude_w_inds)

all_exclude_drop_w_inds = [];
for w_ind = 1:length(exclude_drop);
f  = cellfun(@(x) strfind(x,exclude_drop{w_ind}),n.labels,'UniformOutput',0);
ff = cellfun(@isempty,f) == 0;
ff_inds = find(ff);
all_exclude_drop_w_inds = [all_exclude_drop_w_inds;ff_inds];
end
all_include_w_inds(all_include_w_inds == ff_inds) = [];



%% Find co-occuring terms
r = find(sum(n.features(:,inc.inds),2));
c = ismember([1:size(n.labels,1)],[inc.inds;1]) == 0;
t.inds = sum(n.features(r,c),1);
t.lbls = n.labels(c);
[Y,I] = sort(t.inds','descend');
sort_cooccuring.inds = t.inds(I)';
sort_cooccuring.words = t.lbls(I)';
%tt = t.lbls(I)
%tt(1:10)
i = 1:10;
ct = table(sort_cooccuring.inds(i),sort_cooccuring.words(i)','VariableNames',{'Magnitude' 'Term'});
disp(sprintf('Top %d Co-occurint terms',max(i)))
disp(ct)
%%

inc.inds






