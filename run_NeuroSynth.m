%CMD+ENTER to run portion of the script 
%% Part 0: Load up datasets and specify parameters (RUN ONCE)
clear all
% Load Thangs
p.root_folder = '/Users/aidasaglinskas/Desktop/Neurosynth/';

disp('loading database & feature variables')
load(fullfile(p.root_folder,'Neurosynth_all.mat'))
load(fullfile(p.root_folder,'NeuroSynth_labels.mat'))
%load('/Users/aidasaglinskas/Desktop/Neurosynth/face_inds.mat')
disp('loaded')
% Parameters
p.name_of_MRI_template = 'single_subj_T1.nii';
p.meta_analysis_output_fn = 'meta.nii';
p.meta_analysis_output_fn_smoothed = 'smeta.nii';
k = 3; % smoothing kernel for meta analysis;
%% Part 1; Choose word list (RUN UNTIL SATISFIED)
% choose words to include and exclude until satisfied, proceed

words = {'face' 'facial' 'person' 'people'};
to_drop = {'al' 'surface' 'rs' 'interpersonal'};

% % % % %
clear t tt list inds
l = cellfun(@num2str,Neurosynth_all.labels,'UniformOutput',0); %list of all words
for w = 1:length(words);
t.o = arrayfun(@(x) findstr(l{x},words{w}),1:length(l),'UniformOutput',0);
tt{w} = cellfun(@isempty,t.o) == 0;
inds{w} = find(tt{w});
end
list.words = l([inds{:}]);
list.inds = [inds{:}];
% Words to drop
if isempty(to_drop) == 0
clear w_d d_temp d
for w_d = 1:length(to_drop)
d_temp = cellfun(@(x) strmatch(x,to_drop{w_d}),l([inds{:}]),'UniformOutput',0);
d{w_d} = find(cellfun(@isempty,d_temp) == 0);
end
disp('Dropping these from list');
disp(list.words([d{:}]));
list.words([d{:}]) = [];
list.inds([d{:}]) = [];
end
disp('FINAL LIST')
disp(list.words)
% save('/Users/aidasaglinskas/Desktop/Neurosynth/face_inds.mat','list')
%% Part 2: Get the coordinates, map to brain, smooth etc (RUN ONCE)
temp_fn = fullfile(p.root_folder,p.name_of_MRI_template);
ofn = fullfile(p.root_folder,p.meta_analysis_output_fn);
ofsn = fullfile(p.root_folder,p.meta_analysis_output_fn_smoothed);

n_reduced = Neurosynth_all.features(:,[1 list.inds]);
% Find studies containing the keywords;
temp.p = sum(n_reduced(:,2:size(n_reduced,2)),2);
studies_matching_ind = find(temp.p);
studies_matching_IDs = n_reduced(studies_matching_ind,1);
disp(sprintf('%d Studies found matching criteria',length(studies_matching_IDs)))
% Find indeces, coords, name of journal etc, matching query
%t.m = ismember(Neurosynth_all.Database.id,studies_matching_IDs);
db_inds = find(ismember(Neurosynth_all.Database.id,studies_matching_IDs));

% Get Rid of non MNI coords; 
t.p1 = cellfun(@(x) strmatch(x,'MNI'),Neurosynth_all.Database.space(db_inds),'UniformOutput',0);
t.p2 = find(cellfun(@isempty,t.p1) == 1);
disp(sprintf('%d studies dropped due to non-MNI coordinates',length(unique(Neurosynth_all.Database.title(db_inds(t.p2))))));
db_inds(t.p2) = [];
disp(sprintf('%d studies included in meta-analysis',length(unique(Neurosynth_all.Database.title(db_inds)))));
meta.db_inds = db_inds;
meta.x = Neurosynth_all.Database.x(meta.db_inds);
meta.y = Neurosynth_all.Database.y(meta.db_inds);
meta.z = Neurosynth_all.Database.z(meta.db_inds);
% map to a brain
disp('Mapping to brain')
im = load_nii(temp_fn); %load
im.img(:) = deal(0); %wipe
% Transformation Matrix from mm to vx
m = [im.hdr.hist.srow_x;im.hdr.hist.srow_y;im.hdr.hist.srow_z;[0 0 0 1]];
m(:,4) = m(:,4) - sum(m(:,1:3),2); % adjustment fuckery
try
a = cosmo_fmri_dataset(temp_fn);
%fprintf('looks like you have cosmo installed..\n..lemme run this additional check real quick')
mat = a.a.vol.mat;
if all(m(:) == mat(:));
%    fprintf('\ntransformation matrix looks good nevermind')
else warning('hmm.. transformation matrix might be fucked, check if the two below are similar')
disp('Transformation matrix from nifti header')
disp(m)
disp('CosmoMVPA''s transformation matrix')
disp(mat)
disp('Fuck it, I''ll use cosmo''s instead')
m = mat;
end
catch
end
%
disp('Mapping Coordinates')
shitty_cords = [];
for i = 1:length(meta.db_inds)
c_mm = [meta.x(i) meta.y(i) meta.z(i)];
c_vx = inv(mat)*[c_mm 1]';
c_vx = round(c_vx(1:3));
if sum(c_vx > 0) ~= 3
    c_vx(c_vx<0) = 1;
    shitty_cords(end+1) = i;
end
if sum(c_vx > [90 109 90]') > 0
    c_vx(find(c_vx > [90 109 90]')) = size(im.img,find(c_vx > [90 109 90]'));
    shitty_cords(end+1) = i;
end
im.img(c_vx(1),c_vx(2),c_vx(3)) = im.img(c_vx(1),c_vx(2),c_vx(3)) + 1;
end
save_nii(im,ofn)
p_k = [k k k];
spm_smooth(ofn,ofsn,p_k)
disp('done')
%%
%unique(im.img(im.img>0))
%tabulate(im.img(im.img>0))

% Play ground
%b = zeros(size(im.img));
%spm_smooth(ofn,b,p)


% Poisson distrib assumptions, random and independent
%lambda = mean(X)
%variance == mean
%P(X=x) = (lambda^x - 2.718^-lambda) / factorial(x)



% x = poissrnd(4,20,1);
% pd = fitdist(x,'poisson');
% pd.NLogL

