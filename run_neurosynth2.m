%CMD+ENTER to run portion of the script 
%% Part 0: Load up datasets and specify parameters (RUN ONCE)
clear all
clc
tic
% EDIT: POINT SCRIPT TO THE FOLDER
p.root_folder = '/Users/aidasaglinskas/Desktop/Matlab_Neurosynth';

addpath(genpath(p.root_folder))
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

u_space = find(cellfun(@(x) strcmp(x,'UNKNOWN'),Neurosynth_all.Database.space));
u_space_studies = unique(Neurosynth_all.Database.id(u_space));
disp(sprintf('%s percent of studies have UNKNOWN space, dropping them',num2str(length(unique(Neurosynth_all.Database.id(u_space))) / length(unique(Neurosynth_all.Database.id)) * 100)))
Neurosynth_all.features(find(ismember(Neurosynth_all.features(:,1),u_space_studies)),:) = [];
Neurosynth_all.Database(u_space,:) = [];
%% Part 1; Choose word list (RUN UNTIL SATISFIED)
% choose words to include and exclude until satisfied, proceed
% EDIT: CHOOSE WORDS TO INCLUDE EXCLUDE
words = {'face' 'facial' 'person' 'people'} %{'face' 'facial' 'person' 'people'};
to_drop = {'al' 'surface' 'rs' 'fa' 'interpersonal'} %{'al' 'surface' 'rs' 'interpersonal'};

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
d = cellfun(@(x) x',d,'UniformOutput',0)';d = [d{:}]';
disp(list.words([d]));
list.words([d]) = [];
list.inds([d]) = [];
end
disp('FINAL LIST')
disp(list.words)
% save('/Users/aidasaglinskas/Desktop/Neurosynth/face_inds.mat','list')
%% Part 2: Get the coordinates, map to brain, smooth etc (RUN ONCE)
temp_fn = fullfile(p.root_folder,p.name_of_MRI_template);
ofn = fullfile(p.root_folder,p.meta_analysis_output_fn);
ofsn = fullfile(p.root_folder,p.meta_analysis_output_fn_smoothed);
shitty_coords = [];

% override for all meta
disp('Overriden for all meta')
list.inds = 2:size(Neurosynth_all.features,2)


n_reduced = Neurosynth_all.features(:,[1 list.inds]);
% for all meta 
                    %%%% plot n_reduced
%                     f = figure(1);
%                     mat = n_reduced(:,2:end);
%                     imagesc(mat');
%                     f.CurrentAxes.YTick = 1:size(mat,2);
%                     f.CurrentAxes.YTickLabel = list.words;
%                     f.CurrentAxes.FontSize = 14;
%                     f.CurrentAxes.FontWeight = 'bold';
%                     f.CurrentAxes.YTickLabelRotation = 0;

% Find studies containing the keywords;
temp.p = sum(n_reduced(:,2:size(n_reduced,2)),2); % sum of all feature vectors 
studies_matching_ind = find(temp.p); %find studies with non zero temp.p
studies_matching_IDs = n_reduced(studies_matching_ind,1); % % get ID of those studies 

%db_inds = find(ismember(Neurosynth_all.Database.id,studies_matching_IDs));
%b_corrds = find(cellfun(@(x) strcmp(x,'UNKNOWN'),Neurosynth_all.Database.space(db_inds)));
%%Neurosynth_all.Database.space(db_inds(b_corrds));
%disp(sprintf('%d out of %d coordinates reported have UNKNOWN space (%s percent)',length(db_inds(b_corrds)),length(db_inds),num2str(length(db_inds(b_corrds)) / length(db_inds)* 100)))
%db_inds(b_corrds) = []; % drop coords reported in unknown space; 




disp(sprintf('%d Studies found matching criteria',length(studies_matching_IDs)))
study_weights = temp.p(studies_matching_ind);
    % Find indeces, coords, name of journal etc, matching query
    %t.m = ismember(Neurosynth_all.Database.id,studies_matching_IDs);
    %db_inds = find(ismember(Neurosynth_all.Database.id,studies_matching_IDs));
                % Load Template Brain
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


          % pre allocate
          im.img = zeros([size(im.img) length(studies_matching_IDs)]);
          %im.img_smoothed = zeros([size(im.img) length(studies_matching_IDs)]);
                
disp('Mapping')
tic
for study_ind = 1:length(studies_matching_IDs)
    report_vect = [[1:round(length(studies_matching_IDs)/10):length(studies_matching_IDs)] 1:100:length(studies_matching_IDs)];
    if ismember(study_ind,report_vect);disp(sprintf('%s percent done',num2str([(study_ind / length(studies_matching_IDs) * 100)]))); end
    
inds = find(Neurosynth_all.Database.id == studies_matching_IDs(study_ind));

Coords = [Neurosynth_all.Database.x(inds) Neurosynth_all.Database.y(inds) Neurosynth_all.Database.z(inds)];
for c = 1:length(inds)

if strcmp(Neurosynth_all.Database.space(inds(1)),'MNI')
c_mm = Coords(c,:);
else
c_mm = tal2mni(Coords(c,:));
end
c_vx = inv(m)*[c_mm 1]';
c_vx = round(c_vx(1:3));


if sum(c_vx <= 0) >= 1 | [c_vx(1) > size(im.img,1) c_vx(2) > size(im.img,2) c_vx(3) > size(im.img,3)] % if voxel coords are negaive or outside the template; 
shitty_coords(end+1,:) = c_vx';
else
im.img(c_vx(1),c_vx(2),c_vx(3),study_ind) = study_weights(study_ind); % adds weight to the coordinate

end %ends shitty coord if
end%end coords
end %ends study loop
disp('Done mapping')
toc
%%
disp('Smoothing Individual Studies')
im.img_smoothed = zeros(size(im.img));

%spm_smooth(P,Q,s,dtype)
%v = zeros(size((im.img(:,:,:,s))));
disp('Smoothing')
for s = 1:size(im.img,4); 
    if ismember(s,round([1:size(im.img,4)/10:size(im.img,4)])); disp(sprintf('%s percent done',num2str(s / size(im.img,4) * 100)));end
    v = zeros(size(im.img(:,:,:,s)));
spm_smooth(im.img(:,:,:,s),v,[k k k]);
im.img_smoothed(:,:,:,s) = v;
end
disp('done')
%figure;image3(im.img_smoothed(:,:,:,100))
toc
pwd
%%
opts.plot_rnd = 0;
if opts.plot_rnd == 1
t_s = 112
f = figure(5)
subplot(1,2,1)
image3(im.img(:,:,:,t_s))
disp('1')
drawnow
subplot(1,2,2)
image3(im.img_smoothed(:,:,:,t_s))
disp('2')
drawnow
end
%%
% disp('Exporting')
% s = 100;
% im_copy = im;
% im_copy = rmfield(im_copy,'img_smoothed')
% im_copy = rmfield(im_copy,'img')
% im_copy.img = im.img_smoothed(:,:,:,s)
% 
% t_ofn = fullfile('/Users/aidasaglinskas/Desktop/scans/',[datestr(datetime) '.nii']);
% save_nii(im,t_ofn)
% disp('done')
%%

% % Get Rid of non MNI coords; 
% t.p1 = cellfun(@(x) strmatch(x,'UNKNOWN'),Neurosynth_all.Database.space(db_inds),'UniformOutput',0);
% t.pp = find(cellfun(@isempty,t.p1) == 1);
% unique(Neurosynth_all.Database.space(db_inds(t.pp)))
% 
% 
% 
% t.p2 = find(cellfun(@isempty,t.p1) == 1);
% 
% 
% 
% disp(sprintf('%d studies dropped due to non-MNI coordinates',length(unique(Neurosynth_all.Database.title(db_inds(t.p2))))));
% db_inds(t.p2) = [];
% disp(sprintf('%d studies included in meta-analysis',length(unique(Neurosynth_all.Database.title(db_inds)))));
% meta.db_inds = db_inds;
% meta.x = Neurosynth_all.Database.x(meta.db_inds);
% meta.y = Neurosynth_all.Database.y(meta.db_inds);
% meta.z = Neurosynth_all.Database.z(meta.db_inds);
% % map to a brain
% disp('Mapping to brain')
% im = load_nii(temp_fn); %load
% im.img(:) = deal(0); %wipe
% % Transformation Matrix from mm to vx
% m = [im.hdr.hist.srow_x;im.hdr.hist.srow_y;im.hdr.hist.srow_z;[0 0 0 1]];
% m(:,4) = m(:,4) - sum(m(:,1:3),2); % adjustment fuckery
% try
% a = cosmo_fmri_dataset(temp_fn);
% %fprintf('looks like you have cosmo installed..\n..lemme run this additional check real quick')
% mat = a.a.vol.mat;
% if all(m(:) == mat(:));
% %    fprintf('\ntransformation matrix looks good nevermind')
% else warning('hmm.. transformation matrix might be fucked, check if the two below are similar')
% disp('Transformation matrix from nifti header')
% disp(m)
% disp('CosmoMVPA''s transformation matrix')
% disp(mat)
% disp('Fuck it, I''ll use cosmo''s instead')
% m = mat;
% end
% catch
% end
% %
% disp('Mapping Coordinates')
% shitty_cords = []; % coords that are either out of the template or negative in voxel space
% for i = 1:length(meta.db_inds)
% c_mm = [meta.x(i) meta.y(i) meta.z(i)];
% c_vx = inv(m)*[c_mm 1]';
% c_vx = round(c_vx(1:3));
% if sum(c_vx > 0) ~= 3 % are the coords that are negative in voxel space?
%     c_vx(c_vx<0) = 1;
%     shitty_cords(end+1) = i;
% end
% if sum(c_vx > [90 109 90]') > 0 % are they out of range/ 
%     c_vx(find(c_vx > [90 109 90]')) = size(im.img,find(c_vx > [90 109 90]'));
%     shitty_cords(end+1) = i;
% end
% im.img(c_vx(1),c_vx(2),c_vx(3)) = im.img(c_vx(1),c_vx(2),c_vx(3)) + 1;
% end
% disp(sprintf('%d coordinate entries were outside the brain temp',length(shitty_cords)))
% disp(round([meta.x(shitty_cords) meta.y(shitty_cords) meta.z(shitty_cords)]))
% save_nii(im,ofn)
% p_k = [k k k];
% spm_smooth(ofn,ofsn,p_k)
% disp('done')