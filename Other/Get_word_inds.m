clear all
load('/Users/aidasaglinskas/Desktop/Neurosynth/NeuroSynth_labels.mat')

% PARAMS, words
words = {'face' 'facial'};
to_drop = {'al' 'surface'};

l = cellfun(@num2str,l,'UniformOutput',0);
for w = 1: length(words);
t = arrayfun(@(x) findstr(l{x},words{w}),1:length(l),'UniformOutput',0);
tt{w} = cellfun(@isempty,t) == 0;
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
save('/Users/aidasaglinskas/Desktop/Neurosynth/face_inds.mat','list')

