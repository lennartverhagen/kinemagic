function convert_Marius2Kinemagic

load('/home/action/lenver/Data/Marius/orig/polhemusdata_1.mat');

% init data structure
data = [];
ntrl = length(polhemusdata);
for t = 2:21
    nsamp = size(polhemusdata(t).speed_x,1);
    data.vel{t-1} = reshape([polhemusdata(t).speed_x polhemusdata(t).speed_y polhemusdata(t).speed_z]',[1 3 nsamp]);
    data.time{t-1} = (0:(nsamp-1))/240;
end
data.label = {'hand'};
data.derivaxis = {'x','y','z'};
data.dimord = 'marker_axis_time';

% init cfg structure
cfg = [];
cfg.exp = 'Marius';
cfg.dataset = 'test';
cfg.subj = '01';
cfg.sessions = {''};
cfg.dir.root = '/home/action/lenver/Data/Marius';
cfg.dir.preproc = '/home/action/lenver/Data/Marius/preproc';
cfg.trl = ones(20,1);
cfg.vars = {'ones'};

% save
save(fullfile(cfg.dir.preproc,'01_test_CFG'),'cfg');
save(fullfile(cfg.dir.preproc,'01_test_DATA'),'data');