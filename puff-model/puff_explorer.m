
%%
% This script launches an interactive explorer which allows you to manipulate
% parameters in this model

close all
clearvars


dm = dataManager;
load(dm.getPath('69a8e45311519ae0fd22db1a17d4cc46'))

Model = TwoTubesX;

use_this = 1;

Model.Stimulus = fd(use_this).stimulus;
Model.Response = fd(use_this).response;

Model.manipulate;

set(gca,'YLim',[0 1])
xlabel('Time')