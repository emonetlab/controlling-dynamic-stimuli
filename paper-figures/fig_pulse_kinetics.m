

%% Pulse Kinetics
% In this figure, I show the diversity of pulse kinetics and fit a model to reproduce these kinetics


% specify model to use
addpath('../puff-model/')
Model = TwoTubes;

dm = dataManager;
load(dm.getPath('0b91a233caac557fbcbfa03f5324f9ca'))


% remove PID baseline from all the data
for i = 1:length(alldata)
	for j = 1:size(alldata(i).PID,2)
		alldata(i).PID(:,j) = alldata(i).PID(:,j) - mean(alldata(i).PID(1:500,j));
	end
end

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
clear ax fd
c = 30;
sp = 1:30;

for i = length(alldata):-1:1
	ax(i) = subplot(5,6,sp(i)); hold on; c = c - 1;
	title(alldata(i).odour_name)
	S = alldata(i).PID;
	conc = alldata(i).conc;
	if ~any(isnan(conc))
		% more than one concenration, pick the max
		S = S(:,conc == max(conc));

	end

	% normalise
	for j = 1:size(S,2)
		S(:,j) = S(:,j)/max(S(:,j));
	end

	time = 1e-3*(1:length(S));
	plot(time,S,'Color',[.5 .5 .5])
	set(ax(i),'XLim',[0.5 3],'YLim',[-.1 1.05])

	all_S = S;

	% save for fitting
	S = mean(S,2);
	z = 3e3;
	all_S = all_S(1:z,:);
	all_S = all_S(1:10:end,:);
	S = S(1:z); % ten seconds enough
	S = S(1:10:end); % now at 10 ms timestep
	S = S/max(S);
	S(S<0) = 0;
	fd(i).response = S;
	P = 0*S;
	t = 1:length(P);
	P1 = 1 - exp(-(t-98)/2.5);
	P2 = exp(-(t-147)/3.8);
	P(100:150) = P1(100:150);
	P(151:end) = P2(151:end);
	fd(i).stimulus = P;
	fd(i).response(1:90) = NaN;



end


% show the simulations too

savename = [class(Model) '.fitparams'];
load(dm.getPath('f1bd124ccba607c8230dfda54fd4bedb'),'-mat')


for i = length(fd):-1:1

	% pick the best r^2
	this_r2 = all_r2(i,:);
	[~,pick_me] = max(this_r2);


	Model.Stimulus = fd(i).stimulus;
	Model.Parameters = p(i,pick_me);
	try
		Model.evaluate;
		plot(ax(i),time(1:10:3e3),Model.Prediction,'r')
	catch
		
	end
	axis(ax(i),'off')
end

figlib.pretty('FontSize',13)