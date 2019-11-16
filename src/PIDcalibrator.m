classdef PIDcalibrator < handle

	properties
		% details about the odorant
		odorant_name
		odorant_volume % in uL
		formula_weight % in g/mol
		density % in g/mL

		% details about the experiment 
		flow_rate_range = linspace(50,500,6);
		MFC_units = 100; % mL/V
		current_flow_rate

		% data
		calibration_data

		% graphics handles
		handles

		% fits 
		flow_rate_to_total_signal_fit
		PID_to_odorant_flux_fit
		PID_to_odorant_flux_fit_r2

	end % end props

	methods
		function [] = calibrate(pc,make_plot)
			
			if nargin < 2
				make_plot = true;
			end

			if make_plot
				pc.handles.calibration_figure = figure('outerposition',[0 0 901 930],'PaperUnits','points','PaperSize',[901 930]); hold on
				subplot(4,3,1:3); hold on
				try
					o = imread([ fileparts(fileparts(which(mfilename)))  filesep 'images' filesep 'calib.png']);
					imagesc(o);
					axis ij
					axis image
					axis off
				catch
				end
			end
			

			pc.buildCalibrationCurve;
			n_moles = ((pc.density/1000)*pc.odorant_volume)/pc.formula_weight;

			if make_plot
				pc.handles.ax.depletion_curves = subplot(4,3,5:6); hold on
				c = parula(length(pc.calibration_data)+1);
				L = {};
				[~,idx] = sort(pc.flow_rate_range);
				for i = 1:length(pc.calibration_data)
					plot(pc.handles.ax.depletion_curves,pc.calibration_data(i).elapsed_time,pc.calibration_data(i).PID,'Color',c(i,:))
					L{i} = [mat2str(pc.flow_rate_range(i)) 'mL/min'];

				end
				legend(L)
				pc.handles.ax.depletion_curves.XLim(1) = 0;
				pc.handles.ax.depletion_curves.YLim(1) = 0;
				xlabel(pc.handles.ax.depletion_curves,'Time (s)')
				ylabel(pc.handles.ax.depletion_curves,'PID (V)')

				pc.handles.ax.air_required = subplot(4,3,4); hold on
				plot([pc.flow_rate_range]/60,[pc.calibration_data.air_req_to_deplete]/1e3,'k+')
				xlabel(pc.handles.ax.air_required,'Flow rate (mL/s)')
				ylabel(pc.handles.ax.air_required,['Air required to' char(10) 'deplete odor (L)'])
				pc.handles.ax.air_required.XLim(1) = 0;
				pc.handles.ax.air_required.YLim(1) = 0;

				pc.handles.ax.total_signal = subplot(4,3,7); hold on
				xx = linspace(0,1.2*max([pc.flow_rate_range]/60),100);
				plot(pc.handles.ax.total_signal,xx,pc.flow_rate_to_total_signal_fit(xx),'r');
				plot(pc.handles.ax.total_signal,[pc.flow_rate_range]/60,[pc.calibration_data.total_signal],'k+')
				xlabel(pc.handles.ax.total_signal,'Flow rate (mL/s)')
				ylabel(pc.handles.ax.total_signal,'Total signal (V s)')
				pc.handles.ax.total_signal.XLim(1) = 0;
				pc.handles.ax.total_signal.YLim(1) = 0;

				pc.handles.ax.cum_signal = subplot(4,3,8:9); hold on
				for i = 1:length(pc.calibration_data)
					y = pc.calibration_data(i).odour_flux;
					plot(pc.calibration_data(i).elapsed_time,1e3*mean(diff(pc.calibration_data(i).elapsed_time))*cumsum(y),'Color',c(i,:))
				end
			
				plot(pc.handles.ax.cum_signal,[0 pc.handles.ax.cum_signal.XLim(2)],1e3*[n_moles n_moles],'k--')
				xlabel(pc.handles.ax.cum_signal,'Time (s)')
				ylabel(pc.handles.ax.cum_signal,'Cumulative odorant (m mol)')
				pc.handles.ax.cum_signal.XLim(1) = 0;
				pc.handles.ax.cum_signal.YLim(1) = 0;

				pc.handles.ax.flow_rate_to_flux = subplot(4,3,10); hold on
				plot(pc.handles.ax.flow_rate_to_flux,[pc.flow_rate_range]/60,[pc.calibration_data.mean_flux]*1e6,'k+')
				xlabel(pc.handles.ax.flow_rate_to_flux,'Flow rate (mL/s)')
				ylabel(pc.handles.ax.flow_rate_to_flux,'Odorant flux (\mu mol/s)')
				pc.handles.ax.flow_rate_to_flux.XLim(1) = 0;
				pc.handles.ax.flow_rate_to_flux.YLim(1) = 0;

				pc.handles.ax.PID_to_flux = subplot(4,3,11); hold on
			end
			% this uses the peak values
			% x = cellfun(@max,{pc.calibration_data.PID});
			% y = [pc.calibration_data.peak_flux]*1e6;

			% now we use mean values -- form 10% to 90% of odor consumed
			y = [pc.calibration_data.mean_flux]*1e6;
			x = [pc.calibration_data.mean_PID];

			if make_plot
				plot(pc.handles.ax.PID_to_flux,x,y,'r+')
				xlabel(pc.handles.ax.PID_to_flux,'PID (V)')
				ylabel(pc.handles.ax.PID_to_flux,'Odorant flux (\mu mol/s)')
			end
			ff = fit(x(:),y(:),'poly1','Lower',[0 0],'Upper',[Inf 0]);
			xx = linspace(0,max(x),100);
			pc.PID_to_odorant_flux_fit = ff;
			if make_plot
				clear l
				l(1) = plot(pc.handles.ax.PID_to_flux,xx,ff(xx),'r');
				

				pc.handles.ax.PID_to_flux.XLim(1) = 0;
				pc.handles.ax.PID_to_flux.YLim(1) = 0;

				figlib.pretty('FontSize',14);

				pc.handles.ax.PID_to_flux.Box = 'off';

				% move some boxes
				pc.handles.ax.total_signal.Position(2) = .31;
				pc.handles.ax.cum_signal.Position(2) = .31;
				pc.handles.ax.flow_rate_to_flux.Position(2) = .08;
				pc.handles.ax.PID_to_flux.Position(2) = .08;

				[~,lh] = figlib.label('x_offset',-0.01,'y_offset',0);
				lh(1).Position(1:2) = [.25 .935];
				uistack(lh(1),'top')
			end

			
			pc.PID_to_odorant_flux_fit_r2 = (corr(pc.PID_to_odorant_flux_fit(x),y(:)))^2;

		end % end calibrate 


		function [] = buildCalibrationCurve(pc)


			% sort 
			[pc.flow_rate_range,idx] = sort(pc.flow_rate_range);
			pc.calibration_data = pc.calibration_data(idx);


			% compute the number of moles we have depleted
			n_moles = ((pc.density/1000)*pc.odorant_volume)/pc.formula_weight;

			% reformat x axis to a more sensical format
			if isfield(pc.calibration_data(1),'elapsed_time')
			else
				for i = 1:length(pc.calibration_data)
					
					t = pc.calibration_data(i).t;
					t = t - t(1);
					t = datevec(t);
					pc.calibration_data(i).elapsed_time = t(:,4)*60*60 + t(:,5)*60 + t(:,6);
				end
			end

			% remove minimum
			for i = 1:length(pc.calibration_data)
				pc.calibration_data(i).PID = pc.calibration_data(i).PID - min(pc.calibration_data(i).PID);
			end

			% amount of air required to deplete 99% of odour
			for i = 1:length(pc.calibration_data)
				pc.calibration_data(i).air_req_to_deplete = (pc.flow_rate_range(i)/60)*(pc.calibration_data(i).elapsed_time(find(cumsum(pc.calibration_data(i).PID) > sum(pc.calibration_data(i).PID)*.99,1,'first')));
			end

			% now compute the integrated PID signal under these different flow rate conditions
			for i = 1:length(pc.calibration_data)
				pc.calibration_data(i).total_signal = mean(diff(pc.calibration_data(i).elapsed_time))*sum((pc.calibration_data(i).PID)); % this is the time integrated PID signal
			end

			x = pc.flow_rate_range/60; % now in mL/s
			pc.flow_rate_to_total_signal_fit = fit(x(:),corelib.vectorise([pc.calibration_data.total_signal]),'smoothingspline');

			for i = 1:length(pc.calibration_data)
				pc.calibration_data(i).odour_flux = (pc.calibration_data(i).PID*(n_moles)/pc.flow_rate_to_total_signal_fit(pc.flow_rate_range(i)/60));
				pc.calibration_data(i).peak_flux = max(pc.calibration_data(i).odour_flux);
				at = find((pc.calibration_data(i).elapsed_time - pc.calibration_data(i).elapsed_time(1)) > 400,1,'first');

				% compute the mean flux from 10% to 90% of the odor consumed
				cum_flux = cumsum(pc.calibration_data(i).odour_flux);
				cum_flux = cum_flux/cum_flux(end);
				a = find(cum_flux > .1,1,'first');
				z = find(cum_flux > .9,1,'first');

				pc.calibration_data(i).mean_flux = mean(pc.calibration_data(i).odour_flux(a:z));
				pc.calibration_data(i).mean_PID = mean(pc.calibration_data(i).PID(a:z));

			end


		end % end buildCalibrationCurve

		function [] = save(pc)
			% shut down all figures
			pc.cleanup_handles;
			savename = [pc.odorant_name '_calib.mat'];
			save(savename,'pc')

		end % end save

		function cleanup_handles(nsb,~,~)
			if isempty(nsb.handles)
				return
			end
			fn = fieldnames(nsb.handles);
			for i = 1:length(fn)
				try
					delete(nsb.handles.(fn{i}));
				catch
				end
				try
					rmfield(nsb.handles,fn{i});
				catch
				end
			end
			try
				nsb.handles = [];
			end
			try
				rmfield(nsb,'handles');
			end
		end % end cleanup_handles

		function [t, PID] = deplete(pc)

			% calibrates PID using the depletion technique 
			% make sure you have the experiment set up correctly:
			% output1: Analogue Output to main air MFC
			% output2: anaogue output to secondary air MFC
			% and nothing else
			%
			% make sure kontroller is configures so that PID and MFC are acquired 

			% figure out which flow rate to use
			try
				for i = 1:length(flow_rate_range)
					temp = pc.calibration_data(i).PID;
				end
			catch
				pc.current_flow_rate = pc.flow_rate_range(i);
			end


			% make control paradigms
			ControlParadigm.Name = 'Zero';
			ControlParadigm.Outputs = zeros(2,10);
			ControlParadigm.Outputs(1,:) = 2;


			% measure the zero point
			kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',ones(3,1),'w',10);
			first_blank = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',ones(5,1),'w',10);

			PID = [];
			MFC = [];
			t = [];
			this_PID = Inf;
			i = 1;

			% make the figure and the plots
			pc.handles.depeletion_curve = figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700],'NumberTitle','off','Name',['Depletion curve at ' oval(pc.current_flow_rate) 'mL/min']); hold on; 
			h = plot(NaN,NaN,'k+');

			h2 = plot(NaN,NaN,'r');
			h3 = plot(NaN,NaN,'r--');

			% make control paradigms
			ControlParadigm.Name = 'Flow';
			ControlParadigm.Outputs = ones(2,10);
			ControlParadigm.Outputs(1,:) = 2;
			ControlParadigm.Outputs(2,:) = pc.current_flow_rate/pc.MFC_units;
			c = 1;
			while (this_PID > (mean(mean([first_blank.PID]')) + 3*std(mean([first_blank.PID]'))) || c < 10) && isvalid(pc.handles.depeletion_curve)
			    clc
			    try
			        data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10);
			        this_PID = mean(data.PID(:));
			        this_MFC = mean(data.MFC(:));
			        PID = [PID; this_PID];
			        MFC = [MFC; this_MFC];
			        t = [t; now];
			        set(h,'XData',t,'YData',PID);

			        % update the baseline too
			        h2.XData = [min(t) max(t)];
			        h3.XData = [min(t) max(t)];
			        h2.YData = [mean(mean([first_blank.PID]')) mean(mean([first_blank.PID]'))];
			        h3.YData = mean(mean([first_blank.PID]'))+[3*std(mean([first_blank.PID]')) 3*std(mean([first_blank.PID]'))];

			        datetick('x')
			        c = c + 1;

			        save('temp.mat','t','PID','MFC')
			    catch
			        warning('Something went wrong when trying to measure the stimulus')
			    end
			end % end while loop

			% shut it down

			% make control paradigms
			ControlParadigm.Name = 'Zero';
			ControlParadigm.Outputs = zeros(2,10);
			kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10);
			ControlParadigm.Outputs(1,:) = 2;

			% export into the main object
			load_here = find(pc.flow_rate_range == pc.current_flow_rate);
			pc.calibration_data(load_here).t = t;
			pc.calibration_data(load_here).PID = PID;
			pc.calibration_data(load_here).MFC = MFC;

			pc.save;
		end % end deplete
	end % end methods

end % end classdef 