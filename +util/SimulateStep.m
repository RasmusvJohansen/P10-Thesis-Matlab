function [decoupledResults,coupledResults] = SimulateStep(param, saveName, applyYLim,simulationTime_s)
    
    if (nargin < 4)
    simtime = 10*60;
    else
    simtime = simulationTime_s;
    end
    stepTime = 200;
    ref_stop = 21;
    ref_start = 20;
    % Decoupled sim
    for i=1:4
        x0 = [param.ctrl.T_wOP(i); param.ctrl.T_aOP(i); 0];
        [t,x] = ode23s(@(t,x)model.decoupledDynamics(t,x,i,param,stepTime),[0 simtime],x0);
        decoupledResults(i).time = t;
        decoupledResults(i).Tw = x(:,1);
        decoupledResults(i).Ta = x(:,2);
        decoupledResults(i).T = x(:,3);
    end
    %Coupled sim
    x0 = [param.ctrl.T_wOP.'; param.ctrl.T_aOP.'; zeros(4,1)];
    [t,x] = ode23s(@(t,x)model.coupledDynamics(t,x,param,stepTime),[0 simtime],x0);
    coupledResults.time = t;
    coupledResults.Tw = x(:,1:4);
    coupledResults.Ta = x(:,5:8);
    coupledResults.T = x(:,9:12);

    
    step_x = linspace(0,simtime,simtime);
    step_y = [ones(1,stepTime)*ref_start,ones(1,simtime-stepTime)*ref_stop];

    figure()
    f = tiledlayout('vertical');
    for i=1:4
        nexttile
        hold on
        plot(decoupledResults(i).time,decoupledResults(i).Ta-273.15,'blue')
        plot(coupledResults.time,coupledResults.Ta(:,i)-273.15,['red','--'])
        plot(step_x,step_y,"k--")
        if i==1
            legend('Decoupled','Coupled',"Ref",Location='northeast')
        end
        if (applyYLim)
            ylim([19.8,21])
        end
        hold off
    end
    xlabel(f,'Time [s]')
    ylabel(f,'Air temperature [^{o}C]')
    saveas(gcf,strcat('Images/',saveName),'epsc')
end

