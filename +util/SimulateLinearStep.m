function [decoupledResults,coupledResults] = SimulateLineareStep(param, saveName, applyYLim,simulationTime_S,ref)
    if nargin < 4
        simtime = 10*60;
    else
        simtime = simulationTime_S;
    end
    stepTime = 200;
    ref_start = 20;
    if nargin < 5
        ref = 1;
    end
    ref_stop = ref_start+ref;
    x0 = [param.ctrl.T_wOP.'; param.ctrl.T_aOP.'; zeros(4,1)];
    % Decoupled model
    [t,x] = ode15s(@(t,x)model.LinearDynamics(t,x,param,1,stepTime,ref),[0 simtime],x0);
    decoupledResults.time = t;
    decoupledResults.Tw = x(:,1:4);
    decoupledResults.Ta = x(:,5:8);
    decoupledResults.T = x(:,9:12);

    % Coupled model
    [t,x] = ode15s(@(t,x)model.LinearDynamics(t,x,param,0,stepTime,ref),[0 simtime],x0);
    coupledResults.time = t;
    coupledResults.Tw = x(:,1:4);
    coupledResults.Ta = x(:,5:8);
    coupledResults.T = x(:,9:12);
    q = zeros(4,length(t));
    for(i = 1:length(t))
        [~,q(:,i)]=model.LinearDynamics(t(i),x(i,:),param,0,stepTime,ref);
    end
    coupledResults.q=q.';

    step_x = linspace(0,simtime,simtime);
    step_y = [ones(1,stepTime)*ref_start,ones(1,simtime-stepTime)*ref_stop];

    figure()
    f = tiledlayout('vertical');
    for i=1:4
        nexttile
        hold on
        plot(decoupledResults.time,decoupledResults.Ta(:,i)-273.15,'blue')
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

