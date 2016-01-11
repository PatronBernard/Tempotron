function [V]=Tempotron(n_afferents, t)
    %n_a    Number of afferents, e.g. 5 or 10
    %t      Vector with discretized time, e.g. t=linspace(0,0.5,1000);  
    
    w=rand(10); %For now we just generate random weights
    V_thresh=1; %Threshold voltage
    
    tau=0.015;  %Membrane integration time constant
    
    V_sub=zeros(size(t)); %Submembrane voltage
    
    %Generate a random pattern and calculate the output with random
    %weights.
    for i=1:n_afferents
        Kn(i,:)=reshape(K(t,rand(3)/2),1,length(t));
        V_sub=V_sub+w(i)*Kn(i,:);
    end
    
    %See if the threshold voltage was crossed
    [~,idx]=find(V_sub > V_thresh);
    V=V_sub;
    
    %If so, reduce the voltage smoothly to 0.
    if ~isempty(idx)
        idx=min(idx);
        V(idx:end)=V_sub(idx).*exp(-(t(idx:end)-t(idx))/tau);
    end
    if nargout==0
        close all; 
        %Plot the individual afferent activity in an image
        subplot(2,1,1)
        imagesc(Kn); colormap('hot');
        set(gca,'YTick',1:n_afferents);
        xlabel('t (ms)'); ylabel('n');
        
        %Plot the subthreshold membrane voltage and the neuron's response
        subplot(2,1,2)
        plot(t,V,'k','LineWidth',2);   
        
        hold on
        plot(t,V_sub,'Color',[0.5 0.5 0.5])
        legend('Output voltage','Subthreshold membrane voltage');
        
        plot(t,V_thresh.*ones(size(t)),'k-.');
        xlabel('t (ms)'); ylabel('V (AU)');
        hold off;
        
        set(gcf,'Position',[125 189 1115 789])
    end
end

function V=K(t,t_spikes)
    %For some time interval <t> (an N-element vector) this function will
    %generate exponentially decaying spikes at times <t_spikes>. 
    
    %Model parameters
    tau=0.015;      %Membrane integration time constant
    tau_s=tau/4;    %Synaptic current time constant
    V0=2.12;        %Normalization (see paper)
    
    %The step function is required to acquire the correct form of the
    %spikes (i.e. 0 everywhere except after some spike). 
    V=zeros(1,length(t));
    for i=1:length(t_spikes)
        V=V+heaviside(t-t_spikes(i)).*(exp(-(t-t_spikes(i))/tau)-exp(-(t-t_spikes(i))/tau_s));
    end
    V=V0.*V;
end

%Written by Jan Morez
%Computational Neuroscience, Miniproject 2016
%University of Antwerp