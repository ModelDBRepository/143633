function varargout = generate_stim_poisson_harmaline(genes,tspan,ngen)
 
    rand('state',100*sum(clock));
    randn('state',100*sum(clock));
    
    if nargin<3
        ngen = 1;
    end

     if nargin<1
		genes = 1;
    end
	if nargin<2
		tspan = [0 2000 12000 14000 0.01];
    end		
	
	if length(genes)==1
		Ni = genes;   % Number of individuals
		Ng = 10;			% Number of genes
		
		ipi_min =   1;  % ms
		ipi_max =  20;  % ms
		G = rand(Ni,Ng)*(ipi_max-ipi_min) + ipi_min;
	else
		if iscell(genes)
			Ni = length(genes);
			Ng = 0;
			for n=1:Ni
				Ng = max(Ng,length(genes{n}));
            end
			for n=1:Ni
				gene = genes{n};
				if length(gene<Ng);
					gene(end+1:Ng) = mean(gene);
                end
				G(n,:) = gene;
            end
		else
			[Ni,Ng] = size(genes);
			G = genes;
        end

        if isempty(genes)
            G = (1000/130)*ones(1,10);
            [Ni,Ng] = size(G);
        end
        
    end

    mean_rate = [185*ones(1,6) 0.000001];
    dist_type = [0 1 2 3 4 5 0]; % 0 = regular, 1 = Log-normal, 2 = UniPeak, 3 = Bimodal, 4 = Absence, 5 = Presence
    Ni = length(mean_rate);
    for i=1:Ni
        G(i,:) = (1000/mean_rate(i))*ones(1,Ng);
    end
    

    [POISSON,HARM,STIM_trains] = create_poisson_stim(G,tspan,dist_type);
    
    eval(['save Genes_' num2str(ngen) '.mat G STIM_trains']) 
    
    varargout{1} = tspan;
	varargout{2} = G;
    varargout{3} = STIM_trains;
    varargout{4} = POISSON;
    varargout{5} = HARM;

 
 
nTC = 490;
nCER = 233;
nCTX = 233;
nRN = 233;
nTIN = 113;
NTOT = nTC+nCER+nCTX+nRN+nTIN;

Drop_ax1 = ones(5,100);
load Drop1_ind.mat

for nc=1:100
    
    miniDrop1 = find(Drop1((nc-1)*NTOT+1:(nc)*NTOT)==1);
    if isempty(miniDrop1)
    else
        if min(miniDrop1) <= nTC;
            Drop_ax1(1,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER & miniDrop1 > nTC));
            Drop_ax1(2,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER + nCTX & miniDrop1 > nTC + nCER));
            Drop_ax1(3,nc) = 0;
        else
        end
        if ~isempty(find(miniDrop1 <= nTC + nCER + nCTX + nRN & miniDrop1 > nTC + nCER + nCTX));
            Drop_ax1(4,nc) = 0;
        else
        end
        if max(miniDrop1) > nTC + nCER + nCTX + nRN;
            Drop_ax1(5,nc) = 0;
        else
        end
    end
end

load PHI_by_cell.mat

for n_cell = 1:100;
    f_name = (['phi_pop1_cell_' num2str(n_cell) '.dat']);
    fid = fopen(f_name,'wb');
    fwrite(fid,PHI_NRN1{n_cell,1},'double');
    fclose(fid);
end;

Drop_ax1 = Drop_ax1(1:500)
fid = fopen('Drop_ax1.dat','wb')
fwrite(fid,Drop_ax1,'double');
fclose(fid);
 
end

function varargout = create_poisson_stim(G,tspan,dist_type)
% % tspan is a vector [x y z], where:
%   x = start time (usually 0)
%   y = end time (in ms)
%   z = dt
% % G is a matrix of interstimulus intervals, where each row is for a different stimulus condition
%   e.g., G = [10 10 10; 20 20 20] would create two constant frequency stimulus trains, the first at 100 Hz
%           and the second at 50 Hz.
%   
%   or, G = [10 10 10 10 10 10 10 100] would create a pause train with 5 ISIs at 100 Hz and 1 at 10 Hz.
%
% 	
% size of G = N x M, where N = number of different stimulus trains and M =
% 	number of ISIs in repeated sequence.
% tspan = [t_start t_stim_on t_stim_off t_end dt]
 
    t_start = tspan(1);
    t_stim_on = tspan(2);
    t_stim_off = tspan(3);
    t_end = tspan(4); % Times in ms
    dt = tspan(5); % Adjust dt so that PW is a multiple of dt
    PW = 0.1; % 100 us pulses
    t = dt:dt:t_end;
    p_rate = 20; % rate of poisson process
    p_tref = 0.001; % 
%     delay = 1.2; 
    Template = 2*ones(1,PW/dt); % Monophasic: uses amplitude of 2 to activate exc and inh axons
    
    ratio = 10;
    Stim_template = [ones(1,PW/dt) -(1/ratio)*ones(1,PW*ratio/dt)]; % Biphasic
    
    
    HARM = zeros(1,round(t_end/dt));
    POISSON = zeros(1,round(t_end/dt));

 
    % Generate Poisson pulse times
    [poisson] = create_poisson(t_end/1000,p_rate,0.001);
    poisson = 1000*poisson;


    % Generate Harmaline pulse times
    isi_intra = 7; % 6
    isi_long = 95; % 55
    num_intra = 11; % 10
    total_int = isi_intra*num_intra+isi_long;
    h = [isi_intra*ones(1,num_intra) isi_long];
    repeats = floor(t_end/total_int);
    remainder = total_int*(t_end/total_int - repeats);
    portion_rep = floor(remainder/isi_intra);


    if portion_rep > 10
        portion_rep = 10;
    else
    end

    h_full = [repmat(h, [1,repeats]) isi_intra*ones(1,portion_rep)];
    harm = cumsum(h_full);


     for i=1:length(harm);
        A=round(harm(i)/dt);
        beg = A;
        last= A+length(Template)-1;
        HARM(beg:last) = Template;
    end

    for i=1:length(poisson);
        A=round(poisson(i)/dt);
        beg = A;
        last= A+length(Template)-1;
        POISSON(beg:last) = Template;
    end


    % Initialize Stim Trains
    STIM_trains = cell(size(G,1),1);
    
    for i=[1 7]
    
        s = sum(G(i,:));
        ipi = repmat(G(i,:),[1,round(t_end/s)+2]);
%         subplot(3,2,i)
        hist(ipi,25)
        stim = cumsum(ipi)-3;
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end

    % Create Uniform Distribution Stimulus
    for i=2
        delta = 0.312929;
        LU = log10(90) + 2*delta*rand(round(t_end),1);
        F = 10.^LU;
        ipi = 1000./F;
        subplot(3,2,i)
        hist(ipi,25)
        stim = t_stim_on + cumsum(ipi);
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end
    
    % Create UniPeak Distribution Stimulus
    for i=3
        delta_tall = 2*0.312929;
        LUtall = log10(185) - delta_tall + 2*(delta_tall)*rand(round(t_end),1);
        X = rand(round(t_end),1);
        p185 = 0.2833887;
        LUtall(X<p185) = log10(185);

        F = 10.^LUtall;
        ipi = 1000./F;
        subplot(3,2,i)
        hist(ipi,25)
        stim = t_stim_on + cumsum(ipi);
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end

    % Create Bimodal Distribution Stimulus
    for i=4
        delta = 0.312929;
        delta_tall = 2*0.312929;
        LUlow = log10(185) - delta_tall + delta*rand(round(t_end),1);
        LUup = log10(185) + delta + delta*rand(round(t_end),1);
        X = rand(round(t_end),1);

        LUbim = zeros(t_end,1);
        for j=1:t_end;
            if X(j)<0.5
                LUbim(j) = LUlow(j);
            else
                LUbim(j) = LUup(j);
            end
        end
            
        F = 10.^LUbim;
        ipi = 1000./F;
        subplot(3,2,i)
        hist(ipi,25)
        stim = t_stim_on + cumsum(ipi);
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end

    % Create Absence Distribution Stimulus
    for i=5    
        ipi = repmat([5.070496576*ones(35,1); 50.70496576],round(t_end/36),1);
        subplot(3,2,i)
        hist(ipi,25)
        stim = cumsum(ipi)-3;
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end
    
    % Create Presence Distribution Stimulus
    for i=6    
        ipi = repmat([7.009943539*ones(25,1); 3.504971769*ones(15,1)],round(t_end/40),1);
        subplot(3,2,i)
        hist(ipi,25)
        stim = cumsum(ipi)-3;
        stim = stim(stim < t_stim_off & stim > t_stim_on);
        eval(['STIM_' num2str(i) ' = zeros(1,round(t_end/dt));'])
        for j=1:length(stim)
            A=round(stim(j)/dt);
            beg = A;
            last= A+length(Stim_template)-1;
            eval(['STIM_' num2str(i) '(beg:last) = Stim_template;'])
        end

    eval(['STIM = STIM_' num2str(i) ';'])
    eval(['clear STIM_' num2str(i) ';'])
    STIM = STIM(1:length(t));
    file = ['Stim_' num2str(i) '.dat'];
    fid = fopen(file,'wb');
    fwrite(fid,STIM,'double');
    fclose(fid);
    STIM_trains{i} = STIM;

    end
   
    
    % Write Harmaline Data to File for use in NEURON
    fid = fopen('Harmaline.dat','wb');
    fwrite(fid,HARM,'double');
    fclose(fid);


    % Write Poisson Data to File for use in NEURON
    fid = fopen('Poisson.dat','wb');
    fwrite(fid,POISSON,'double');
    fclose(fid);
    
    varargout{1} = POISSON;
    varargout{2} = HARM;
    varargout{3} = STIM_trains;
end




function IPI = create_skip(N,n_skip,n_during);
    if nargin < 3
        n_during = 27;
    end

    if nargin < 2
        n_skip = 6;
    end

    if nargin < 1
        N = 10000;
    end
    
    intra_ipi = ((1000^(n_during+1))/((130^(n_during+1))*(n_skip+1)))^(1/(n_during+1));
    skip_ipi = (n_skip + 1) * intra_ipi;
    ipi = [intra_ipi * ones(n_during,1); skip_ipi];
    IPI = repmat(ipi,[ceil(N/(n_during+1)) 1]);
  
end




  
function train = create_poisson(T,r,t_ref); 
	% function train = create_poisson(T,r,t_ref); 
	% This function is a homogeneous poisson spiker that can be used to
	% create one train of spikes based on the rate and parameters (T, r, t_ref) provided.
	
	%Initialize spikes, the spike time vector%%%%%%%%%%%%%%%%%%%
	spikes(1,1) = 0;
	k=1;
	
	%Create a series of spikes for the rate given %%%%%%%%%%%%%%%%%%%%
	while spikes(k,1) < T;
		spikes(k+1,1) = spikes(k,1)-(log(rand))/r; % This line introduces the exponentially distributed ISIs with average rate r.
	
		if spikes(k+1,1) - spikes(k,1) >= t_ref % Throw out the spike if it's isi < refractory period.
			k=k+1;
		else 
		end
	
	end
	spikes(k,1) = 0; % remove the spike that is after the Time period has ended
	
	%Create a vector (train) that holds all necessary spike time information    
	train(:,1) = nonzeros(spikes(2:end-1)); % No real spike at t = 0 or after T
end