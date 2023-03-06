clear all

N = 100;
c = 30;
cp = 10;
cn = 3;
Tmax = 10;
num = 50;
movie = 1;
%% Simulation

[p_list, R_list] = triadic_percolation_simulation_poisson_orbit(N, c, cp, cn, Tmax, num);
scatter(p_list, R_list)


%% Theory

[p_list, R_list] = triadic_percolation_theory_poisson_orbit(c, cp, cn, Tmax, num);
scatter(p_list, R_list)

%% Movie: Illustrate the time series of the dynamic in a 2-d space.

movie_simulation(N, c, cp, cn, Tmax)

%%

function [p_list, R_list] = triadic_percolation_theory_poisson_orbit(c, cp, cn, Tmax, num)
    % Theoretical calculation of the orbit diagram of the dynamic

    % c: Average structural degree of the Poisson network
    % cp: Average positive regulatory degree of the Poisson regulatory network
    % cn: Average negative regulatory degree of the Poisson regulatory network
    % Tmax: max duration of the dynamic
    % num: Number of sample used to plot the orbit diagram

    R_list = [];
    p_list = [];

   

    for p = 0:0.01:1
        SN2 = 0.9;
        pL = 0.11;   % Initial condition

        disp(p)
        for i = 1:Tmax
            
            for nrun = 1:100
                SN2 = (1 - exp(-c * pL * SN2));
            end
            pL = p * abs(exp(-cn * SN2) * (1 - exp(-cp * SN2)));
            
            R(i) = SN2;
            
        end
        R_list = [R_list; R(end-num+1:end)']; 
        p_list = [p_list; p * ones(num,1)];
    end


end


function [p_list, R_list] = triadic_percolation_simulation_poisson_orbit(N, c, cp, cn, Tmax, num)

    % Monte Carlo simulation of the orbit diagram of the dynamic

    % N: Number of nodes
    % c: Average structural degree of the Poisson network
    % cp: Average positive regulatory degree of the Poisson regulatory network
    % cn: Average negative regulatory degree of the Poisson regulatory network
    % Tmax: max duration of the dynamic
    % num: Number of sample used to plot the orbit diagram

    % Adjacency matrix of a Poisson network

    a=rand(N,N);
    a=a<(c/(N-1));
    a=triu(a,1);
    a=a+a';
    [I,J,V]=find(triu(a));

    % Regulatory network

    L = sum(sum(a)) / 2;    % Number of links
    rand_reg = rand(N, L);
    
    adj_pos = rand_reg < cp/N;
    adj_neg =(rand_reg > cp/N) .* (rand_reg < (cp+cn)/N);
    

    % The dynamic
    p_list = [];
    R_list = [];    % Store p and corresponding R

    for p = 0:0.02:1
        disp(p);

        pL0 = 0.11;  % Initial conidition

        xL=rand(L,1);
        state=(xL<pL0);  % initial states of link
        [linkID, ~, ~] = find(state);   
        adj = sparse(I(linkID), J(linkID), V(linkID), N, N); % Generate a new network with active links
        adj = adj + adj';
        G=graph(adj);
    
        [bin,binsizes]=conncomp(G); 
        % bin[i]: id of component that node i belongs to
        % binsizes[i]: size of component i
    
        s = zeros(N,1);  % State of the nodes

        if max(binsizes) > 1    % Giant component
            [~,J0,~]=find(binsizes==max(binsizes));
            id=J0(1);   % id of the giant component
            s=(bin==id)';
        end

        R(1)=max(binsizes) / N;
        for it=2:Tmax
   
            xL=rand(L,1);
            state=(xL<p).*((adj_pos'*s)>0).*((adj_neg'*s)==0);  % state of link
            
            [linkID, ~, ~] = find(state);   
            adj = sparse(I(linkID), J(linkID), V(linkID), N, N); % Generate a new network with active links
            adj = adj + adj';
            G = graph(adj);
            [bin,binsizes]=conncomp(G);
            s = zeros(N,1);
            if max(binsizes) > 1    % Giant component: sqrt(N)?
                [~,J0,~]=find(binsizes==max(binsizes));
                id=J0(1);   % id of the giant component
                s=(bin==id)';
            end
            R(it)=max(binsizes) / N;


        end
        p_list = [p_list; p * ones(num, 1)];
        R_list = [R_list; R(end-num+1:end)'];

    end

end


function movie_simulation(N, c, cp, cn, Tmax)
    
    % Produce a movie of the time series of the dynamic.

    % N: Number of nodes
    % c: Average structural degree of the Poisson network
    % cp: Average positive regulatory degree of the Poisson regulatory network
    % cn: Average negative regulatory degree of the Poisson regulatory network
    % Tmax: max duration of the dynamic
   
    a=rand(N,N);
    a=a<(c/(N-1));
    a=triu(a,1);
    a=a+a';
    [I,J,V]=find(triu(a));

    x = rand(N,1);      % x-coordinate of nodes
    y = rand(N,1);      % y-coordinate of nodes. The coordinates are only for illustration purpose.

    % Regulatory network

    L = sum(sum(a)) / 2;    % Number of links
    rand_reg = rand(N, L);
    
    adj_pos = rand_reg < cp/N;
    adj_neg =(rand_reg > cp/N) .* (rand_reg < (cp+cn)/N);

    pL0 = 0.11;  % Initial conidition
    p = 0.8; 
    xL=rand(L,1);
    state=(xL<pL0);  % initial states of link
    [linkID, ~, ~] = find(state);   
    adj = sparse(I(linkID), J(linkID), V(linkID), N, N); % Generate a new network with active links
    adj = adj + adj';
    G=graph(adj);

    [bin,binsizes]=conncomp(G); 
    % bin[i]: id of component that node i belongs to
    % binsizes[i]: size of component i

    s = zeros(N,1);  % State of the nodes

    if max(binsizes) > 1    % Giant component
        [~,J0,~]=find(binsizes==max(binsizes));
        id=J0(1);   % id of the giant component
        s=(bin==id)';
    end
    R(1)=max(binsizes) / N;

    for it=2:Tmax
        xL=rand(L,1);
        state=(xL<p).*((adj_pos'*s)>0).*((adj_neg'*s)==0);  % state of link
        
        [linkID, ~, ~] = find(state);   
        adj = sparse(I(linkID), J(linkID), V(linkID), N, N); % Generate a new network with active links
        adj = adj + adj';
        G = graph(adj);
        [bin,binsizes]=conncomp(G);
        s = zeros(N,1);
        if max(binsizes) > 1    % Giant component: sqrt(N)?
            [~,J0,~]=find(binsizes==max(binsizes));
            id=J0(1);   % id of the giant component
            s=(bin==id)';
        end
        R(it)=max(binsizes) / N;
        [Is, Js, Vs] = find(bin == id);
        [IsN,JsN,Vs]=find(s'==0);
        
        hold off
        plot(x(Js),y(Js),'o','MarkerFaceColor','#77AC30', 'MarkerSize',12)
        hold on
       
        ChildrenBefore=get(gca,'children');
        gplot(adj, [x y])
      
        ChildrenAfter=get(gca,'children');
        NewChildren=setdiff(ChildrenAfter,ChildrenBefore);
        set(intersect(findall(gcf,'type','line'),NewChildren),'LineWidth',2,'Color','#4DBEEE')
       
        plot(x(Js),y(Js),'o','MarkerFaceColor','#77AC30', 'MarkerSize',12)
    
        plot(x(JsN),y(JsN),'o','Color','k', 'MarkerSize',12)
        axis off
        drawnow;
        pause(0.4);
   end


end