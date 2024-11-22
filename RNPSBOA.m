
function[Best_score,Best_pos,SBOA_curve]=RNPSBOA(SearchAgents_no,Max_iter,lowerbound,upperbound,dimension,fitness)

lowerbound=ones(1,dimension).*(lowerbound);                              
upperbound=ones(1,dimension).*(upperbound);     


%% INITIALIZATION
r = 0.7; 
X_chaos = rand(SearchAgents_no, dimension);
for i = 1:SearchAgents_no
    for j = 1:dimension
        if X_chaos(i, j) < 0.5
            X_chaos(i, j) = cos(pi * (r * sin(pi * X_chaos(i, j)) + 2 * (1 - r) * X_chaos(i, j) - 0.5));
        else
            X_chaos(i, j) = cos(pi * (r * sin(pi * X_chaos(i, j)) + 2 * (1 - r) * (1 - X_chaos(i, j)) - 0.5));
        end
    end
end
X = repmat(lowerbound, SearchAgents_no, 1) + X_chaos .* repmat((upperbound - lowerbound), SearchAgents_no, 1);
for i =1:SearchAgents_no
    M=X(i,:);
    fit(i)=fitness(M);
end
%% main loop
for t=1:Max_iter
    CF=(1-t/Max_iter)^(2*t/Max_iter);
    %% update the global best (fbest)
    [best , location]=min(fit);
    if t==1
        Bast_P=X(location,:);                                           
        fbest=best;                                          
    elseif best<fbest
        fbest=best;
        Bast_P=X(location,:);
    end
    
    %% The secretary bird's predation strategy
    for i=1:SearchAgents_no
        if t<Max_iter/3  % Secretary bird search prey stage
            Epsilon = 0.5;
            q = 2;
            selected_index_q = randperm(SearchAgents_no, q);
            Xq = X(selected_index_q, :);
            Xqmean = mean(Xq);
           A = randperm(SearchAgents_no);
            R1 = A(1);

            if rand < Epsilon
           Rn=size(X,1);
           X_random_1=randi([1,Rn]);
           X_random_2=randi([1,Rn]);
           R2=rand(1,dimension);
            X1=X(i,:)+(X(X_random_1,:)-X(X_random_2,:)).*R2; 
            
            else
            X1 = X(i, :) + (Xqmean - X(R1, :)) .* rand; % Eq. (4)
            
            end 
        elseif t>Max_iter/3 && t<2*Max_iter/3  % Secretary bird approaching prey stage
              
               RB=randn(1,dimension);
               X1=Bast_P+exp((t/Max_iter)^4)*(RB-0.5).*(Bast_P-X(i,:));
              
        else       % Secretary bird attacks prey stage
            q = 0.5;
            r = rand; 
            if r > q
            RL=0.5*Levy(dimension);
            X1=Bast_P+CF*X(i,:).*RL;
          
            else
            q = 2;
            selected_index_q = randperm(SearchAgents_no, q);
            Xq = X(selected_index_q, :);
            Xqmean = mean(Xq);  
             R=0.02*(1-t/Max_iter);% Eq.(6)
             X1= X(i,:)+ (-R+2*R*rand(1,dimension)).*(Xqmean - X(i,:));% Eq.(7)
           
            end
        end
         X1= max(X1,lowerbound);
            X1 = min(X1,upperbound);

        f_newP1 = fitness(X1);
        if f_newP1 <= fit (i)
            X(i,:) = X1;
            fit (i)=f_newP1;
        end
    end
     
%% Secretary Bird's escape strategy
    r=rand;
    k=randperm(SearchAgents_no,1);
    Xrandom=X(k,:);
    for i=1:SearchAgents_no
        if r<0.5
            %% C1: Secretary birds use their environment to hide from predators
            if rand >0.9
            RB=rand(1,dimension);
            X2= Bast_P+(1-t/Max_iter)*(1-t/Max_iter)*(2*RB-1).*X(i,:);% Eq.(5) S1
            
            else
            r = rand; 
            b = 4; 
           if r <= 0.4
              X2 = X(i,:) + (upperbound - X(i,:)) * r * (1 - t/Max_iter)^b;
              
           elseif r > 0.4 && r <= 0.8
              X2 = X(i,:) - (X(i,:) - lowerbound) * r * (1 - t/Max_iter)^b;
             
           elseif r > 0.8 && r <= 0.9
              X2 = X(i,:) + (upperbound - X(i,:)) * r * (t/Max_iter)^b;
              
           else
              X2 = X(i,:) - (X(i,:) - lowerbound) * r * (t/Max_iter)^b;
              
           end
            end
            
      
        else
            %% C2:  Secretary birds fly or run away from the predator
            if rand>0.5
            K=round(1+rand(1,1));
            R2=rand(1,dimension);
            X2=X(i,:)+ R2.*(Xrandom-K.* X(i,:)); %  Eq(5) S2
            
            else
            Epsilon = unifrnd(0, 1);
            F=0.5;
            rand_leader_index1 = floor(SearchAgents_no * rand() + 1);
            rand_leader_index2 = floor(SearchAgents_no * rand() + 1);
            X_rand1 = X(rand_leader_index1, :);
            X_rand2 = X(rand_leader_index2, :);
            step2 = X_rand1 - X_rand2;
            if rand<0.5
                X2 = X(i, :) + Epsilon .* step2; % Eq.(11)4.3
              
            else
                X2 = X(i, :) + F .* Levy(dimension) .* step2;
                
            end
            end
        end
         X2= max(X2,lowerbound);X2 = min(X2,upperbound);
         f_newP2 = fitness(X2); %Eq (6)
        if f_newP2 <= fit (i)
            X(i,:) = X2;
            fit (i)=f_newP2;
        end

    end %

    best_so_far(t)=fbest;
    average(t) = mean (fit);

end 
Best_score=fbest;
Best_pos=Bast_P;
SBOA_curve=best_so_far;
end


