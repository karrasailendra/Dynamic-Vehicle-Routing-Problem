clc;
clear;
close all;
% prompt = 'Enter the number of passengers = ';
P = 30; % input(prompt); 
% prompt = 'Enter the number of cars = ';
C = 7; % input(prompt);  
pos=100*rand(P*2,2);
number_of_parents=2;
number_of_children=4;
opttourlength=[];
itern = 200;
tourlengthmat=[];
for i=1:number_of_parents
parent(i).mat=zeros(C,P);
end
l=0;

 for i=1:P
     for l = 1:number_of_parents
     k=randi(C);
     parent(l).mat(k,i)=1;
     end
 end

for iter = 1 : itern 
 s1=sum(parent(1).mat');
 s2=sum(parent(2).mat');

out = vertisplit(parent(1).mat,parent(2).mat);
child(1).mat=out.child1;
child(2).mat=out.child2;
child(3).mat=parent(1).mat;
child(4).mat=parent(2).mat;

for i= 1:number_of_children % check for second iteration for child 3 and 4
    while max(sum(child(i).mat'))>5
         child(i).mat =  correction(child(i).mat);
         l=l+1;
%          disp(max(sum(child(i).mat')))
    end
end

for i = 1:number_of_children
    child(i).tour = tournum(child(i).mat,C,P);
    for k = 1: C
    temp = routetoxy(cell2mat(struct2cell(child(i).tour(k))),pos,P);
%     child(i).xy(:,1,k)=temp(:,1);
%     child(i).xy(:,2,k)=temp(:,2);
    opttourlength(k)= Antcolony(temp);
    end
    child(i).tourlength=sum(opttourlength);
end
for o = 1:number_of_children
tourlengthmat=[tourlengthmat child(o).tourlength];
end

inctour=sort(tourlengthmat);
for p = 1:number_of_children
if inctour(1)==tourlengthmat(p)
    nextchild(1)=p;
end
if inctour(2)==tourlengthmat(p)
    nextchild(2)=p;
end
end

tempchild(1).mat=child(nextchild(1)).mat;
tempchild(2).mat=child(nextchild(2)).mat;
parent(1).mat=tempchild(1).mat;
parent(2).mat=tempchild(2).mat;

besttourlength(iter)=inctour(1);
disp(iter)
end
%%

x=1:iter;
plot(x,besttourlength)
%%


%%  calculating route for ACO
function col = tournum(child,C,P)
    idx = find(child == 1);
    idx = idx';

    for i=1:C
    car(i).tour = [];
    end
    for l=1:P
    r = rem(idx(l),C);
    if r==0
        r=C;
    end
    car(r).tour= [car(r).tour idx(l)];
%     car(r,length(car(r,:))+1)=idx(l);
    end
    for k = 1:C
        [row,col(k).tour]=ind2sub(size(child),car(k).tour);
    end
  
 end
%%
function out = correction(child)
    s1=sum(child');
    maxk=max(s1);
    maxpass=find(s1 == maxk);
    idn=find(child(maxpass(1),:) == 1); % put randomperm for multiple maximum 
    pk=randperm(length(idn));
    child(maxpass(1),idn(pk(1)))=0;
    mink=min(s1);
    minpass=find(s1 == mink);
    minp=randperm(length(minpass));
    child(minpass(minp(1)), idn(pk(1)))=1;
    out = child;
end

function out = vertisplit(parent1,parent2)
    S=size(parent1);
    if mod(S(2),2) == 0
       k = S(2)/2;
    else
        k = (S(2)-1)/2;
    end
    child1=zeros(S);
    child2=zeros(S);
    child1(:,1:k)=parent1(:,1:k);
    child1(:,k+1:S(2))= parent2(:,k+1:S(2));
    child2(:,1:k)=parent2(:,1:k);
    child2(:,k+1:S(2))= parent1(:,k+1:S(2));
    out.child1=child1;
    out.child2=child2;
end
%% ACO function

function out = Antcolony(pos)
x = pos(:,1);
y = pos(:,2);
n = numel(x);
N=n;
    for i=1:n-1
        for j=i+1:n  
            D(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            D(j,i)=D(i,j);
        end
    end
    val.n=n;
    val.x=x;
    val.y=y;
    val.D=D;

Cost=@(tour) TourLength(tour,val);

nVar=val.n;
% conditions
it = 5;      
Ant  = val.n;      % check for number of `ants
Q     = 100;
tau0=10^-6;	
alpha=1;        
beta=5;         
rho=0.5;       
e = 5;


% Initialization
eta=1./val.D;          
tau=tau0*ones(nVar,nVar); 
Bestpath=zeros(it,1);    
empty_ant.Tour=[];
empty_ant.Cost=[];
ant=repmat(empty_ant,Ant,1);
BestSol.Cost=inf;


% ACO Main Loop

for it=1:it
    
    % Move Ants
    for k=1:Ant        
        % pool for requests 
        nonpool = N+1:2*N;
        ant(k).Tour=randi([1 N]); 
        pool = [1:N ant(k).Tour(end)+N];
        nonpool(nonpool==ant(k).Tour(end)+N)=[];
        for l=2:2*N           
            i=ant(k).Tour(end);            
            P=tau(i,:).^alpha.*eta(i,:).^beta;            
            P(ant(k).Tour)=0; 
            P(nonpool)=0;
            P=P/sum(P);           
            j=RouletteWheel(P);          
            ant(k).Tour=[ant(k).Tour j];
            if j<=N
            pool = [pool j+N]; 
            end
            nonpool(nonpool==j+N)=[];
        end       
        ant(k).Cost=Cost(ant(k).Tour);       
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
        end     
    end
    %disp(BestSol)
    tau=(1-rho)*tau;
    
    %  Pheromone Update
    for k=1:Ant
        tour=ant(k).Tour;
        bst = BestSol.Tour;
        tour=[tour tour(1)]; 
        bst=[bst bst(1)];
        for l=1:nVar
            i=tour(l);
            j=tour(l+1);
            ind = find(bst==i,1);
            if (j==bst(ind+1))
                tau(i,j)=tau(i,j)+ Q/ant(k).Cost + (e*Q/BestSol.Cost);
            else
                tau(i,j)=tau(i,j)+Q/ant(k).Cost;
            end
        end
    end
    

    out=BestSol.Cost; 
%     Bestpath(it)=BestSol.Cost;
%     disp(['Iteration ' num2str(it) ': Best path length = ' num2str(Bestpath(it))]);
%     figure(1);
%      tour=[BestSol.Tour BestSol.Tour(1)];
%     plot(val.x(tour),val.y(tour),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
%     
%     xlabel('x');
%     ylabel('y');
%     grid on;
%     drawnow;
    
end
end

%%
function L=TourLength(tour,model)

    n=numel(tour);

    tour=[tour tour(1)];
    
    L=0;
    for i=1:n
        L=L+model.D(tour(i),tour(i+1));
    end

end
function j=RouletteWheel(P)

    r=rand;
    
    C=cumsum(P);
    
    j=find(r<=C,1,'first');

end
function out=routetoxy(route,pos,P)
k =route;
for i= 1:length(k)
    xp(i)=pos(k(i),1);
    yp(i)=pos(k(i),2);
    xd(i)=pos(k(i)+P,1);
    yd(i)=pos(k(i)+P,2);
end

x = [xp xd];
y = [yp yd];
out(:,1)=x';
out(:,2)=y';
end