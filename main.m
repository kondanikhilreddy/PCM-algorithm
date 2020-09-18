% Possibilistic C-means algo
clear,clc,close all;
tic;                                         %Auto time count and display
data = load('ex7data2.mat');    %Data load
data= data.X;

%data(:,1) = data(:,1) ./(max(data(:,1)) - min(data(:,1))) ; %normalization
%data(:,2) = data(:,2) ./(max(data(:,2)) - min(data(:,2))) ; % if required

%{
%Random data generation
data1 = randn(200,1)+3; % Initialized the Data Variables
data2 = randn(200,1)+7;

data(1:100,1:2) = [data1(1:100),data2(1:100)];
data(101:200,1:2) = [data2(101:200),data1(101:200)];
%}

[centers,U] = pcm(data,3);

%plot
index1 = find(U(1,:) == max(U));
index2 = find(U(2, :) == max(U));
index3 = find(U(3, :) == max(U));
line(data(index1,1), data(index1,2), 'lineStyle', 'none','marker', 'o','color', 'g');
line(data(index2,1), data(index2,2), 'LineStyle', 'none','Marker', 'o','Color', 'r');
line(data(index3,1), data(index3,2), 'LineStyle', 'none','Marker', 'o','Color', 'b');
line(centers(:,1),centers(:,2),'lineStyle','none','markersize',21,'marker', '.','color', 'k');
toc

function [centers,U_newest] = pcm(dataset,clusters)
data_points = size(dataset,1);
distances=ones(data_points,clusters);

% membership matrix  is initialized
U = randn(clusters,data_points) + 3;
val = sum(U);
U = U./val;

%parameters
m = 2;                     %fuzzifier
a= 1/(m-1);              %power raised 
K=0.01;                   %Factor for intial neta 
iter = 100;                %max iterations
epsilon = 0.0001;      %condtion to break
alpha=0.2;               % alpha cut 

neta= ones(1,clusters);
mf = U.^m;
centers =[0,0; max(dataset(:,1)),0; 0,max(dataset(:,2))];
 for k=1:clusters
     neta(1,k) = (mf(k,:) * (pdist2(dataset,centers(k,:)) .^2));
 end
 neta = neta./(sum(mf,2))';
 neta = ((ones(data_points,1) * neta) .^a)*K ; %neta for all clusters

 %Updating of memebership matrix
 for i = 1:iter    
    for k=1:clusters
        distances(:,k) =  pdist2(dataset,centers(k,:)) .^(2*a) ;
    end
    den = neta + distances;
    U_new = (neta ./den)';
    mf = U_new.^m;
    denominator = ones(1,2) .* sum(mf,2) ;
   centers = (mf * dataset) ./ denominator;    %center update
    if max(max(abs(U_new - U))) < epsilon  % break condition
        break;
    end
    U = U_new;
end
%Repeating the updating with right estimate of neta
d=zeros(1,clusters);
neta= ones(1,clusters);
for b = 1:clusters                    
    for c= 1: data_points
        if U_new(b,c) >=alpha  %alpha cut
            mf(b,c)=1;
            d(1,b)= d(1,b)+1;
        end
       if U_new(b,c) < alpha    %alpha cut
           mf(b,c)=0;
       end
    end
end
 for k=1:clusters
     neta(1,k) = (mf(k,:) * (pdist2(dataset,centers(k,:)) .^2));
 end
 neta = neta./d;
 neta = (ones(data_points,1) * neta) .^a; %The appropriate neta

 %Updating the membership matrix again
for i = 1:iter     
    for k=1:clusters
        distances(:,k) =  pdist2(dataset,centers(k,:)) .^(2*a) ;
    end
    den = neta + distances;
    U_newest = (neta ./den)';
    mf = U_newest.^m;
    denominator = ones(1,2)  .* sum(mf,2) ;
    centers = (mf * dataset) ./ denominator;                %center update
    if max(max(abs(U_newest - U_new))) < epsilon  % break condition 
        break;
    end
    U_new = U_newest;
end
end

