% November 4, 2020
% From: https://courses.washington.edu/matlab1/ModelFitting.html#3

%% One parameter example: Weber's Law 

% Weber's law states that the ability for a subject to notice an increase 
% in stimulus intensity is proportional to the starting, or baseline 
% intensity. That is, if x is the stimulus intensity, the increment threshold
% is kx, where k is the 'Weber fraction'. This fraction is our one parameter.


% baseline weights
x = [.5,1,1.5,2,2.5,3];

% increment thresholds
y = [  0.0619    0.0888    0.1564    0.1940    0.2555    0.2890];

% plot the raw data
figure(1)
clf
h1=plot(x,y,'ko','MarkerFaceColor','k');
set(gca,'YLim',[0,0.35]);
set(gca,'XLim',[0,3.5]);
set(gca,'XTick',x)
xlabel('Baseline weight (Kg)');
ylabel('Increment threshold (Kg)');


% write a function containing the model that will predict the data
% this is contained in the folder "scratch"
pred = WebersLaw(p,x);

% plot the prediction
hold on
h2=plot(x,pred,'r-');
legend([h1,h2],{'Data','First Guess'},'Location','NorthWest');

%% now compare the model prediction to the data

%The SSE between our model and the data for our first guess is:

err = WebersLawErr(p,x,y);





