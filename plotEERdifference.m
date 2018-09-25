% intbits = [4,8,12,16,20,24,28,32,4,8,12,16,20,24,28,32,4,8,12,16,20,24,28,32, ... 
%     4,8,12,16,20,24,28,32,4,8,12,16,20,24,28,32,4,8,12,16,20,24,28,32, ...
%     4,8,12,16,20,24,28,32,2,3,4,5,6,7,8,12,16,20,24,28,32, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7, ...
%     2,3,5,6,7];
% fracbits = [4,4,4,4,4,4,4,4,8,8,8,8,8,8,8,8,12,12,12,12,12,12,12,12, ...
%     16,16,16,16,16,16,16,16,20,20,20,20,20,20,20,20,24,24,24,24,24,24,24,24, ...
%     28,28,28,28,28,28,28,28,32,32,32,32,32,32,32,32,32,32,32,32,32, ...
%     4,4,4,4,4, ...
%     8,8,8,8,8, ...
%     12,12,12,12,12, ...
%     16,16,16,16,16, ...
%     20,20,20,20,20, ...
%     24,24,24,24,24, ...
%     28,28,28,28,28];
% % EER value without reducing bits
% EER = 0.398; 
% 
% % EER values after reducing bits
% EERsafter = [0.4845,0.3995,0.4045,0.4105,0.4040,0.4075,0.4030,0.4065,0.4885,...
%     0.4050,0.4020,0.4025,0.4005,0.4025,0.4085,0.4055, 0.4970,0.4125, 0.4050,...
%     0.4035,0.4000,0.4030,0.4015,0.3995,0.5005,0.4040,0.3985,0.4060,0.3995, ...
%     0.3955,0.4030,0.4090,0.4905,0.3990,0.4075,0.4055,0.4055,0.4005,0.4055,...
%     0.4050,0.4920,0.3980,0.4110,0.4050,0.4050,0.4045,0.4110,0.4055,0.4855,...
%     0.4035,0.4040,0.4050, 0.4080,0.4060,0.3980,0.3990,0.5120,0.5030,0.4935,0.4725,0.4490,0.3995,0.4055,0.4035,...
%     0.3985,0.4030,0.4005,0.4000,0.4130, ...
%     0.5040,0.5030,0.4645,0.4550,0.4015, ...
%     0.5095,0.5010,0.4670,0.4535,0.4015, ...
%     0.5060,0.5030,0.4595,0.4580,0.4045, ...
%     0.5085,0.5140,0.4680,0.4515,0.4045, ...
%     0.5060,0.5100,0.4675,0.4620,0.3985, ...
%     0.5120,0.5040,0.4615,0.4645,0.4065, ...
%     0.5000,0.5030,0.4630,0.4600,0.4035];
% 
% EERdifferences = EERsafter-EER;



intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [4,4,4,4,4,4,4,4,4,4,4,4,4];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5040, 0.5030, 0.4845,0.4645 ,0.4550 , 0.4015,0.3995,0.4045,0.4105,0.4040,0.4075,0.4030,0.4065];

EERdifferences = EERsafter-EER;



grid on
% 
% pts = [EERs; features];
% 
% plot3(pts(:,1), pts(:,2), pts(:,3));

p = plot3(intbits, fracbits, EERdifferences, '-o');
hold on



intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [8,8,8,8,8,8,8,8,8,8,8,8,8];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5095,0.5010,0.4885,0.4670,0.4535,0.4015,0.4050,0.4020,0.4025,0.4005,0.4025,0.4085,0.4055];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');


intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [12,12,12,12,12,12,12,12,12,12,12,12,12];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5060,0.5030,0.4970,0.4595,0.4580,0.4045,0.4125, 0.4050,0.4035,0.4000,0.4030,0.4015,0.3995];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');


intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [16,16,16,16,16,16,16,16,16,16,16,16,16];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5085,0.5140,0.5005,0.4680,0.4515,0.4045,0.4040,0.3985,0.4060, 0.3995,0.3955,0.4030,0.4090];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');


intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [20,20,20,20,20,20,20,20,20,20,20,20,20];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5060,0.5100,0.4905,0.4675,0.4620,0.3985,0.3990,0.4075,0.4055,0.4055,0.4005,0.4055,0.4050];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');

intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [24,24,24,24,24,24,24,24,24,24,24,24,24];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5120,0.5040, 0.4920,0.4615, 0.4645,0.4065,0.3980,0.4110,0.4050,0.4050,0.4045,0.4110,0.4055];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');

intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [28,28,28,28,28,28,28,28,28,28,28,28,28];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5000,0.5030,  0.4855,0.4630,0.4600,0.4035,0.4035,0.4040,0.4050,0.4080,0.4060,0.3980,0.3990];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');


intbits = [2,3,4,5,6,7,8,12,16,20,24,28,32];
fracbits = [32,32,32,32,32,32,32,32,32,32,32,32,32];
% EER value without reducing bits
EER = 0.398; 

% EER values after reducing bits
EERsafter = [0.5120,0.5030, 0.4935,0.4725,0.4490, 0.3995,0.4055,0.4035,0.3985,0.4030,0.4005,0.4000,0.4130];

EERdifferences = EERsafter-EER;
p = plot3(intbits, fracbits, EERdifferences, '-o');


ax = gca;
ax.ZGrid = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
% set(gca,'yticklabel',[])
set(gca,'FontSize', 24)
% p.LineStyle = 'none'
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
% plot3(temp1, temp2, temp3, '-o');

% legend('Using doubles','Using integers');

% plot3([0.36, 0.24], [0.45,0.235], [0.31,0.12], '-o');
xlabel('integer bits');
ylabel('fraction bits');
zlabel('EER differences');