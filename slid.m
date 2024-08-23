% Sliders

function [SL1 SL2 fig] = slid( txt, Pos, nh, RN )
%[SL1 SL2 fig] = slid('h-counter & root counter for i',[1300 70 600 130],nh,RN); 
Pos = [1300 70 600 130];
fig = uifigure('Name','h-counter & root counter for i','Position',Pos);
%pn  = uipanel(fig);
SL1 = uislider(fig,'Position',[20 100 300 3]);
SL1.Limits = [1 nh];
SL2 = uislider(fig,'Position',[20 45 550 3]);
SL2.Limits = [0 RN];
fig = uifigure('Name',txt,'Position',Pos);
%pn  = uipanel(fig);
SL1 = uislider(fig,'Position',[20 100 300 3]);
SL1.Limits = [1 nh];
SL2 = uislider(fig,'Position',[20 45 550 3]);
SL2.Limits = [0 RN];