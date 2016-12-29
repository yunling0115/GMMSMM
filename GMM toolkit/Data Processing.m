% Data Processing
clear all;
close all;
load DSHOCKS.DAT;
load DERETS.DAT;
load DEMZS.DAT;
% name all variables
% (1) DERETS
yyyy = DERETS(:,1);
rm_rf = DERETS(:,2);
smb = DERETS(:,3);
hml = DERETS(:,4);
momwml = DERETS(:,5);
revlmw = DERETS(:,6);
junk_cb = DERETS(:,7);
gb_rf = DERETS(:,8);
roaspd = DERETS(:,9);
invspd = DERETS(:,10);
% (2) DEMZS
rf = DEMZS(:,2);
dp = DEMZS(:,3);
% (3) DSHOCKS
cgrowth_shock = DSHOCKS(:,2);
expgrowth_shcok = DSHOCKS(:,3);
vol_shock = DSHOCKS(:,4);
% (4) All Factors
AF = [DERETS(:,2:10), DEMZS(:,2:3), DSHOCKS(:,2:4)];
fprintf('Data Processing Done \n ');

