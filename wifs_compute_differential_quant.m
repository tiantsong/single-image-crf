function wifs_compute_differential_quant(fileinfo)

% ----------- constant -------------

isLimitQ = true;

% ------- vieville bilateral ----------

if ischar(fileinfo.imname)
    R_image = imread(fileinfo.imname);
    R_image = im2double(R_image(:,:,fileinfo.color_id));
else
    R_image = im2double(imname(:,:,fileinfo.color_id));
end


outputOrder = 2;
polyOrder = 4;
half_winN = 4;

tic
[vRo, vRx, vRy, vRxx, vRyy, vRxy] = derivatives_xy(R_image,'vieville',outputOrder,polyOrder);
toc

clear R_image

compute_differential_quantities(fileinfo, vRo, vRx, vRy, vRxx, vRyy, vRxy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function compute_differential_quantities(fileinfo,Ro, Rx, Ry, Rxx, Ryy, Rxy)

LR = Rxx + Ryy;
LRm = Rxx - Ryy;
Rx2 = Rx.^2;
Ry2 = Ry.^2;

% ------- invariant -------

den = (Rx2 + Ry2).^2;
den(den == 0) = eps;

DI1 = (Rx2.*(LR - 2*Rxy) + Ry2.*(LR + 2*Rxy) + 2*Rx.*Ry.*LRm)./den;
DI2 = ((Rx2 - Ry2).*LRm + 4*Rx.*Ry.*Rxy)./den;
DI3 = (Rx2.*(LR + 2*Rxy) + Ry2.*(LR - 2*Rxy) - 2*Rx.*Ry.*LRm)./den;

save(fileinfo.exdataname,'Rx','Ry','DI1','DI2','DI3');

clear LR LRm Rx2 Ry2 den

% ------- Q ----------

Q = (DI1 + DI2 + DI3)/3;
Q = 1./(1-Q.*Ro);

% --------------------

Qbar = exp_5_calibrate_Q(Q,Ro,'gammacali','constant');
Qbar = max(min(Qbar,1),0);

SM = sqrt(((DI3 - DI2)/0.9).^2 + (DI1-DI3).^2 + ((DI2-DI1)/0.9).^2);

stdQ = 1./(1+7.5*Qbar);
SM = SM.*(Ro.^0.7)./stdQ;



clear Qbar stdQ

% -------- normal grad by magnitude ----------

grad_mag = sqrt(Rx.^2 + Ry.^2);
grad_mag(grad_mag == 0) = 1;

Cx = Rx./grad_mag;
Cy = Ry./grad_mag;

clear grad_mag;


% -------- Cave -------

nL = 3;
H = fspecial('gaussian',2*nL+1,1);

C_avex = imfilter(Cx,H,'symmetric','same');
C_avey = imfilter(Cy,H,'symmetric','same');
C_ave = 1 - sqrt(C_avex.^2 + C_avey.^2);

clear C_avex C_avey



% ====== compute weight ========

ylimit = [0,1.6,40];
ylimit_rq = [0,1.8,40];
xlimit = [0,1,80];



BR.Q = Q;
BR.Ro = Ro;
BR.C_ave = C_ave;
BR.SM = SM/4.3;

clear Q Ro C_ave SM

BR.SM = abs(BR.SM);

% -------- calibrate Q -------

BR.Q = exp_5_calibrate_Q(BR.Q,BR.Ro,fileinfo.Qtype,fileinfo.Qorder);

% -------- selection ----------------

H = fspecial('gaussian',7,0.75);
BR.SMN = imfilter(BR.SM,H,'symmetric','same');

% -------- weight ------------

sm_std = 1;
C_std = 0.15/2;

W = gaussian_function(BR.SMN,sm_std).*gaussian_function(BR.C_ave,C_std);

W_threshold = gaussian_function(2*sm_std,sm_std)*gaussian_function(2*C_std,C_std);

I = W > W_threshold & ...
    BR.Q > ylimit(1) & BR.Q < ylimit(2);

% -------- remove selection at border -----------

nL = 20;
I([1:nL end-nL+1:end],:) = false;
I(:,[1:nL end-nL+1:end]) = false;

fprintf('#I = %d\n',sum(sum(I)));

R = BR.Ro;
Q = BR.Q;

save(fileinfo.dataname,'R','Q','W','I');