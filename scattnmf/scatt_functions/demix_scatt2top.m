function [speech1, speech2, xest1, xest2] = demix_scatt2top(mix, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, epsf, filts, options, param1, param2, Npad)

	%%% we will use a simpler algo: 
	%%% 1. compute |W1 x|
	[S2, S1, P, Plp, Php, scratch] = audioscatt_fwd_haar(pad_mirror(mix',Npad), filts, options);
	
		size(S2)
		[S2r,norm1] = renorm_spect_data(S2, stds2, epsf);
	
	H = full(mexLasso(S2r,[Dnmf21,Dnmf22],param2));
	Srec21 = Dnmf21*H(1:size(Dnmf21,2),:);
	Srec22 = Dnmf22*H(size(Dnmf21,2)+1:end,:);
	
		Srec21u = unrenorm_spect_data(Srec21,stds2,norm1);
		Srec22u = unrenorm_spect_data(Srec22,stds2,norm1);

	rec21 = audioreconstruct2(Srec21u, Plp, Php, filts, scratch);
	rec22 = audioreconstruct2(Srec22u, Plp, Php, filts, scratch);

	%%% 3. use mask to obtain Ui = |W1 hat{x_i} | 
	eps = 1e-6;
        V_ap = rec21.^2 +rec22.^2 + eps;
        U1 = ((rec21.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
        U2 = ((rec22.^2)./(V_ap)).*S1(:,1:size(V_ap,2));

	%%% 2. estimate D1iz1i 
		S1r = renorm_spect_data(U1,stds1,epsf);
		S2r = renorm_spect_data(U2,stds1,epsf);
        H1=  full(mexLasso(S1r,Dnmf11,param1));
        H2=  full(mexLasso(S2r,Dnmf12,param1));
        rec11 = Dnmf11*H1;
        rec12 = Dnmf12*H2;

	%%% 3. use mask to obtain Ui = |W1 hat{x_i} | 
	eps = 1e-6;
        V_ap = rec11.^2 +rec12.^2 + eps;
        U1 = ((rec11.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
        U2 = ((rec12.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
		
	speech1 = audioreconstruct1(U1, options, filts{1}, P);
	speech2 = audioreconstruct1(U2, options, filts{1}, P);


	%%%%first level estimator

		S1r = renorm_spect_data(S1,stds1, epsf);
        H =  full(mexLasso(S1r,[Dnmf11,Dnmf12],param1));
        rec11 = Dnmf11*H(1:size(Dnmf11,2),:);
        rec12 = Dnmf12*H(size(Dnmf11,2)+1:end,:);
	
	eps = 1e-6;
        V_ap = rec11.^2 +rec12.^2 + eps;
        U1 = ((rec11.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
        U2 = ((rec12.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
		
	xest1 = audioreconstruct1(U1, options, filts{1}, P);
	xest2 = audioreconstruct1(U2, options, filts{1}, P);

	speech1 = speech1';
	speech2 = speech2';
	xest1 = xest1';
	xest2 = xest2';

	speech1 = speech1(1:length(mix));
	speech2 = speech2(1:length(mix));
	xest1 = xest1(1:length(mix));
	xest2 = xest2(1:length(mix));
