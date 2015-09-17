function [speech1, speech2, xest1, xest2] = demix_scatt2(mix, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, epsf, filts, options, param1, param2, Npad)

	%%% we will use a simpler algo: 
	%%% 1. compute |W1 x|
	[S2, S1, P] = audioscatt_fwd_haar(pad_mirror(mix',Npad), filts, options);
	
	%%% 2. estimate D1iz1i 
	if 1
		S1r = renorm_spect_data(S1,stds1,epsf);
	end
        H =  full(mexLasso(S1r,[Dnmf11,Dnmf12],param1));
        rec11 = Dnmf11*H(1:size(Dnmf11,2),:);
        rec12 = Dnmf12*H(size(Dnmf11,2)+1:end,:);
	
	%%% 3. use mask to obtain Ui = |W1 hat{x_i} | 
	eps = 1e-6;
        V_ap = rec11.^2 +rec12.^2 + eps;
        U1 = ((rec11.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
        U2 = ((rec12.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
		
	speech1 = audioreconstruct1(U1, options, filts{1}, P);
	speech2 = audioreconstruct1(U2, options, filts{1}, P);

	xest1 = speech1;
	xest2 = speech2;

	niters=getoptions(options,'demix_iters',1);
	for l=1:niters

	%%% 4. compute new Ui by minimizing || |W2 Ui | - D2i z2i || and using the phase of each Ui to go down.
	[S21, S11, P1, P1lp, P1hp, scratch1] = audioscatt_fwd_haar(pad_mirror(speech1,Npad), filts, options);
	[S22, S12, P2, P2lp, P2hp, scratch2] = audioscatt_fwd_haar(pad_mirror(speech2,Npad), filts, options);

	if 1
		[S21r,norm1] = renorm_spect_data(S21, stds2,epsf);
		[S22r,norm2] = renorm_spect_data(S22, stds2,epsf);
	end
	
	H1 = full(mexLasso(S21r,Dnmf21,param2));
	Srec21 = Dnmf21*H1;
	H2 = full(mexLasso(S22r,Dnmf22,param2));
	Srec22 = Dnmf22*H2;
	
	if 1
		Srec21u = unrenorm_spect_data(Srec21,stds2,norm1);
		Srec22u = unrenorm_spect_data(Srec22,stds2,norm2);
	end

	rec21 = audioreconstruct2(Srec21u, P1lp, P1hp, filts, scratch1);
	rec22 = audioreconstruct2(Srec22u, P2lp, P2hp, filts, scratch2);

	%%%5. denoise in the first level dictionary
	if 1
	rec1r = renorm_spect_data(rec21,stds1);
	rec2r = renorm_spect_data(rec22,stds1);
	H1 = full(mexLasso(rec1r,Dnmf11, param1));
	H2 = full(mexLasso(rec2r,Dnmf12, param1));
	rec21 = Dnmf11*H1;
	rec22 = Dnmf12*H2;
	end
	%%% 6. reapply the mask to fullfil constraints and output xi.
	eps = 1e-6;
        V_ap = rec21.^2 +rec22.^2 + eps;
        U1b = ((rec21.^2)./(V_ap)).*S1(:,1:size(V_ap,2));
        U2b = ((rec22.^2)./(V_ap)).*S1(:,1:size(V_ap,2));

		
	speech1 = audioreconstruct1(U1b, options, filts{1}, P);
	speech2 = audioreconstruct1(U2b, options, filts{1}, P);

        %%%%%%%%%%%%%%%%
	end




