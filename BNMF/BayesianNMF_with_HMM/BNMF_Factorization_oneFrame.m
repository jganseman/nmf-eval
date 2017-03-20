function [hMixed,MagOut]=BNMF_Factorization_oneFrame(hMixed,Y,n,snr_in, spec_scale,init_noise_data)
% This function is the main processing routine
% where the MMSE estimates for speech are computed
% (in a frame-by-frame routine and the noise dictionary is updated if required )
%
% Nasser Mohammadiha Nov 2012
%
%         reference:
%         [1] N. Mohammadiha, P. Smaragdis and A. Leijon, “Supervised and
%         Unsupervised Speech Enhancement Using Nonnegative Matrix
%         Factorization,” IEEE Trans. Audio, Speech, and Language Pro-
%         cess., vol. 21, no. 10, pp. 2140–2151, oct. 2013.


%persistent variable are used to implement the sliding window based noise
%basis update since their value shoud be preserved from each frame to
%the next one
persistent noise_data local_buffer_cashed newNoiseInBuffer


size_local_buffer=15;%size of buffer m in Section III.B in the paper. required for online learning
update_length_noise_data=5;%q in Section III.B in the paper.required for online learning

nmf_V=hMixed.OutputDistr;
UserData=hMixed.UserData;
noise_name=UserData{2};

%SNR-alpha curve. obtained emprically, see fig. 5 in the paper.
%interpolate between two closest neighbours to find smoothing factor
if strcmp(noise_name,'online')
    snrss=[-7.5 -5 -2.5 0 2 4 5 6.5 8 10 12.5 15];
    alpha=[.99 .99 .97 	.95 .9 .8 .7 .6 .5 .5 .4 .2];
else
    snrss=[-10:2.5:5 6 7.5 9 10];
    alpha=[.995 .99 .99 .97 .95 .9 .7 .6 .4 .2 0];
end

if snr_in>=snrss(end)
    alpha_smoothing=alpha(end);
elseif snr_in<=snrss(1)
    alpha_smoothing=alpha(1);
else
    e=abs(snrss-snr_in);
    [es,ind_es]=sort(e,'ascend');
    if es(1)==0
        alpha_smoothing=alpha(ind_es(1));
    elseif snr_in<snrss(ind_es(1))
        weights(1)=1/(es(1)+eps);
        alphaa(1)=alpha(ind_es(1));
        weights(2)=1/(abs(snrss(ind_es(1)-1)-snr_in)+eps);
        alphaa(2)=alpha(ind_es(1)-1);
        weights=weights/sum(weights);
        alpha_smoothing=sum(weights.*alphaa);
    elseif snr_in>snrss(ind_es(1))
        weights(1)=1/(es(1)+eps);
        alphaa(1)=alpha(ind_es(1));
        weights(2)=1/(abs(snrss(ind_es(1)+1)-snr_in)+eps);
        alphaa(2)=alpha(ind_es(1)+1);
        weights=weights/sum(weights);
        alpha_smoothing=sum(weights.*alphaa);
    end
end


speech_order= nmf_V.speech_order;
Gain_min=.1;
Xf=Y(:,1);
if n>1
    a_speech=UserData{1}(1);    a_noise=UserData{1}(2); %shape parameters for the peiors for activatioans
    sm_vector=[a_speech a_noise]; %a_noise for noise and a_speech for speech
    
    Xf=round(Xf);
    nmf_V.X=Xf(:,1);
    Bays_learn_V(nmf_V,sm_vector,nmf_V.V_mean,rand(nmf_V.order,1));
    
    if strcmp(noise_name,'online')
        %-------cash data for online noise estimate--------
        if isempty(local_buffer_cashed)
            local_buffer_cashed=zeros(size(Y,1),size_local_buffer);
        end
        if isempty(newNoiseInBuffer)
            newNoiseInBuffer=0;
        end
        local_buffer_cashed=circshift(local_buffer_cashed,[0 -1]);
        local_buffer_cashed(:,end)=Y(:,1);
        newNoiseInBuffer=newNoiseInBuffer+1;
        if newNoiseInBuffer==size_local_buffer
            newNoiseInBuffer=0;
            %find spectra with minima energy in the local buffer
            en=sum(local_buffer_cashed,1);
            [~,ind]=sort(en,'ascend');
            ind_save=sort(ind(1:update_length_noise_data),'ascend');
            new_noise_data=local_buffer_cashed(:,ind_save);
            
            %update the main buffer with new_noise_data
            noise_data=circshift(noise_data,[0 -update_length_noise_data]);
            noise_data(:,end-update_length_noise_data+1:end)=new_noise_data;
            if n>((size(noise_data,2))/update_length_noise_data)*size_local_buffer %to be sure all the buffers have relevant data
                %the online basis should be obtained from a data that is
                %scaled properly (considering the given spec_scale) to
                %reduce rounding effect
                var_noise_data=mean(noise_data(:).^2);
                noise_data=(noise_data/sqrt(var_noise_data))*spec_scale;
                scale=sqrt(mean(noise_data(:).^2)/var_noise_data);
                
                %update noise dictionary
                updateNoiseBasis(nmf_V,noise_data,5,50,scale);
                nmf_V.X=Xf(:,1);
            end
        end
    end
    V=nmf_V.Ev;
    speech_mask_bays=nmf_V.Sp_mask;
    %update the mean of the priors for the next time frame
    V_mean=nmf_V.V_mean;
    V_mean(1:speech_order)=cal_mean (V_mean(1:speech_order),V(1:speech_order),0);
    V_mean(speech_order+1:end)=cal_mean (V_mean(speech_order+1:end),V(1+speech_order:end,:),alpha_smoothing);
    nmf_V.V_mean=V_mean;
else
    if n==1
        noise_data=init_noise_data;
    end
    %in the first time frame apply regular multiplicative NMF
    init_V=rand(nmf_V.order,1);
    [~,V,cost]=NMF.SNMF_mult({Xf,nmf_V.Et,init_V,2000},...
        {'costfn','kl','threshold',5e-3,'updateBasis',0,'min_it',150});
    speech_mask_bays=nmf_V.Et(:,1:speech_order)*V(1:speech_order)./(nmf_V.Et*V);
    
    nmf_V.X=Xf(:,1);
    nmf_V.V_mean=V+.1;
    
end
%---------------------Filter the noisy signal--------------------------
speech_mask_bays=max(speech_mask_bays,Gain_min);
MagOut=Y.*speech_mask_bays;%
end

function [updated_mean]=cal_mean (V_mean,V,alpha)
updated_mean=alpha*V_mean+(1-alpha)*V+.1;
end

