function Parms =  BSS_EVAL...
  (wav_truth_signal, wav_truth_noise, wav_pred_signal, wav_pred_noise, wav_mix)
% Evaluate performance using BSS Eval 2.0

%% evaluate
if length(wav_pred_noise)==length(wav_truth_noise)
    sep = [wav_pred_noise , wav_pred_signal]';
    orig = [wav_truth_noise , wav_truth_signal]';
else
    minlength=min( length(wav_pred_noise), length(wav_truth_noise) );

    sep = [wav_pred_noise(1:minlength) , wav_pred_signal(1:minlength)]';
    orig = [wav_truth_noise(1:minlength) , wav_truth_signal(1:minlength)]';
end

for i = 1:size( sep, 1)
   [e1,e2,e3] = bss_decomp_gain( sep(i,:), i, orig);
   [sdr(i),sir(i),sar(i)] = bss_crit( e1, e2, e3);
end


[e1,e2,e3] = bss_decomp_gain( wav_mix', 1, wav_truth_signal');
[sdr_,sir_,sar_] = bss_crit( e1, e2, e3);

Parms.SDR=sdr;
Parms.SIR=sir;
Parms.SAR=sar;
Parms.NSDR=Parms.SDR-sdr_;
