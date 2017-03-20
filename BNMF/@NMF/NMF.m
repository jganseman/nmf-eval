classdef NMF<handle
    % NMF class implements algorithms for Bayesian NMF. Both basic and Bayesian NMF are addressed.
    % notation: X=TV where T has 'order' columns and X has W*K dimensions
    % The class allows for training a Bayesian NMF using a mathod based on the Variational
    % Bayes. The class offers different constructors (using function NMF )to just copy another class, combine two NMF classes,
    % or just create a new class. To set the parameters setParameters methods
    % are defined. Additionally, in the source separation or noise
    % reduction when basis are fixed, activations can be learned using the
    % provided methods. Finally, the online learning is made possible using
    % the implemented methods. See demo.m to see how these functions are
    % used.
    % to see the list of properties and methods, type 'NMF' in the command window and then
    % in the last line click on 'Methods'
    % to see this help, type 'help NMF' in the command window
    % NMF is a handle class not a vlaue class.
    %
    %
    % ref:
    % [1] N. Mohammadiha, P. Smaragdis and A. Leijon, “Supervised and
    % Unsupervised Speech Enhancement Using Nonnegative Matrix
    % Factorization,” IEEE Trans. Audio, Speech, and Language Pro-
    % cess., vol. 21, no. 10, pp. 2140–2151, oct. 2013.
    % [2] N. Mohammadiha, J. Taghia, and A. Leijon, “Single channel speech
    % enhancement using Bayesian NMF with recursive temporal updates
    % of prior distributions,” in Proc. IEEE Int. Conf. Acoustics, Speech,
    % and Signal Process. (ICASSP), mar. 2012, pp. 4561–4564.
    % [3] A. T. Cemgil, “Bayesian inference for nonnegative matrix factorisation
    % models,” Computational Intelligence and Neuroscience, vol. 2009,
    % 2009, article ID 785152, 17 pages.
    %
    % Nasser Mohammadiha april 2010
    % Last revise: Nov 2013
    
    properties
        X=[];%input nonnegative matrix to be factorized
        V_mean=[]; %setting of this parameter must be done very carefully,
        %this is actually to set a value for prior parameters
    end
    properties(GetAccess='public',SetAccess='protected')
        order=[];%order of factorization, basis matrix T  has 'order' columns
        speech_order=[];%number of basis of speech signal (or first source) in modeling mixtures
        bound=-inf;%lower bound for variational bayes (VB)
        LogEvidence=[];%loverbound for LogLikelihood during optimization
        Et;% MMSE estimate of T, mean of posterior of T
        Ev;% MMSE estimate of V, mean of posterior of V
        
        
        max_it=200; %maximum number of iterations in VB
        num_initial_it=10; %number of iteration used in the initialization algorithm
        
        Sp_mask=[];%is used for bayesian speech enhancement, sp_com=Sp_mask.*input
        signal_name='';
    end
    properties(GetAccess='public',SetAccess='protected',Hidden=true)%
        At;%prior shape hyperparameters for T. prior for T is Gamma(t,At,Bt./At)
        Bt;%prior scale hyperparameters for T
        Av;%prior shape hyperparameters for V. prior for V is Gamma(v,Av,Bv./Av)
        Bv;%prior scale hyperparameters for V
        
        Alphat;%posterior shape for T
        Betat;%posterior scale for T. q(t)=Gamma(t;Alphat,Betat)
        Alphav;%posterior shape for V
        Betav;%posterior scale for V. q(v)=Gamma(v;Alphav,Betav)
        update_boundFlag=0; %it is time consuming to update bound in each iteration, if set to zero will report just the final result.
        M=[];%Mask for lost data
        AtTying='tie_all'; %Tying type for At. prior for T is Gamma(t,At,Bt./At)
        BtTying='tie_all'; %Tying type for Bt. prior for V is Gamma(v,Av,Bv./Av)
        AvTying='tie_all'; %Tying type for Av. prior for T is Gamma(t,At,Bt./At)
        BvTying='tie_all'; %Tying type for Bv. prior for V is Gamma(v,Av,Bv./Av)
        %                      'free', : Learn parameters for each cell (prone to overfitting)
        %                      'rows', : learn a single parameter for each row
        %                      'cols',: learn a single parameter for each column
        %                      'tie_all': learn a single parameter
        %                      'clamp': Don't learn the hyperparameter
        
        update_hyper_delay=30; %start updating hyperparameters in the BNMF after update_hyper_delay iterations
        update_hyper=1; %update hyperparameters in the BNMF
    end
    properties(Access='protected',Hidden=true)%internal properties
        Lt;%means of logs for T
        Lv;%means of logs for V
        St;%Source sufficient statistics for T
        Sv;%Source sufficient statistics for V
        coefAt=1e-3% initial value At=coefAt*uniform(0,1)
        coefBt=10%initial value Bt=coefBt*uniform(0,1)
        coefAv=1e-3%initial value Av=coefAv*uniform(0,1)
        %coefBv is data dependent
    end
    properties(Access='protected',Hidden=true,Dependent = true)%internal properties
        coefBv;%E(Bv)=coefBv*.5
        W;%X is W*K matrix
        K;
    end
    
    methods(Access='public')
        function nmf=NMF(arg1,arg2,arg3)
            switch nargin
                case 0
                case 1
                    if isa(arg1,'NMF')
                        nmf=arg1;%just copy it
                    else
                        nmf.X=arg1;
                    end
                case 2
                    if isa(arg1,'NMF')
                        nmf.order=arg1.order;
                        nmf.Et=[arg1.Et];
                        nmf.max_it=arg1.max_it; %maximum number of iterations
                        
                        nmf.Alphat=[arg1.Alphat];
                        nmf.Betat=[arg1.Betat];
                        nmf.At=[arg1.At];
                        nmf.Bt=[arg1.Bt];
                        nmf.Alphav=nmf.coefAv*ones(nmf.order,1);%uninformative prior
                        nmf.Betav=nmf.coefBv*ones(nmf.order,1);
                        
                        
                        nmf.update_boundFlag=arg1.update_boundFlag; %it is time consuming to update bound in each iteration,
                        %if set to zero will report just the final result.
                        nmf.Lt=[arg1.Lt];
                        
                        nmf.Ev=arg2;
                    else
                        nmf.X=arg1;
                        nmf.order=arg2;
                    end
                case 3
                    if isa(arg1,'NMF')&& isa(arg2,'NMF')
                        %concatenate two input objects to keep only the information for basis vectors
                        %and give it as the output nmf. Required for source
                        %separation
                        nmf.order=arg1.order+arg2.order;
                        nmf.X=arg3;
                        nmf.speech_order=arg1.order;
                        nmf.Et=[arg1.Et arg2.Et];
                        nmf.max_it=arg1.max_it; %maximum number of iterations
                        
                        nmf.Alphat=[arg1.Alphat arg2.Alphat];
                        nmf.Betat=[arg1.Betat arg2.Betat];
                        nmf.Alphav=nmf.coefAv*ones(nmf.order,1);%uninformative prior
                        nmf.Betav=nmf.coefBv*ones(nmf.order,1);
                        
                        nmf.Ev=nmf.Alphav.*nmf.Betav;
                        
                        nmf.update_boundFlag=arg1.update_boundFlag;
                        nmf.Lt=[arg1.Lt arg2.Lt];
                    end
            end
        end
        %----------------------------------------------------
        function adjust_ShapeparamBasis(nmf,shape_max)
            %set a shape parameter to posterior distributions
            ratio=shape_max/max(nmf.Alphat(:));
            nmf.Alphat=nmf.Alphat*ratio;
            nmf.Betat=nmf.Betat/ratio;
            nmf.Et=nmf.Alphat .* nmf.Betat;
            nmf.Lt=exp(psi(0,nmf.Alphat+eps)).*nmf.Betat+eps;
        end
        %----------------------------------------------------------
        function nmf=setParameters(nmf,varargin) %initialize internal variables for input/output
            nmf.M=ones(size(nmf.X));
            [nmf.AtTying nmf.BtTying nmf.AvTying nmf.BvTying nmf.update_hyper ...
                nmf.update_hyper_delay nmf.max_it At ...
                Bt Av Bv nmf.M nmf.update_boundFlag nmf.num_initial_it nmf.signal_name] = parse_optionlist(nmf,varargin{:});
            
            %initialize hyperparameters
            nmf.At=nmf.coefAt*ones(nmf.W,nmf.order);
            nmf.Bt=nmf.coefBt*ones(nmf.W,nmf.order);
            nmf.Av=nmf.coefAv*ones(nmf.order,nmf.K);
            nmf.Bv=nmf.coefBv*ones(nmf.order,nmf.K);
            
            
            %initialize the posterior means by applying regular NMF
            [nmf.Et,nmf.Ev]=nmf.SNMF_mult({nmf.X,rand(nmf.W,nmf.order),rand(nmf.order,nmf.K),nmf.num_initial_it},...
                {'costfn','kl','threshold',.01,'updateBasis',1});
            
            nmf.Lt=nmf.Et;
            nmf.Alphat=nmf.At;
            nmf.Betat=nmf.Bt;
            nmf.St=zeros(nmf.W,nmf.order);
            
            nmf.Lv=nmf.Ev;
            nmf.Alphav=nmf.Av;
            nmf.Betav=nmf.Bv;
            nmf.Sv=zeros(nmf.order,nmf.K);
        end
        %--------------------------------------------------------------
        function nmf=setParameters_V(nmf,varargin) %initialize internal variables for input/output
            %use this when basis is fixed and only V is supposed to be
            %initialized
            nmf.M=ones(size(nmf.X));
            [nmf.AtTying nmf.BtTying nmf.AvTying nmf.BvTying nmf.update_hyper ...
                nmf.update_hyper_delay nmf.max_it At ...
                Bt Av Bv nmf.M nmf.update_boundFlag nmf.num_initial_it nmf.signal_name] = parse_optionlist(nmf,varargin{:});
            
            nmf.Av=Av;
            nmf.Bv=Bv;
            %initialize the posterior means
            
            [Ett,nmf.Ev]=nmf.SNMF_mult({nmf.X,nmf.Et,nmf.Ev,nmf.num_initial_it},...
                {'costfn','kl','threshold',.01,'updateBasis',0});
            
            nmf.Ev=nmf.Ev+.5; %numerical problems
            nmf.Lv=nmf.Ev;
            nmf.Alphav=nmf.Av;
            nmf.Betav=nmf.Bv;
            nmf.Sv=zeros(nmf.order,nmf.K);
        end
        %---------------------------------------------------------------
        function Bays_learn_V(nmf,sm_vector,Bv,init)
            %learn only V when basis T is fixed. Use the given prior
            %information for temporal modeling
            a_speech=sm_vector(1);
            a_noise=sm_vector(2);
            I_speech=nmf.speech_order;
            
            %set the prior for distribution
            Av=zeros(nmf.order,1);
            Av(1:I_speech)=a_speech;
            Av(I_speech+1:end)=a_noise;
            
            nmf.Ev=init;
            if sm_vector(1)==0
                Av(1:I_speech)=nmf.coefAv*ones(I_speech,1);%uninformative prior
                CC=max(mean(mean(nmf.X)),100);
                Bv(1:I_speech)=CC*ones(I_speech,1);
            end
            if sm_vector(2)==0
                LAv=length(Av(I_speech+1:end));
                Av(I_speech+1:end)=nmf.coefAv*ones(LAv,1);%uninformative prior
                CC=max(mean(mean(nmf.X)),100);
                Bv(I_speech+1:end)=CC*ones(LAv,1);
            end
            
            setParameters_V(nmf,'max_it',200,...
                'update_boundFlag',0,'Av',Av,'Bv',Bv,'num_initial_it',1);
            learn_V(nmf);%the actual VB algorithm
        end
        %-------------------------------------------------------------
        function nmf=train(nmf)
            %learn NMF basis and activations by Variational Bayes
            LogEvidence=zeros(nmf.max_it,1);
            for it=1:nmf.max_it
                nmf=update_statistics(nmf);
                if (nmf.update_boundFlag ~=0 || it==nmf.max_it)
                    [nmf,bound]=update_bound(nmf);%if required and in the last iteration calc the bound
                    LogEvidence(it)=bound; %#ok<AGROW>
                end
                nmf=update_LogMeans(nmf);
                
                if((it>nmf.update_hyper_delay) && (nmf.update_hyper))
                    nmf=update_hyp(nmf);
                end
                nmf.LogEvidence=LogEvidence;
            end
        end
        %-----------------------------------------------------------------
        function updateNoiseBasis(nmf,noise_data,max_it,shapes_basis,scale)
            %online learning of noise basis
            I_speech=nmf.speech_order;
            
            Av=nmf.coefAv*ones(nmf.order,size(noise_data,2));%uninformative prior
            CC=max(mean(mean(noise_data)),100);
            Bv=CC*ones(nmf.order,size(noise_data,2));
            At=nmf.Alphat;
            Bt=nmf.Et;
            
            ratio=shapes_basis/max(max(nmf.Alphat(:,I_speech+1:end)));
            At(:,I_speech+1:end)=nmf.Alphat(:,I_speech+1:end)*ratio;
            Bt(:,I_speech+1:end)=nmf.Et(:,I_speech+1:end);%Gam(a,b/a)
            
            nmf.X=noise_data;
            
            %--------------set parameters--------
            nmf.M=ones(size(nmf.X));
            nmf.At= At;
            nmf.Av= Av;
            nmf.Bt= Bt;
            nmf.Bv= Bv;
            nmf.max_it=max_it;
            
            %initialize the posterior means
            A=nmf.Et;
            A(:,nmf.speech_order+1:end)=rand(size(A(:,nmf.speech_order+1:end)));
            V=rand(nmf.order,nmf.K)+1;
            %initialize the noise basis using regular NMF
            [A(:,nmf.speech_order+1:end),V(nmf.speech_order+1:end,:)]=...
                nmf.SNMF_mult({noise_data,A(:,nmf.speech_order+1:end),V(nmf.speech_order+1:end,:),4},...
                {'costfn','kl','threshold',.01,'updateBasis',1});
            
            %initialize the noise basis+activation of speech and noise using regular NMF
            for l=1:2
                V=(V.*(A'*(nmf.X./(A*V+eps)))./(max(repmat(A'*ones(size(nmf.X,1),1),1,size(nmf.X,2)),1e-3)));
                
                A_n = (A.*((nmf.X./(A*V + eps))*V')./(max(repmat(sum(V,2)',size(nmf.X,1),1),1e-3)));
                A_n = A_n*diag(1./(sum(A_n,1) + eps)); %normalize A such that each basis vector has a norm1=1
                A(:,nmf.speech_order+1:end)=A_n(:,nmf.speech_order+1:end);
            end
            nmf.Ev=V+realmin;
            nmf.Et=A+realmin;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nmf.update_hyper=0;
            nmf.Lt=nmf.Et;
            nmf.Alphat=nmf.At;
            nmf.Betat=nmf.Et./nmf.Alphat;
            nmf.St=zeros(nmf.W,nmf.order);
            nmf.Lv=nmf.Ev;
            nmf.Alphav=nmf.Av;
            nmf.Betav=nmf.Ev./nmf.Alphav;
            nmf.Sv=zeros(nmf.order,nmf.K);
            %--------------------------
            
            for it=1:nmf.max_it
                nmf.St=nmf.Lt.*(((nmf.X.*nmf.M)./(nmf.Lt*nmf.Lv))*nmf.Lv');
                nmf.Sv=nmf.Lv.*(nmf.Lt'*((nmf.X.*nmf.M)./(nmf.Lt*nmf.Lv)));
                
                %use nmf to calculate statistics
                nmf.Alphat(:,nmf.speech_order+1:end)=nmf.At(:,nmf.speech_order+1:end)+nmf.St(:,nmf.speech_order+1:end);
                nmf.Betat(:,nmf.speech_order+1:end)=1./(nmf.At(:,nmf.speech_order+1:end)./(nmf.Bt(:,nmf.speech_order+1:end)+eps) + ...
                    nmf.M*nmf.Ev(nmf.speech_order+1:end,:)'+eps);
                nmf.Et(:,nmf.speech_order+1:end)=nmf.Alphat(:,nmf.speech_order+1:end) .* nmf.Betat(:,nmf.speech_order+1:end);
                nmf.Alphav=nmf.Av+nmf.Sv;
                nmf.Betav=1./(nmf.Av./(nmf.Bv+eps) + nmf.Et'*nmf.M+eps);
                nmf.Ev=nmf.Alphav .* nmf.Betav;
                
                nmf=update_LogMeans_TV(nmf);
                if(nmf.update_hyper)
                    Av=nmf.Av;Av(nmf.speech_order+1:end,2:end)=200;
                    Bv=nmf.Bv;Bv(nmf.speech_order+1:end,2:end)=nmf.Ev(nmf.speech_order+1:end,1:end-1);
                    nmf.Av=Av;nmf.Bv=Bv;
                end
            end
            
            nmf.X=[];
            nmf.Alphav=nmf.Alphav(:,end);
            nmf.Betav=nmf.Betav(:,end);
            
            nmf.Ev=scale*nmf.Alphav.*nmf.Betav;
            nmf.Lt=nmf.Lt;
        end
        %-----------------------------------------------------------------
        function nmf=learn_V(nmf)
            %learn V by VB where posterior distribution of basis T is fixed
            exit=0;
            it=1;
            LogEvidence=zeros(nmf.max_it,1);
            V_b=zeros(nmf.order,1);
            while(~exit)
                nmf=update_statistics_V(nmf);
                nmf=update_LogMeans_V(nmf);
                %convergence
                V_c=sum(abs(nmf.Ev-V_b))/(sum(nmf.Ev))*100;
                if (it>50 && V_c<.5)
                    exit=1;
                elseif(it==nmf.max_it)
                    exit=1;
                end
                V_b=nmf.Ev;
                it=it+1;
            end
            [nmf,bound]=update_bound_V(nmf);%if required and in the last iteration calc the bound
            nmf.LogEvidence=bound;%LogEvidence;
            %calculate the mmse estimator for speech and noise components
            Lt=nmf.Lt+eps;
            Lv=nmf.Lv+eps;
            Lt_speech=nmf.Lt(:,1:nmf.speech_order)+eps;
            Lv_speech=nmf.Lv(1:nmf.speech_order)+eps;
            nmf.Sp_mask=(Lt_speech*Lv_speech)./(Lt*Lv);
        end
        %--------------------------------------------------------------------------
    end
    
    methods(Access='protected',Hidden=true)
        function [varargout] = parse_optionlist(nmf, varargin)
            % parse_optionlist
            % [varargout] = parse_optionlist(defaults, varargin)
            % Usage Example:
            % parse_optionlist(parse_optionlist({'a',1;'b','2';'c',3}, varargin{:});
            % Nasser Mohammadiha march 2010
            
            names =  {'AtTying';'BtTying';'AvTying';'BvTying';'update_hyper';...
                'update_hyper_delay';'max_it';...
                'At';'Bt';'Av';'Bv';'M';'update_boundFlag';'num_initial_it';'signal_name'};
            
            varargout = {nmf.AtTying;nmf.BtTying;nmf.AvTying;nmf.BvTying...
                ;nmf.update_hyper;nmf.update_hyper_delay;nmf.max_it;...
                nmf.At;nmf.Bt;nmf.Av;nmf.Bv;nmf.M; nmf.update_boundFlag; nmf.num_initial_it; nmf.signal_name};
            
            for i=1:2:length(varargin),
                idx = strcmp(varargin{i},names);
                if isempty(idx), error(['Unknown Option : ' varargin{i}]); end;
                varargout{idx} =  varargin{i+1};
            end;
        end
        %-----------------------------------------------------------
        [a] = SolveNewton(C, ai);
        % solve log(a)-psi(a)+1-C to find a using
        % newton-raphson method;
        % the function body is in @NMF/private
        %-----------------------------------------------------------
    end
    methods(Access='private')
        function Nnmf=update_statistics(nmf)
            %update statistics of V and T in VB iterations
            Nnmf=nmf;
            Lt=nmf.Lt+eps;Lv=nmf.Lv+eps;X=nmf.X;M=nmf.M;
            Nnmf.St=Lt.*(((X.*M)./(Lt*Lv))*Lv');
            Nnmf.Sv=Lv.*(Lt'*((X.*M)./(Lt*Lv)));
            
            %use Nnmf to calculate statistics
            At=nmf.At+eps;Bt=nmf.Bt+eps;Ev=nmf.Ev+eps;
            Nnmf.Alphat=At+Nnmf.St;
            Nnmf.Betat=1./(At./Bt + M*Ev');
            Nnmf.Et=Nnmf.Alphat .* Nnmf.Betat;
            Av=nmf.Av+eps;Bv=nmf.Bv+eps;Et=Nnmf.Et+eps;
            Nnmf.Alphav=Av+Nnmf.Sv;
            Nnmf.Betav=1./(Av./Bv + Et'*M);
            Nnmf.Ev=Nnmf.Alphav .* Nnmf.Betav;
            if(sum(sum(isnan(Nnmf.St))) || sum(sum(isnan(Nnmf.Sv))) || ...
                    sum(sum(isnan(Nnmf.Alphat))) ||sum(sum(isnan(Nnmf.Alphav))) ||sum(sum(isnan(Nnmf.Betat))) ||...
                    sum(sum(isnan(Nnmf.Betav))))
                error('NaN is occured in update_statistics!')
            end
        end
        %----------------------------------------------------------
        function Nnmf=update_statistics_V(nmf)
            %update statistics of V in VB iterations
            Nnmf=nmf;
            Lt=nmf.Lt+eps;Lv=nmf.Lv+eps;X=nmf.X;M=nmf.M;
            Nnmf.Sv=Lv.*(Lt'*((X.*M)./(Lt*Lv)));
            Av=nmf.Av+eps;Bv=nmf.Bv+eps;Et=Nnmf.Et+eps;
            Nnmf.Alphav=Av+Nnmf.Sv;
            Nnmf.Betav=1./(Av./Bv + Et'*M);
            Nnmf.Ev=Nnmf.Alphav .* Nnmf.Betav;
            if(sum(sum(isnan(Nnmf.St))) || sum(sum(isnan(Nnmf.Sv))) || ...
                    sum(sum(isnan(Nnmf.Alphat))) ||sum(sum(isnan(Nnmf.Alphav))) ||sum(sum(isnan(Nnmf.Betat))) ||...
                    sum(sum(isnan(Nnmf.Betav))))
                error('NaN is occured in update_statistics!')
            end
        end
        %-----------------------------------------------------------
        function Nnmf=update_LogMeans(nmf)
            %update additional statistics required in VB iterations
            Nnmf=nmf;
            Nnmf.Lt=exp(psi(0,Nnmf.Alphat+eps)).*Nnmf.Betat+eps;
            Nnmf.Lv=exp(psi(0,Nnmf.Alphav+eps)).*Nnmf.Betav+eps;
            if(sum(sum(isnan(Nnmf.Lt))) || sum(sum(isnan(Nnmf.Lv))))
                error('NaN is occured in update_LogMeans!')
            end
        end
        %-----------------------------------------------------------
        function Nnmf=update_LogMeans_V(nmf)
            %update additional statistics required in VB iterations
            Nnmf=nmf;
            Nnmf.Lv=exp(psi(0,Nnmf.Alphav+eps)).*Nnmf.Betav+eps;
            if(sum(sum(isnan(Nnmf.Lt))) || sum(sum(isnan(Nnmf.Lv))))
                error('NaN is occured in update_LogMeans!')
            end
        end
        %--------------------------------------------------------------
        function Nnmf=update_LogMeans_TV(nmf)
            %update statistics required in online learning
            Nnmf=nmf;
            Nnmf.Lt(:,nmf.speech_order+1:end)=exp(psi(0,Nnmf.Alphat(:,nmf.speech_order+1:end)+eps)).*Nnmf.Betat(:,nmf.speech_order+1:end)+eps;
            Nnmf.Lv=exp(psi(0,Nnmf.Alphav+eps)).*Nnmf.Betav+eps;
            if(sum(sum(isnan(Nnmf.Lt))) || sum(sum(isnan(Nnmf.Lv))))
                error('NaN is occured in update_LogMeans!')
            end
        end
        %-----------------------------------------------------------
        function [Nnmf,bound]=update_bound(nmf)
            %update VB bound
            Nnmf=nmf;
            Lt=nmf.Lt;Lv=nmf.Lv;X=nmf.X;At=nmf.At;Av=nmf.Av;Bt=nmf.Bt;Bv=nmf.Bv;
            Et=nmf.Et;Ev=nmf.Ev;Alphat=nmf.Alphat;Alphav=nmf.Alphav;Betat=nmf.Betat;Betav=nmf.Betav;
            
            t1=sum(sum(-Et*Ev-gammaln(X+1)));
            t2=sum(sum(-X .* (((Lt.*log(Lt))*Lv + Lt*(Lv.*log(Lv)))./(Lt*Lv)-log(Lt*Lv))));
            t3=sum(sum(-(At./Bt).*Et-gammaln(At)+At.*log(At./Bt)));
            t4=sum(sum(Alphat.*(log(Betat)+1)+gammaln(Alphat)));
            t5=sum(sum(-(Av./Bv).*Ev-gammaln(Av)+Av.*log(Av./Bv)));
            t6=sum(sum(Alphav.*(log(Betav)+1)+gammaln(Alphav)));
            
            bound=t1+t2+t3+t4+t5+t6;
            Nnmf.bound=bound;
        end
        %-----------------------------------------------------------
        function [Nnmf,bound]=update_bound_V(nmf)
            %update VB bound
            Nnmf=nmf;
            Lt=nmf.Lt;Lv=nmf.Lv;X=nmf.X;At=nmf.At;Av=nmf.Av;Bt=nmf.Bt;Bv=nmf.Bv;
            Et=nmf.Et;Ev=nmf.Ev;Alphat=nmf.Alphat;Alphav=nmf.Alphav;Betat=nmf.Betat;Betav=nmf.Betav;
            if isempty(At)
                At=nmf.coefAt*ones(size(Et));
                Bt=nmf.coefBt*ones(size(Et));
            end
            t1=sum(sum(-Et*Ev-gammaln(X+1)));
            t2=sum(sum(-X .* (((Lt.*log(Lt))*Lv + Lt*(Lv.*log(Lv)))./(Lt*Lv)-log(Lt*Lv))));
            t3=0; t4=0;
            t5=sum(sum(-(Av./Bv).*Ev-gammaln(Av)+Av.*log(Av./Bv)));
            t6=sum(sum(Alphav.*(log(Betav)+1)+gammaln(Alphav)));
            bound=t1+t2+t3+t4+t5+t6;
            Nnmf.bound=bound;
        end
        %-----------------------------------------------------------
        function Nnmf=update_hyp(nmf)
            %update hyperparameters
            Nnmf=nmf;
            At=nmf.At;Et=nmf.Et;Bt=nmf.Bt;Lt=nmf.Lt;
            Zt=Et./Bt-log(Lt./Bt+eps);
            switch nmf.AtTying
                case 'free'
                    Nnmf.At = SolveNewton(Zt, At);
                case 'rows'
                    Nnmf.At = SolveNewton(sum(Zt,1)/nmf.W, At);
                case 'cols'
                    Nnmf.At = SolveNewton(sum(Zt,2)/nmf.order, At);
                case 'tie_all'
                    Nnmf.At = SolveNewton(sum(Zt(:))/(nmf.W*nmf.order), At);
            end
            switch nmf.BtTying
                case 'free'
                    Nnmf.Bt = Et;
                case 'rows'
                    Nnmf.Bt = repmat(sum(At.*Et, 1)./sum(At,1), [nmf.W 1]);
                case 'cols'
                    Nnmf.Bt = repmat(sum(At.*Et, 2)./sum(At,2), [1 nmf.order]);
                case 'tie_all'
                    Nnmf.Bt = sum(sum(At.*Et))./sum(At(:)).*ones(nmf.W, nmf.order);
            end;
            Av=nmf.Av;Ev=nmf.Ev;Bv=nmf.Bv;Lv=nmf.Lv;
            Zv=Ev./Bv-log(Lv./Bv+eps);
            switch nmf.AvTying
                case 'free'
                    Nnmf.Av = SolveNewton(Zv, Av);
                case 'rows'
                    Nnmf.Av = SolveNewton(sum(Zv,1)/nmf.order, Av);
                case 'cols'
                    Nnmf.Av = SolveNewton(sum(Zv,2)/nmf.K, Av);
                case 'tie_all'
                    Nnmf.Av = SolveNewton(sum(Zv(:))/(nmf.order*nmf.K), Av);
            end;
            
            switch nmf.BvTying
                case 'free'
                    Nnmf.Bv = Ev;
                case 'rows'
                    Nnmf.Bv = repmat(sum(Av.*Ev, 1)./sum(Av,1), [nmf.order 1]);
                case 'cols'
                    Nnmf.Bv = repmat(sum(Av.*Ev, 2)./sum(Av,2), [1 nmf.K]);
                case 'tie_all'
                    Nnmf.Bv = sum(sum(Av.*Ev))./sum(Av(:)).*ones(nmf.order, nmf.K);
            end;
            if(sum(sum(isnan(Nnmf.At))) || sum(sum(isnan(Nnmf.Av)))||sum(sum(isnan(Nnmf.Bt))) || sum(sum(isnan(Nnmf.Bv))))
                error('NaN is occured in update_hyper!')
            end
        end
    end
    %set and get methods
    methods
        function nmf=set_property(nmf,feature_name,value)
            Features={'Ev','Et','update_boundFlag','bound'};
            ind=find(strcmp(Features(:),feature_name));
            switch ind
                case 1
                    nmf.Ev=value;
                case 2
                    nmf.Et=value;
                case 3
                    nmf.update_boundFlag=value;
                case 4
                    nmf.bound=value;
            end
        end
        function value=get_property(nmf,feature_name)
            eval(['value=nmf.' feature_name ';'])
        end
        function W=get.W(nmf)%dependent property
            W=size(nmf.X,1);
        end
        function K=get.K(nmf)%dependent property
            K=size(nmf.X,2);
        end
        function coefBv=get.coefBv(nmf)%dependent property
            coefBv=max(mean(mean(nmf.X))/(nmf.coefBt),nmf.coefBt);%
        end
        function nmf=set_Alphat(nmf,value)
            nmf.Alphat=value;
            if ~isempty(nmf.Betat)
                nmf.Et=nmf.Alphat.*nmf.Betat;
            else
                warning('nmf.Et is not valid!,nmf.Betat is empty! ')
            end
        end
        function nmf=set_Betat(nmf,value)
            nmf.Betat=value;
            if ~isempty(nmf.Alphat)
                nmf.Et=nmf.Alphat.*nmf.Betat;
            else
                warning('nmf.Et is not valid! nmf.Alphat is empty!')
            end
        end
        function nmf=set_Alphav(nmf,value)
            nmf.Alphav=value;
            if ~isempty(nmf.Betav)
                nmf.Ev=nmf.Alphav.*nmf.Betav;
            else
                warning('nmf.Ev is not valid! nmf.Betav is empty!')
            end
        end
        function nmf=set_Betav(nmf,value)
            nmf.Betav=value;
            if ~isempty(nmf.Alphav)
                nmf.Ev=nmf.Alphav.*nmf.Betav;
            else
                warning('nmf.Ev is not valid! nmf.Alphav is empty!')
            end
        end
        
    end
    methods(Access='public',Static=true)
        %-----------------------------------------------------------
        [A,X,cost] = SNMF_mult(varargin,Dataoptions)
        % multiplicative rules for KL-NMF
        %function body is located at @NMF
    end
    
end
