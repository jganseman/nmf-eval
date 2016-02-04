function [filtout, energy]=generate_dualfilters(filters)

filtout=filters;
%compute LP
  W=filters.psi;
  Phi=filters.phi;
  
  J=size(W);
  L=size(W{1},2);
	if L >1
  [N,M] = size(Phi);
  
  energy=zeros(N,M);
  for j=1:J(2)
    if ~isempty(W{j})
      O=size(W{j});
      for o=1:O(2)
        energy = energy + abs(W{j}{o}).^2;
      end
    end
  end
  energy = energy + abs(Phi).^2;
  
  energym = fliplr(flipud(energy));
  energym = circshift(energym,[1 1]);
  energy = 0.5*(energy+energym);

 Wout = W;
Phiout = Phi;

	for j=1:J(2)
		if ~isempty(W{j})
			for l=1:L
				supp=(abs(W{j}{l})>eps);
				Wout{j}{l} = supp.*(conj(W{j}{l})./(energy));
			end
		end
	end
	supp=(abs(Phi)>eps);
	Phiout=supp.*((Phi)./(energy));

  filtout.dpsi = Wout;
  filtout.dphi = Phiout;


	else 
	%1D case

  [N,M] = size(Phi);
  
  energy=zeros(N,M);
  for j=1:J(2)
        energy = energy + abs(W{j}).^2;
  end
  energy = energy + abs(Phi).^2;
  
  energym = fliplr(flipud(energy));
  energym = circshift(energym,[1 0]);
  energy = 0.5*(energy+energym);

 Wout = W;
Phiout = Phi;

	for j=1:J(2)
		supp=(abs(W{j})>eps);
		Wout{j} = supp.*(conj(W{j})./(energy));
	end
	supp=(abs(Phi)>eps);
	Phiout=supp.*((Phi)./(energy));

  filtout.dpsi = Wout;
  filtout.dphi = Phiout;


end
