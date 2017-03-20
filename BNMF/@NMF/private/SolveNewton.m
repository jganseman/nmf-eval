function [a] = SolveNewton(C, ai)
% solve log(a)-psi(a)+1-C to find a using newton-raphson method;

if nargin<2,
    ai = 0.1*ones(size(c));
end;
[M N] = size(ai);
[Mc Nc] = size(C);

if (M==Mc && N == Nc), % No tie
    a = ai+eps;
    cond = 1;

elseif (Mc==1 && Nc>1),  % Tie rows
    cond = 2;

    a = ai(1,:)+eps;

elseif (Mc>1 && Nc==1),  % Tie cols
    cond = 3;
    a = ai(:,1)+eps;

elseif (Mc==1 && Nc==1),  % Tie all
    cond = 4;
    a = ai(1,1)+eps;
end

%solve
flag=1;it=1;
while(flag)

    Delta=(log(a)-psi(0,a)+1-C)./(1./a-psi(1,a));
    IsNeg=1;
    while(sum(sum(IsNeg)))
        an=a-Delta;
        IsNeg=an<0;
        Delta=Delta-(Delta/2).*IsNeg;
    end

    %check conditions
    if(it>10)
        flag=0;
    end
    it=it+1;
    if(sum(sum(an<=0)))
        error('negative number is occured in SolveNewton: perhaps the input C in Inf')
    end
    a=an;
end


switch cond,
    %      case 1, do nothing
    case 2, % tie rows
        a = repmat(a, [M 1]);
    case 3, % tie cols
        a = repmat(a, [1 N]);
    case 4,
        a = a*ones(M, N);
end;
end

