function curvedata = eval_sk_curve_function(data,useStandardRes,input_R,isSimple)

if ~exist('useStandardRes','var')
    useStandardRes = true;
end

if ~exist('input_R','var') | isempty(input_R)
    if exist('useStandardRes','var')
        if useStandardRes
            load('B.mat','B');
            R = B;
            R(1) = 0.00001;
        else
            R = 0.00001:0.0001:1;
        end
    else
        R = 0.00001:0.0001:0.99999;
    end
else
    R = input_R;
end


gammaexp = data.k;
gammaexp = flipud(gammaexp(:)); % since gammaexp is in a decreasing order

polyExp = 0;
for i = 1:length(gammaexp)
    polyExp = polyExp + (R.^(i-1))*gammaexp(i);
end

polyExp_I = polyExp == 0;
polyExp(polyExp_I) = 0.0001;

r = R.^(1./polyExp);

I = R < data.b;
r(I) = data.a*R(I);



if exist('isSimple','var') & isSimple
    Q = 0;
    return
end


% -------------

order = length(gammaexp) - 1;

switch order
    case 0
        a = gammaexp(1);
        Q = a*ones(1,length(r));
        
    case 1

        lR = log(R);
        R2 = R.^2;
        R3 = R.^3;
        a1 = gammaexp(2);
        a0 = gammaexp(1);

        den = ((a1^3*lR - 2*a1^3).*R3 + ...
            (a1^2 + a1^2*lR.^2 - 4*a1^2*a0 - 2*a1^2*lR).*R2 + ...
            (-a0^2*a1*lR + 2*a1*a0 - 2*a1*a0^2 - 2*a1*lR*a0).*R ...
            + a0^2);
        den_I = den == 0;
        den(den_I) = 0.0001;
        
        Q = - (a1*R + a0).^2.*((a1*lR - a1).*R-a0)./den;

    case 2

        % gammaexp was already flipped
        a2 = gammaexp(3);
        a1 = gammaexp(2);
        a0 = gammaexp(1);


        lR = log(R);
        lR2 = lR.^2;
        R2 = R.^2;
        R3 = R.^3;
        R4 = R.^4;
        R5 = R.^5;
        R6 = R.^6;


        den = (4*lR*a2^3 - 4*a2^3).*R6 + (-10*a2^2*a1 + 7*lR*a2^2*a1).*R5 + ...
            (4*lR*a2*a1^2 - 8*a2^2*a0 + a2^2 + 4*lR2*a2^2 - 8*a2*a1^2 - 4*lR*a2^2).*R4 + ...
            (-2*a1^3 - 6*lR*a2*a1 + 4*lR2*a2*a1 + lR*a1^3 - 12*a2*a1*a0 + 2*a2*a1 - 2*lR*a2*a1*a0).*R3 + ...
            (-4*lR*a2*a0^2 - 4*a1^2*a0 - 4*a2*a0^2 + lR2*a1^2 - 4*lR*a2*a0+2*a2*a0 + a1^2 - 2*lR*a1^2).*R2 + ...
            (-a0^2*lR*a1 + 2*a1*a0 - 2*lR*a1*a0 - 2*a1*a0^2).*R + a0^2;

        num1 = a2*R2 + a1*R + a0;

        num2 = (2*lR*a2 - a2).*R2 + (lR*a1-a1).*R - a0;

        Q = -(num1.^2).*num2./den;

    otherwise

        fprintf('order %d not supported\n',order);
        Q = 0;
        
end

I = R < data.b;
Q(I) = 1;
% Q = smooth(Q,'moving',5);
% Q = Q(:)';

curvedata.R = R;
curvedata.r = r;
curvedata.Q = Q;


% Qc = sqrt(3)/(sqrt(3)-1)*(1-sqrt(1./(2*Q+1)));
