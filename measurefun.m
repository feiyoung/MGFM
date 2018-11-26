function [measure] = measurefun(H, hH, type)
if(~exist('type', 'var'))
    type = 'canonical';
end
q = size(H, 2);
switch type 
    case 'canonical'
        [~, ~, r] = canoncorr(hH,H);
        measure = r(q);
    case 'Fnorm'
        measure = norm(H- hH, 'fro');
end