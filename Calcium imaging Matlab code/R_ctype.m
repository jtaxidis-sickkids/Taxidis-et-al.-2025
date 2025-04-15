function R = R_ctype(R,ct,Npy)

if ct == 1                                                           % If plotting PY
    R = R(1:Npy,:,:);                                               % Keep their rates
elseif ct == 2                                                       % Same for IN
    R = R(Npy+1:end,:,:);
end