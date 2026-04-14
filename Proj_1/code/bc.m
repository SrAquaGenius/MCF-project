function [res] = bc(type, op, side, S, t, K, r, T)
% Boundary Condition: returns the boundary conditions for a given option
% - type: Type of Option: "Am" American or "Eu" European
% - op: Operation type: "Put" or "Call"
% - side: Side boundary: "left" or "right"
% - S: Stock price
% - t: Time value
% - K: Strike price
% - r: Risk free interest rate

res = zeros(size(t));
    
if type == "Eu" && op == "Put" && side == "left"
    res = K*exp(-r*(T-t)) - S;
elseif type == "Eu" && op == "Put" && side == "right"
    res = 0;

elseif type == "Eu" && op == "Call" && side == "left"
    res = 0;
elseif type == "Eu" && op == "Call" && side == "right"
    res = S - K*exp(-r*(T-t));
end

end
