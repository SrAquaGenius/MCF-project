function [res] = payoff(type, op, S, K)
% Payoff: returns the payoff for a given option
% - type: Type of Option: "Am" American or "Eu" European
% - op: Operation type: "Put" or "Call"
% - S: Stock price
% - K: Strike price

res = 0;

if type == "Eu" && op == "Put"
    res = max(K-S, 0);
else if type == "Eu" && op == "Call"
    res = max(S-K, 0);
%how will I define this V_cont.
% we will define equally for American and European option.
end

end
