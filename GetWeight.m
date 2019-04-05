function [W, Wfuse, Wwing, Wpay] = GetWeight(AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

Wfuse = GetWfuse(AR, S);

Wwing = GetWingWeight(AR, S);

Wpay = GetWpay(AR, S);

W = Wfuse + Wwing + Wpay;







