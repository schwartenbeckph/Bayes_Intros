function [KL] = KLD_discrete(P,Q)
%
%Calculate KL divergence KL(P||Q) between two discrete probability
%distributions; P = posterior, Q = prior

KL = sum(P.*(log(P./Q)));