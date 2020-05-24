% Function that applies simple Bayesian belief updating in context-reversal
% gambling task
function [ll, kld, surprise] = Experiment_solution(params,obs,plot)

try
    plot;
catch
    plot = 0;
end

%% Initialise model

% free parameters (will be estimated individually)
% rev_prob = 1-3/60;  % probability of context reversal
% val_cue  = 0.85;    % validity of the cue, needed for likelihood function

rev_prob = params(1);
val_cue  = params(2);

% first, create 'value' of offer based on observed cues:
% v = randi(2,60,1);              % shape;	1=good, 2=bad
% a = randi(2,60,1);              % tone;	1=good, 2=bad

v = obs(:,1);                     % shape;	1=good, 2=bad
a = obs(:,2);                     % tone;	1=good, 2=bad

c = mod(sum([v a],2)',2)' + 1;    % congruency;	1=good, 2=bad

% values = [v a c];
% 
% observed outcomes (wins or losses)
% outcome = zeros(60,1);
% 
% prob = randsample(2,60,true,[0.85 0.15]);
% 
% context = [2*ones(20,1); ones(20,1); 3*ones(20,1)];
% 
% possible_outcomes = [1 2;
%                      2 1]';
% 
% for i=1:length(values)
%     outcome(i) = possible_outcomes(prob(i),values(i,context(i))); % outcome   1=win, 2=loss
% end

outcome = obs(:,3);  % outcome   1=win, 2=loss

% define a simple likelihood function for observing a particular outcome;
% make a matrix with 2 rows (for good and bad offer) and 2 columns (for win and loss)
% entries should be [valcue 1-valcue; 1-valcue valcue]
% val_cue = 0.85;
likelihood = [val_cue 1-val_cue; 
              1-val_cue val_cue];

% initialise prior and posterior
prior     = zeros(length(outcome),3);
posterior = zeros(length(outcome),3);

% rev_prob = 1-3/60;
reversal = [rev_prob      (1-rev_prob)/2 (1-rev_prob)/2; 
           (1-rev_prob)/2 rev_prob       (1-rev_prob)/2; 
           (1-rev_prob)/2 (1-rev_prob)/2 rev_prob];

% initialise regressors for Kullback-Leibler divergence and surprise
kld      = zeros(60,1);
surprise = zeros(60,1);

if plot, figure, end

%% simulate belief updates

for i=1:length(outcome)
   
    if i==1
        prior(i,:) = [1/3 1/3 1/3];
    else
        prior(i,:) = posterior(i-1,:)*reversal;
    end
    
    if plot
    bar(prior(i,:))
    set(gca,'XTickLabel',{'Tone','Shape','Both'}); axis([0,4,0,1]); set(gcf,'color','white'); 
    xlabel('Context','FontSize',12); ylabel('Probability','FontSize',12); title('Belief Updating','FontSize',16);
    
    pause(0.4)
    end
    
    % compute probability of the particular observation in that trial
    prob_obs = [likelihood(v(i),outcome(i)) likelihood(a(i),outcome(i)) likelihood(c(i),outcome(i))];
    
    posterior(i,:) = prior(i,:).*prob_obs/sum(prior(i,:).*prob_obs);
    
    if plot
    bar(prior(i,:))
    set(gca,'XTickLabel',{'Tone','Shape','Both'}); axis([0,4,0,1]); set(gcf,'color','white'); 
    xlabel('Context','FontSize',12); ylabel('Probability','FontSize',12); title('Belief Updating','FontSize',16);
    
    pause(0.4)
    end
    
    kld(i)      = KLD_discrete(posterior(i,:),prior(i,:));  % compute Kullback-Leibler divergence from prior to posterior beliefs
    surprise(i) = -log(sum(prior(i,:).*prob_obs));          % compute (information-theoretic) surprise

end

% check correlation between Kullback-Leibler and surprise!

%% predict behaviour
rule = ~[v-1 a-1 c-1]; % 0 == reject, 1 == accept for the three different contexts {shape,tone,both}

pred_choice = zeros(60,1); % model prediction in terms of probability to accept offer

for i=1:length(outcome)
   pred_choice(i,1) = rule(i,:)*prior(i,:)'; % Bayesian model averaging
end

% Compute nLL und LL
% choice = randi(2,60,1)-1; % 0 == reject, 1 == accept
choice = obs(:,4);

%% compute negative log-likelihood
nll = sum(log(pred_choice.*choice + abs((1-pred_choice)).*abs(~choice)+eps));   
     
prior_rev_prob = [1-3/60 3/60].*100; % prior on reversal probability parameter
      
prior_val_cue  = [0.85 0.15].*100; % prior on cue validity parameter

nll = nll + log(betapdf(rev_prob,prior_rev_prob(1),prior_rev_prob(2)));

nll = nll + log(betapdf(val_cue,prior_val_cue(1),prior_val_cue(2)));

ll = -nll;

