% Script for MLE or MAP parameter estimation

clear
close all

load('behav_results')

subj = {'s02','s03','s04','s05','s07','s08','s09','s10','s11','s12','s13','s14','s15'};

% free parameters (will be estimated individually)
rev_prob = 1-3/60;  % probability of context reversal
val_cue  = 0.85;    % validity of the cue

initparams = [rev_prob val_cue];

LB = [0 0];
UB = [1 1];

options = optimset('Algorithm','interior-point','TolFun',1e-6,'MaxIter',1e5); 

for i=1%:length(subj);
    
    if strcmp(subj{i},'s02') || strcmp(subj{i},'s05')
       sess = {'s1','s3'};
    elseif strcmp(subj{i},'s11')
       sess = {'s1','s2'};
    else
       sess = {'s1','s2','s3'};
    end
    
    for ii=2%:length(sess)

        % simulate trials with visual and auditory cues
        % v = randi(2,60,1); % visual cue;    1=good, 2=bad
        % a = randi(2,60,1); % auditive cue;  1=good, 2=bad
        v = behav.(subj{i}).(sess{ii}).cuevis; v(v==0) = 2; % visual cue;    1=good, 2=bad
        a = behav.(subj{i}).(sess{ii}).cueaud; a(a==0) = 2; % auditive cue;  1=good, 2=bad

        % simulate outcomes
        % outcome = zeros(60,1);
        % prob = randsample(2,60,true,[0.85 0.15]);
        % 
        % possible_outcomes = [1 2;
        %                      2 1]';
        % 
        % for i=1:length(values)
        %     outcome(i) = possible_outcomes(prob(i),values(i,context(i))); % outcome   1=win, 2=loss
        % end
        outcome = behav.(subj{i}).(sess{ii}).outcome; outcome(outcome==0) = 2;  % outcome   1=win, 2=loss

        % simulate choices
        % choice = randi(2,60,1)-1; % 0 == reject, 1 == accept
        choice = behav.(subj{i}).(sess{ii}).response;

        obs = [v a outcome choice];

        % invert model to extimate individual parameters
        params = fmincon(@Experiment_solution,initparams,[],[],[],[],LB,UB,[],options, obs);   

        % based on individual parameters, simulate behaviour of surprise
        % and KL-divergence
        [ll, kld, surprise] = Experiment_solution(params,obs,1);

    end

end

